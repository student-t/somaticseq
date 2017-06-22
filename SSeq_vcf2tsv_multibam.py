#!/usr/bin/env python3

import sys, argparse, math, gzip, os, pysam
import regex as re
import scipy.stats as stats
import genomic_file_handlers as genome
import pileup_reader as pileup
from read_info_extractor import * 
from copy import copy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

input_sites = parser.add_mutually_exclusive_group()
input_sites.add_argument('-myvcf',  '--vcf-format',           type=str,   help='Input file is VCF formatted.', required=False, default=None)
input_sites.add_argument('-mybed',  '--bed-format',           type=str,   help='Input file is BED formatted.', required=False, default=None)
input_sites.add_argument('-mypos',  '--positions-list',       type=str,   help='A list of positions: tab seperating contig and positions.', required=False, default=None)

parser.add_argument('-nprefix', '--tumor-prefixes',   narg='*', type=str,   help='normal prefixes',  required=True, default=None)
parser.add_argument('-tprefix', '--normal-prefixes',  narg='*', type=str,   help='tumor prefixes',   required=True, default=None)
parser.add_argument('-nbams',   '--normal-bam-files', narg='*', type=str,   help='Normal BAM Files', required=True, default=None)
parser.add_argument('-tbams',   '--tumor-bam-files',  narg='*', type=str,   help='Tumor BAM Files',  required=True, default=None)

parser.add_argument('-truth',     '--ground-truth-vcf',       type=str,   help='VCF of true hits',  required=False, default=None)
parser.add_argument('-dbsnp',     '--dbsnp-vcf',              type=str,   help='dbSNP VCF: do not use if input VCF is annotated', required=False, default=None)
parser.add_argument('-cosmic',    '--cosmic-vcf',             type=str,   help='COSMIC VCF: do not use if input VCF is annotated',   required=False, default=None)

parser.add_argument('-ref',     '--genome-reference',         type=str,   help='.fasta.fai file to get the contigs', required=True, default=None)
parser.add_argument('-dedup',   '--deduplicate',     action='store_true', help='Do not consider duplicate reads from BAM files. Default is to count everything', required=False, default=False)

parser.add_argument('-minMQ',     '--minimum-mapping-quality',type=float, help='Minimum mapping quality below which is considered poor', required=False, default=1)
parser.add_argument('-minBQ',     '--minimum-base-quality',   type=float, help='Minimum base quality below which is considered poor', required=False, default=5)
parser.add_argument('-mincaller', '--minimum-num-callers',    type=float, help='Minimum number of tools to be considered', required=False, default=0)

parser.add_argument('-scale',      '--p-scale',               type=str,   help='phred, fraction, or none', required=False, default=None)
parser.add_argument('-outfile',    '--output-tsv-file',       type=str,   help='Output TSV Name', required=False, default=os.sys.stdout)

args = parser.parse_args()


# Rename input:
is_vcf    = args.vcf_format
is_bed    = args.bed_format
is_pos    = args.positions_list

nbam_fn   = args.normal_bam_files
tbam_fn   = args.tumor_bam_files
n_prefix  = args.normal_prefixes
t_prefix  = args.tumor_prefixes

truth     = args.ground_truth_vcf
cosmic    = args.cosmic_vcf
dbsnp     = args.dbsnp_vcf

min_mq    = args.minimum_mapping_quality
min_bq    = args.minimum_base_quality
ref_fa    = args.genome_reference
p_scale   = args.p_scale

outfile   = args.output_tsv_file

# Convert contig_sequence to chrom_seq dict:
fai_file  = ref_fa + '.fai'
chrom_seq = genome.faiordict2contigorder(fai_file, 'fai')

assert len(nbam_fn) == len(tbam_fn) == len(n_prefix) == len(t_prefix)

paired_prefixes = tuple( zip(n_prefix, t_prefix) )
paired_bams = tuple( zip(nbam_fn, tbam_fn) )

# Determine input format:
if is_vcf:
    mysites = is_vcf
elif is_bed:
    mysites = is_bed
elif is_pos:
    mysites = is_pos
else:
    mysites = fai_file
    print('No position supplied. Will evaluate the whole genome.', file=sys.stderr)

# Re-scale output or not:
if p_scale == None:
    print('NO RE-SCALING', file=sys.stderr)
elif p_scale.lower() == 'phred':
    p_scale = 'phred'
elif p_scale.lower() == 'fraction':
    p_scale = 'fraction'
else:
    p_scale = None
    print('NO RE-SCALING', file=sys.stderr)

def rescale(x, original=None, rescale_to=p_scale, max_phred=1001):
    if ( rescale_to == None ) or ( original.lower() == rescale_to.lower() ):
        y = x if isinstance(x, int) else '%.2f' % x
    elif original.lower() == 'fraction' and rescale_to == 'phred':
        y = genome.p2phred(x, max_phred=max_phred)
        y = '%.2f' % y
    elif original.lower() == 'phred' and rescale_to == 'fraction':
        y = genome.phred2p(x)
        y = '%.2f' % y
    return y
    

# Define NaN and Inf:
nan = float('nan')
inf = float('inf')
pattern_chr_position = genome.pattern_chr_position


# Header for the output data, created here so I won't have to indent this line:
out_header = \
'CHROM\t\
POS\t\
ID\t\
REF\t\
ALT\t\
if_dbsnp\t\
COMMON\t\
if_COSMIC\t\
COSMIC_CNT\t\
MaxHomopolymer_Length\t\
SiteHomopolymer_Length\t\
InDel_Length\t\
TrueVariant_or_False'

bam_headers = '{prefix}_DP\t\
{prefix}_REF_MQ\t\
{prefix}_ALT_MQ\t\
{prefix}_Z_Ranksums_MQ\t\
{prefix}_REF_BQ\t\
{prefix}_ALT_BQ\t\
{prefix}_Z_Ranksums_BQ\t\
{prefix}_REF_NM\t\
{prefix}_ALT_NM\t\
{prefix}_NM_Diff\t\
{prefix}_REF_Concordant\t\
{prefix}_REF_Discordant\t\
{prefix}_ALT_Concordant\t\
{prefix}_ALT_Discordant\t\
{prefix}_Concordance_FET\t\
{prefix}_REF_FOR\t\
{prefix}_REF_REV\t\
{prefix}_ALT_FOR\t\
{prefix}_ALT_REV\t\
{prefix}_StrandBias_FET\t\
{prefix}_Z_Ranksums_EndPos\t\
{prefix}_REF_Clipped_Reads\t\
{prefix}_ALT_Clipped_Reads\t\
{prefix}_Clipping_FET\t\
{prefix}_MQ0\t\
{prefix}_Other_Reads\t\
{prefix}_Poor_Reads\t\
{prefix}_REF_InDel_3bp\t\
{prefix}_REF_InDel_2bp\t\
{prefix}_REF_InDel_1bp\t\
{prefix}_ALT_InDel_3bp\t\
{prefix}_ALT_InDel_2bp\t\
{prefix}_ALT_InDel_1bp'


## Running
with genome.open_textfile(mysites) as my_sites, open(outfile, 'w') as outhandle:
        
    my_line = my_sites.readline().rstrip()
    
    nbam    = pysam.AlignmentFile(nbam_fn)
    tbam    = pysam.AlignmentFile(tbam_fn)
    ref_fa  = pysam.FastaFile(ref_fa)
    
    if truth:
        truth = genome.open_textfile(truth)
        truth_line = truth.readline().rstrip()
        while truth_line.startswith('#'):
            truth_line = truth.readline().rstrip()
    
    if cosmic:
        cosmic = genome.open_textfile(cosmic)
        cosmic_line = cosmic.readline().rstrip()
        while cosmic_line.startswith('#'):
            cosmic_line = cosmic.readline().rstrip()

    if dbsnp:
        dbsnp = genome.open_textfile(dbsnp)
        dbsnp_line = dbsnp.readline().rstrip()
        while dbsnp_line.startswith('#'):
            dbsnp_line = dbsnp.readline().rstrip()
    
    
    # Get through all the headers:
    while my_line.startswith('#') or my_line.startswith('track='):
        my_line = my_sites.readline().rstrip()
        
    for nbam_i, tbam_i in paired_prefixes:
        out_header = out_header + '\t' + bam_headers.format(prefix=nbam_i) + '\t' + bam_headers.format(prefix=tbam_i) + '\t' + '{}_SOR'.format( '{}_{}'.format(nbam_i, tbam_i))
    
    # Write out the header
    outhandle.write( out_header  + '\n' )
    
    # Make things '{Feature_1}\t{Feature_2}\t....{Feature_N}'
    features = out_header.split('\t')
    out_string = ''
    for feature_i in features:
        out_string = out_string + '\t{' + feature_i + '}'
    out_string = out_string.lstrip('\t')


    while my_line:
        
        # If VCF, get all the variants with the same coordinate into a list:
        if is_vcf:
            
            my_vcf = genome.Vcf_line( my_line )
            my_coordinates = [(my_vcf.chromosome, my_vcf.position)]
            
            variants_at_my_coordinate = []
            
            alt_bases = my_vcf.altbase.split(',')
            for alt_i in alt_bases:
                vcf_i = copy(my_vcf)
                vcf_i.altbase = alt_i
                variants_at_my_coordinate.append( vcf_i )

            
            # As long as the "coordinate" stays the same, it will keep reading until it's different.
            while my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):

                my_line = my_sites.readline().rstrip()
                my_vcf = genome.Vcf_line( my_line )

                if my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):
                    
                    alt_bases = my_vcf.altbase.split(',')
                    for alt_i in alt_bases:
                        
                        vcf_i = copy(my_vcf)
                        vcf_i.altbase = alt_i
                        variants_at_my_coordinate.append( vcf_i )     
        
        elif is_bed:
            bed_item = my_line.split('\t')
            my_coordinates = genomic_coordinates( bed_item[0], int(bed_item[1])+1, int(bed_item[2]) )
            
        elif is_pos:
            pos_item = my_line.split('\t')
            my_coordinates = genomic_coordinates( pos_item[0], int(pos_item[1]), int(pos_item[1]) )
            
        elif fai_file:
            fai_item = my_line.split('\t')
            my_coordinates = genomic_coordinates( fai_item[0], 1, int(fai_item[1]) )
        
        ##### ##### ##### ##### ##### #####
        for my_coordinate in my_coordinates:
            
            ######## If VCF, can get ref base, variant base, as well as other identifying information ######## 
            if is_vcf:
                
                ref_bases = []
                alt_bases = []
                indel_lengths = []
                all_my_identifiers = []
                
                for variant_i in variants_at_my_coordinate:

                    ref_base = variant_i.refbase
                    first_alt = variant_i.altbase.split(',')[0]
                    indel_length = len(first_alt) - len(ref_base)

                    ref_bases.append( ref_base )
                    alt_bases.append( first_alt )
                    indel_lengths.append( indel_length )
                    
                    # Extract these information if they exist in the VCF file, but they could be re-written if dbSNP/COSMIC are supplied.
                    if_dbsnp  = 1 if re.search(r'rs[0-9]+', variant_i.identifier) else 0
                    if_cosmic = 1 if re.search(r'COS[MN][0-9]+', variant_i.identifier) else 0
                    if_common = 1 if variant_i.get_info_value('COMMON') == '1' else 0
                    num_cases = variant_i.get_info_value('CNT') if variant_i.get_info_value('CNT') else nan
                    
                    if variant_i.identifier == '.':
                        my_identifier_i = set()
                    else:
                        my_identifier_i = variant_i.identifier.split(';')
                        my_identifier_i = set( my_identifier_i )
                    
                    all_my_identifiers.append( my_identifier_i )
                                
            ## If not, 1) get ref_base, first_alt from other VCF files. 
            #          2) Create placeholders for dbSNP and COSMIC that can be overwritten with dbSNP/COSMIC VCF files (if provided)
            else:
                variants_at_my_coordinate = [None] # Just to have something to iterate
                ref_base = first_alt = indel_length = None
                
                # Could be re-written if dbSNP/COSMIC are supplied. If not, they will remain NaN.
                if_dbsnp = if_cosmic = if_common = num_cases = nan

            # Keep track of NumCallers:
            num_callers = 0
            
            #################################### Find the same coordinate in those VCF files ####################################
            if args.ground_truth_vcf: got_truth,   truth_variants,   truth_line   = genome.find_vcf_at_coordinate(my_coordinate, truth_line,   truth,   chrom_seq)
            if args.dbsnp_vcf:        got_dbsnp,   dbsnp_variants,   dbsnp_line   = genome.find_vcf_at_coordinate(my_coordinate, dbsnp_line,   dbsnp,   chrom_seq)
            if args.cosmic_vcf:       got_cosmic,  cosmic_variants,  cosmic_line  = genome.find_vcf_at_coordinate(my_coordinate, cosmic_line,  cosmic,  chrom_seq)
            
            
            # Now, use pysam to look into the BAM file(s), variant by variant from the input:
            for ith_call, my_call in enumerate( variants_at_my_coordinate ):
                
                if is_vcf:
                    # The particular line in the input VCF file:
                    variant_id = ( (my_call.chromosome, my_call.position), my_call.refbase, my_call.altbase )

                    ref_base       = ref_bases[ith_call]
                    first_alt      = alt_bases[ith_call]
                    indel_length   = indel_lengths[ith_call]
                    my_identifiers = all_my_identifiers[ith_call]
                    
                else:
                    variant_id = ( (my_coordinate[0], my_coordinate[1]), ref_base, first_alt )


                ########## Ground truth file ##########
                if args.ground_truth_vcf:
                    if variant_id in truth_variants.keys():
                        judgement = 1
                        my_identifiers.add('TruePositive')
                    else:
                        judgement = 0
                        my_identifiers.add('FalsePositive')
                else:
                    judgement = nan


                ########## dbSNP ########## Will overwrite dbSNP info from input VCF file
                if args.dbsnp_vcf:
                    if variant_id in dbsnp_variants.keys():

                        dbsnp_variant_i = dbsnp_variants[variant_id]
                        if_dbsnp = 1
                        if_common = 1 if dbsnp_variant_i.get_info_value('COMMON') == '1' else 0

                        rsID = dbsnp_variant_i.identifier.split(',')
                        for ID_i in rsID:
                            my_identifiers.add( ID_i )

                    else:
                        if_dbsnp = if_common = 0

                
                ########## COSMIC ########## Will overwrite COSMIC info from input VCF file
                if args.cosmic_vcf:
                    if variant_id in cosmic_variants.keys():

                        cosmic_variant_i = cosmic_variants[variant_id]
                        
                        # If designated as SNP, make it "non-cosmic" and make CNT=nan.
                        if cosmic_variant_i.get_info_value('SNP'):
                            if_cosmic = 0
                            num_cases = nan
                        
                        else:
                            if_cosmic = 1
                            num_cases = cosmic_variant_i.get_info_value('CNT')
                            num_cases = num_cases if num_cases else nan

                        # COSMIC ID still intact:
                        cosmicID = cosmic_variant_i.identifier.split(',')
                        for ID_i in cosmicID:
                            my_identifiers.add( ID_i )
                            
                    else:
                        if_cosmic = num_cases = 0
                
                
                ########################################################################################
                # Tumor BAM file:
                t_reads = tbam.fetch( my_coordinate[0], my_coordinate[1]-1, my_coordinate[1] )
                
                t_ref_read_mq = []
                t_alt_read_mq = []
                t_ref_read_bq = []
                t_alt_read_bq = []
                t_ref_edit_distance = []
                t_alt_edit_distance = []
                
                t_ref_concordant_reads = t_alt_concordant_reads = t_ref_discordant_reads = t_alt_discordant_reads = 0
                t_ref_for = t_ref_rev = t_alt_for = t_alt_rev = T_dp = 0
                t_ref_SC_reads = t_alt_SC_reads = t_ref_notSC_reads = t_alt_notSC_reads = 0
                t_MQ0 = 0
                
                t_ref_pos_from_end = []
                t_alt_pos_from_end = []
                t_ref_flanking_indel = []
                t_alt_flanking_indel = []
                
                t_noise_read_count = t_poor_read_count  = 0
                
                for read_i in t_reads:
                    if not read_i.is_unmapped and dedup_test(read_i):
                        
                        T_dp += 1
                        
                        code_i, ith_base, base_call_i, indel_length_i, flanking_indel_i = position_of_aligned_read(read_i, my_coordinate[1]-1 )
                        
                        if read_i.mapping_quality < min_mq and mean(read_i.query_qualities) < min_bq:
                            t_poor_read_count += 1
                        
                        if read_i.mapping_quality == 0:
                            t_MQ0 += 1
                        
                        # Reference calls:
                        if code_i == 1 and base_call_i == ref_base[0]:
                        
                            t_ref_read_mq.append( read_i.mapping_quality )
                            t_ref_read_bq.append( read_i.query_qualities[ith_base] )
                            
                            try:
                                t_ref_edit_distance.append( read_i.get_tag('NM') )
                            except KeyError:
                                pass
                            
                            # Concordance
                            if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                t_ref_concordant_reads += 1
                            elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                t_ref_discordant_reads += 1
                            
                            # Orientation
                            if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                t_ref_for += 1
                            elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                t_ref_rev += 1
                            
                            # Soft-clipped reads?
                            if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                                t_ref_SC_reads += 1
                            else:
                                t_ref_notSC_reads += 1

                            # Distance from the end of the read:
                            if ith_base != None:
                                t_ref_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                
                            # Flanking indels:
                            t_ref_flanking_indel.append( flanking_indel_i )

                        
                        # Alternate calls:
                        # SNV, or Deletion, or Insertion where I do not check for matching indel length
                        elif (indel_length == 0 and code_i == 1 and base_call_i == first_alt) or \
                             (indel_length < 0  and code_i == 2 and indel_length == indel_length_i) or \
                             (indel_length > 0  and code_i == 3):
                            
                            t_alt_read_mq.append( read_i.mapping_quality )
                            t_alt_read_bq.append( read_i.query_qualities[ith_base] )
                            
                            try:
                                t_alt_edit_distance.append( read_i.get_tag('NM') )
                            except KeyError:
                                pass
                            
                            # Concordance
                            if        read_i.is_proper_pair  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                t_alt_concordant_reads += 1
                            elif (not read_i.is_proper_pair) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                t_alt_discordant_reads += 1
                            
                            # Orientation
                            if (not read_i.is_reverse) and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                t_alt_for += 1
                            elif    read_i.is_reverse  and read_i.mapping_quality >= min_mq and read_i.query_qualities[ith_base] >= min_bq:
                                t_alt_rev += 1
                            
                            # Soft-clipped reads?
                            if read_i.cigar[0][0] == cigar_soft_clip or read_i.cigar[-1][0] == cigar_soft_clip:
                                t_alt_SC_reads += 1
                            else:
                                t_alt_notSC_reads += 1

                            # Distance from the end of the read:
                            if ith_base != None:
                                t_alt_pos_from_end.append( min(ith_base, read_i.query_length-ith_base) )
                                                    
                            # Flanking indels:
                            t_alt_flanking_indel.append( flanking_indel_i )
                        
                        
                        # Inconsistent read or 2nd alternate calls:
                        else:
                            t_noise_read_count += 1
                                
                
                # Done extracting info from tumor BAM. Now tally them:
                t_ref_mq        = mean(t_ref_read_mq)
                t_alt_mq        = mean(t_alt_read_mq)
                t_z_ranksums_mq = stats.ranksums(t_alt_read_mq, t_ref_read_mq)[0]
                
                t_ref_bq        = mean(t_ref_read_bq)
                t_alt_bq        = mean(t_alt_read_bq)
                t_z_ranksums_bq = stats.ranksums(t_alt_read_bq, t_ref_read_bq)[0]
                
                t_ref_NM        = mean(t_ref_edit_distance)
                t_alt_NM        = mean(t_alt_edit_distance)
                t_z_ranksums_NM = stats.ranksums(t_alt_edit_distance, t_ref_edit_distance)[0]
                t_NM_Diff       = t_alt_NM - t_ref_NM - abs(indel_length)
                
                t_concordance_fet = stats.fisher_exact(( (t_ref_concordant_reads, t_alt_concordant_reads), (t_ref_discordant_reads, t_alt_discordant_reads) ))[1]
                t_strandbias_fet  = stats.fisher_exact(( (t_ref_for, t_alt_for), (t_ref_rev, t_alt_rev) ))[1]
                t_clipping_fet    = stats.fisher_exact(( (t_ref_notSC_reads, t_alt_notSC_reads), (t_ref_SC_reads, t_alt_SC_reads) ))[1]
                
                t_z_ranksums_endpos = stats.ranksums(t_alt_pos_from_end, t_ref_pos_from_end)[0]
                
                t_ref_indel_1bp = t_ref_flanking_indel.count(1)
                t_ref_indel_2bp = t_ref_flanking_indel.count(2) + t_ref_indel_1bp
                t_ref_indel_3bp = t_ref_flanking_indel.count(3) + t_ref_indel_2bp + t_ref_indel_1bp
                t_alt_indel_1bp = t_alt_flanking_indel.count(1)
                t_alt_indel_2bp = t_alt_flanking_indel.count(2) + t_alt_indel_1bp
                t_alt_indel_3bp = t_alt_flanking_indel.count(3) + t_alt_indel_2bp + t_alt_indel_1bp
    
                
                # Odds Ratio (just like VarDict, but get from BAM)
                sor_numerator   = (n_alt_for + n_alt_rev) * (t_ref_for + t_ref_rev)
                sor_denominator = (n_ref_for + n_ref_rev) * (t_alt_for + t_alt_rev)
                if sor_numerator == 0 and sor_denominator == 0:
                    sor = nan
                elif sor_denominator == 0:
                    sor = 100
                else:
                    sor = sor_numerator / sor_denominator
                    if sor >= 100:
                        sor = 100
                
                ############################################################################################
                ############################################################################################
                # Homopolymer eval (Make sure to modify for INDEL):
                # The min and max is to prevent the +/- 20 bases from exceeding the reference sequence
                lseq  = ref_fa.fetch(my_coordinate[0], max(0, my_coordinate[1]-20), my_coordinate[1])
                rseq  = ref_fa.fetch(my_coordinate[0], my_coordinate[1]+1, min(ref_fa.get_reference_length(my_coordinate[0])+1, my_coordinate[1]+21) )
                
                # This is to get around buy in old version of pysam that reads the reference sequence in bytes instead of strings
                lseq = lseq.decode() if isinstance(lseq, bytes) else lseq
                rseq = rseq.decode() if isinstance(rseq, bytes) else rseq
                
                seq41_ref = lseq + ref_base  + rseq
                seq41_alt = lseq + first_alt + rseq
                
                ref_counts = genome.count_repeating_bases(seq41_ref)
                alt_counts = genome.count_repeating_bases(seq41_alt)
                
                homopolymer_length = max( max(ref_counts), max(alt_counts) )
                
                # Homopolymer spanning the variant site:
                ref_c = 0
                alt_c = 0
                for i in rseq:
                    if i == ref_base:
                        ref_c += 1
                    else:
                        break
                        
                for i in lseq[::-1]:
                    if i == ref_base:
                        ref_c += 1
                    else:
                        break
                
                for i in rseq:
                    if i == first_alt:
                        alt_c += 1
                    else:
                        break
                        
                for i in lseq[::-1]:
                    if i == first_alt:
                        alt_c += 1
                    else:
                        break
    
                site_homopolymer_length = max( alt_c+1, ref_c+1 )
    
                if my_identifiers:
                    my_identifiers = ';'.join(my_identifiers)
                else:
                    my_identifiers = '.'
                    
                ###
                out_line = out_header.format( \
                CHROM                   = my_coordinate[0],                                       \
                POS                     = my_coordinate[1],                                       \
                ID                      = my_identifiers,                                         \
                REF                     = ref_base,                                               \
                ALT                     = first_alt,                                              \
                if_MuTect               = mutect_classification,                                  \
                if_Strelka              = strelka_classification,                                 \
                if_VarScan2             = varscan_classification,                                 \
                if_JointSNVMix2         = jointsnvmix2_classification,                            \
                if_SomaticSniper        = sniper_classification,                                  \
                if_VarDict              = vardict_classification,                                 \
                if_LoFreq               = lofreq_classification,                                  \
                if_Scalpel              = scalpel_classification,                                 \
                MuSE_Tier               = muse_classification,                                    \
                Strelka_Score           = somatic_evs,                                            \
                Strelka_QSS             = qss,                                                    \
                Strelka_TQSS            = tqss,                                                   \
                VarScan2_Score          = rescale(score_varscan2,      'phred', p_scale, 1001),   \
                SNVMix2_Score           = rescale(score_jointsnvmix2,  'phred', p_scale, 1001),   \
                Sniper_Score            = rescale(score_somaticsniper, 'phred', p_scale, 1001),   \
                VarDict_Score           = rescale(score_vardict,       'phred', p_scale, 1001),   \
                if_dbsnp                = if_dbsnp,                                               \
                COMMON                  = if_common,                                              \
                if_COSMIC               = if_cosmic,                                              \
                COSMIC_CNT              = num_cases,                                              \
                N_DP                    = N_dp,                                                   \
                nBAM_REF_MQ             = '%g' % n_ref_mq,                                        \
                nBAM_ALT_MQ             = '%g' % n_alt_mq,                                        \
                nBAM_Z_Ranksums_MQ      = '%g' % n_z_ranksums_mq,                                 \
                nBAM_REF_BQ             = '%g' % n_ref_bq,                                        \
                nBAM_ALT_BQ             = '%g' % n_alt_bq,                                        \
                nBAM_Z_Ranksums_BQ      = '%g' % n_z_ranksums_bq,                                 \
                nBAM_REF_NM             = '%g' % n_ref_NM,                                        \
                nBAM_ALT_NM             = '%g' % n_alt_NM,                                        \
                nBAM_NM_Diff            = '%g' % n_NM_Diff,                                       \
                nBAM_REF_Concordant     = n_ref_concordant_reads,                                 \
                nBAM_REF_Discordant     = n_ref_discordant_reads,                                 \
                nBAM_ALT_Concordant     = n_alt_concordant_reads,                                 \
                nBAM_ALT_Discordant     = n_alt_discordant_reads,                                 \
                nBAM_Concordance_FET    = rescale(n_concordance_fet, 'fraction', p_scale, 1001),  \
                N_REF_FOR               = n_ref_for,                                              \
                N_REF_REV               = n_ref_rev,                                              \
                N_ALT_FOR               = n_alt_for,                                              \
                N_ALT_REV               = n_alt_rev,                                              \
                nBAM_StrandBias_FET     = rescale(n_strandbias_fet, 'fraction', p_scale, 1001),   \
                nBAM_Z_Ranksums_EndPos  = '%g' % n_z_ranksums_endpos,                             \
                nBAM_REF_Clipped_Reads  = n_ref_SC_reads,                                         \
                nBAM_ALT_Clipped_Reads  = n_alt_SC_reads,                                         \
                nBAM_Clipping_FET       = rescale(n_clipping_fet, 'fraction', p_scale, 1001),     \
                nBAM_MQ0                = n_MQ0,                                                  \
                nBAM_Other_Reads        = n_noise_read_count,                                     \
                nBAM_Poor_Reads         = n_poor_read_count,                                      \
                nBAM_REF_InDel_3bp      = n_ref_indel_3bp,                                        \
                nBAM_REF_InDel_2bp      = n_ref_indel_2bp,                                        \
                nBAM_REF_InDel_1bp      = n_ref_indel_1bp,                                        \
                nBAM_ALT_InDel_3bp      = n_alt_indel_3bp,                                        \
                nBAM_ALT_InDel_2bp      = n_alt_indel_2bp,                                        \
                nBAM_ALT_InDel_1bp      = n_alt_indel_1bp,                                        \
                M2_RPA                  = rpa,                                                    \
                M2_NLOD                 = nlod,                                                   \
                M2_TLOD                 = tlod,                                                   \
                M2_STR                  = tandem,                                                 \
                M2_ECNT                 = ecnt,                                                   \
                M2_HCNT                 = hcnt,                                                   \
                M2_MAXED                = maxED,                                                  \
                M2_MINED                = minED,                                                  \
                SOR                     = sor,                                                    \
                MSI                     = msi,                                                    \
                MSILEN                  = msilen,                                                 \
                SHIFT3                  = shift3,                                                 \
                MaxHomopolymer_Length   = homopolymer_length,                                     \
                SiteHomopolymer_Length  = site_homopolymer_length,                                \
                T_DP                    = T_dp,                                                   \
                tBAM_REF_MQ             = '%g' % t_ref_mq,                                        \
                tBAM_ALT_MQ             = '%g' % t_alt_mq,                                        \
                tBAM_Z_Ranksums_MQ      = '%g' % t_z_ranksums_mq,                                 \
                tBAM_REF_BQ             = '%g' % t_ref_bq,                                        \
                tBAM_ALT_BQ             = '%g' % t_alt_bq,                                        \
                tBAM_Z_Ranksums_BQ      = '%g' % t_z_ranksums_bq,                                 \
                tBAM_REF_NM             = '%g' % t_ref_NM,                                        \
                tBAM_ALT_NM             = '%g' % t_alt_NM,                                        \
                tBAM_NM_Diff            = '%g' % t_NM_Diff,                                       \
                tBAM_REF_Concordant     = t_ref_concordant_reads,                                 \
                tBAM_REF_Discordant     = t_ref_discordant_reads,                                 \
                tBAM_ALT_Concordant     = t_alt_concordant_reads,                                 \
                tBAM_ALT_Discordant     = t_alt_discordant_reads,                                 \
                tBAM_Concordance_FET    = rescale(t_concordance_fet, 'fraction', p_scale, 1001),  \
                T_REF_FOR               = t_ref_for,                                              \
                T_REF_REV               = t_ref_rev,                                              \
                T_ALT_FOR               = t_alt_for,                                              \
                T_ALT_REV               = t_alt_rev,                                              \
                tBAM_StrandBias_FET     = rescale(t_strandbias_fet, 'fraction', p_scale, 1001),   \
                tBAM_Z_Ranksums_EndPos  = '%g' % t_z_ranksums_endpos,                             \
                tBAM_REF_Clipped_Reads  = t_ref_SC_reads,                                         \
                tBAM_ALT_Clipped_Reads  = t_alt_SC_reads,                                         \
                tBAM_Clipping_FET       = rescale(t_clipping_fet, 'fraction', p_scale, 1001),     \
                tBAM_MQ0                = t_MQ0,                                                  \
                tBAM_Other_Reads        = t_noise_read_count,                                     \
                tBAM_Poor_Reads         = t_poor_read_count,                                      \
                tBAM_REF_InDel_3bp      = t_ref_indel_3bp,                                        \
                tBAM_REF_InDel_2bp      = t_ref_indel_2bp,                                        \
                tBAM_REF_InDel_1bp      = t_ref_indel_1bp,                                        \
                tBAM_ALT_InDel_3bp      = t_alt_indel_3bp,                                        \
                tBAM_ALT_InDel_2bp      = t_alt_indel_2bp,                                        \
                tBAM_ALT_InDel_1bp      = t_alt_indel_1bp,                                        \
                InDel_Length            = indel_length,                                           \
                TrueVariant_or_False    = judgement )
                
                # Print it out to stdout:
                outhandle.write(out_line + '\n')
        
        # Read into the next line:
        if not is_vcf:
            my_line = my_sites.readline().rstrip()
        
    ##########  Close all open files if they were opened  ##########
    opened_files = (ref_fa, nbam, tbam, truth, cosmic, dbsnp)
    [opened_file.close() for opened_file in opened_files if opened_file]

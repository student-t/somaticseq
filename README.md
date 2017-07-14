<b>SomaticSeq: An ensemble approach to accurately detect somatic mutations</b>

* See http://bioinform.github.io/somaticseq/ for help and downloads.
* Detailed documentation is included in the package. It's located in [docs/Manual.pdf](docs/Manual.pdf "Documentation").
* Feel free to report issues and/or ask questions at the [Issues](../../issues "Issues") page.
* For a quick description of SomaticSeq, you may watch this 8-minute video:
  [![SomaticSeq Video](SomaticSeqYoutube.png)](https://www.youtube.com/watch?v=MnJdTQWWN6w "SomaticSeq Video")

* This is the specialized branch designed to consolidate multiple sequencing replicates of the same samples 

   SSeq_tsv2vcf.py is modified to move almost all information to the sample columns. Thus, when they are merged using GATK CombineVariants, those information will be retained in the sample columns of the VCF files.

   SSeq_vcf2tsv_multibam.py handles many BAM files. It will also extract sample-specific information from the VCF files produced above. 


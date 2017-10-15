# README #

Immune evasion is a hallmark of cancer. Losing the ability to present productive tumor neoantigens could facilitate evasion from immune predation. 
An integral part of neoantigen presentation is the HLA class I molecule, which presents epitopes to T-cells on the cell surface. Thus, loss of an 
HLA allele, resulting in HLA homozygosity, may be a mechanism of immune escape. However, the polymorphic nature of the HLA locus precludes accurate
copy number calling using conventional copy number tools.  

Here, we present **LOHHLA**, **L**oss **O**f **H**eterozygosity in **H**uman **L**eukocyte **A**ntigen, a computational tool to evaluate HLA loss 
using next-generation sequencing data. 


### What do I need to install? ###

Please ensure a number of dependencies are first installed. These include:

* BEDTools (http://bedtools.readthedocs.io/en/latest/)
* SAMtools (http://samtools.sourceforge.net/)
* R (https://www.r-project.org/about.html)
* Novalign (http://www.novocraft.com/products/novoalign/)
* Picard (http://broadinstitute.github.io/picard/)


### How do I run LOHHLA? ###

LOHHLA is coded in R, and should be executed from the command line (Terminal, in Linux/UNIX/OSX, or Command Prompt in MS Windows). 
Running LOHHLA with no arguments prints the usage information. 

USAGE: Rscript /location/of/LOHHLA/script  [OPTIONS]

OPTIONS:
*   -id / --patientId            patient ID                         type="character"              default=NULL
*   -o  / --outputDir            location of output directory       type="character"              default=NULL
*
              help="location of output directory", metavar="character"),
  make_option(c("-nBAM", "--normalBAMfile"), type="character", default=NULL, 
              help="normal BAM file\n\t\tcan be FALSE to run without normal sample", metavar="character"),
  make_option(c("-BAM", "--BAMDir"), type="character", default=NULL, 
              help="location of all BAMs to test", metavar="character"),
  make_option(c("-hla", "--hlaPath"), type="character", default=NULL, 
              help="location to patient HLA calls", metavar="character"),
  make_option(c("-hlaLoc", "--HLAfastaLoc"), type="character", default="/farm/home/lr-tct-lif/wilson52/installs/polysolver/data/abc_complete.fasta", 
              help="location of HLA FASTA [default= %default]", metavar="character"),
  make_option(c("-cn", "--CopyNumLoc"), type="character", default="FALSE", 
              help="location to patient purity and ploidy output\n\t\tcan be FALSE to only estimate allelic imbalance", metavar="character"),
  make_option(c("-ov", "--overrideDir"), type="character", default='FALSE', 
              help="location of flagstat information if already run [default= %default]", metavar="character"),
  make_option(c("-mc", "--minCoverageFilter"), type="numeric", default=30, 
              help="minimum coverage at mismatch site [default= %default]", metavar="character"),
  make_option(c("-mm", "--numMisMatch"), type="numeric", default=1, 
              help="number of mismatches allowed in read to map to HLA allele [default= %default]", metavar="character"),
  make_option(c("-m", "--mappingStep"), type="logical", default=TRUE, 
              help="does mapping to HLA alleles need to be done [default= %default]", metavar="character"),
  make_option(c("-cu", "--cleanUp"), type="logical", default=TRUE, 
              help="remove temporary files [default= %default]", metavar="character"),
  make_option(c("-no", "--novoDir"), type="character", default='', 
              help="path to novoalign executable [default= %default]", metavar="character"),
  make_option(c("-ga", "--gatkDIR"), type="character", default='', 
              help="path to GATK executable [default= %default]", metavar="character"),
  make_option(c("-ex", "--HLAexonLoc"), type="character", default='', 
              help="HLA exon boundaries for plotting [default= %default]", metavar="character")


### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
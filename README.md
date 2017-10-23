# README #

Immune evasion is a hallmark of cancer. Losing the ability to present productive tumor neoantigens could facilitate evasion from immune predation. 
An integral part of neoantigen presentation is the HLA class I molecule, which presents epitopes to T-cells on the cell surface. Thus, loss of an 
HLA allele, resulting in HLA homozygosity, may be a mechanism of immune escape. However, the polymorphic nature of the HLA locus precludes accurate
copy number calling using conventional copy number tools.  

Here, we present **LOHHLA**, **L**oss **O**f **H**eterozygosity in **H**uman **L**eukocyte **A**ntigen, a computational tool to evaluate HLA loss 
using next-generation sequencing data. 


### What do I need to install to run LOHHLA? ###

Please ensure a number of dependencies are first installed. These include:

* BEDTools (http://bedtools.readthedocs.io/en/latest/)
* SAMtools (http://samtools.sourceforge.net/)
* Novalign (http://www.novocraft.com/products/novoalign/)
* Picard (http://broadinstitute.github.io/picard/)
* R (https://www.r-project.org/about.html)

Within R, the following packages are required:

* seqinr (https://CRAN.R-project.org/package=seqinr)
* Biostrings (http://bioconductor.org/packages/release/bioc/html/Biostrings.html)
* beeswarm (https://CRAN.R-project.org/package=beeswarm)
* zoo (https://cran.r-project.org/package=zoo)
* Rsamtools (http://bioconductor.org/packages/release/bioc/html/Rsamtools.html)

### How do I install LOHHLA? ###

To install LOHHLA, simple clone the repository:

git clone https://nmcgranahan@bitbucket.org/mcgranahanlab/lohhla.git

### How do I run LOHHLA? ###

LOHHLA is coded in R, and can be executed from the command line (Terminal, in Linux/UNIX/OSX, or Command Prompt in MS Windows) directly, 
or using a shell script (see example below).

Running LOHHLA with no arguments prints the usage information. 

USAGE: Rscript /location/of/LOHHLA/script  [OPTIONS]

OPTIONS:

	-id CHARACTER, --patientId=CHARACTER
		patient ID

	-o CHARACTER, --outputDir=CHARACTER
		location of output directory

	-nBAM CHARACTER, --normalBAMfile=CHARACTER
		normal BAM file
		can be FALSE to run without normal sample

	-BAM CHARACTER, --BAMDir=CHARACTER
		location of all BAMs to test

	-hla CHARACTER, --hlaPath=CHARACTER
		location to patient HLA calls

	-hlaLoc CHARACTER, --HLAfastaLoc=CHARACTER
		location of HLA FASTA [default= /farm/home/lr-tct-lif/wilson52/installs/polysolver/data/abc_complete.fasta]

	-cn CHARACTER, --CopyNumLoc=CHARACTER
		location to patient purity and ploidy output
		can be FALSE to only estimate allelic imbalance

	-ov CHARACTER, --overrideDir=CHARACTER
		location of flagstat information if already run [default= FALSE]

	-mc CHARACTER, --minCoverageFilter=CHARACTER
		minimum coverage at mismatch site [default= 30]

	-kmer CHARACTER, --kmerSize=CHARACTER
		size of kmers to fish with [default= 50]

	-mm CHARACTER, --numMisMatch=CHARACTER
		number of mismatches allowed in read to map to HLA allele [default= 1]

	-ms CHARACTER, --mappingStep=CHARACTER
		does mapping to HLA alleles need to be done [default= TRUE]

	-fs CHARACTER, --fishingStep=CHARACTER
		if mapping is performed, also look for fished reads matching kmers of size kmerSize [default= TRUE]

	-ps CHARACTER, --plottingStep=CHARACTER
		are plots made [default= TRUE]

	-cs CHARACTER, --coverageStep=CHARACTER
		are coverage differences analyzed [default= TRUE]

	-cu CHARACTER, --cleanUp=CHARACTER
		remove temporary files [default= TRUE]

	-no CHARACTER, --novoDir=CHARACTER
		path to novoalign executable [default= ]

	-ga CHARACTER, --gatkDir=CHARACTER
		path to GATK executable [default= ]

	-ex CHARACTER, --HLAexonLoc=CHARACTER
		HLA exon boundaries for plotting [default= /camp/lab/swantonc/working/rosentr/data/IMGT/hla.dat]

	-w CHARACTER, --ignoreWarnings=CHARACTER
		continue running with warnings [default= TRUE]

	-h, --help
		Show this help message and exit            
 

### How can I test if LOHHLA is working? ###

Example data is included in the LOHHLA repository. Alter the "example.sh" script to match your local file structure and ensure the requisite dependencies are available / loaded.
The --HLAfastaLoc, --gatkDir, and --novoDir file paths should also be updated to the corresponding locations.
File paths must be full paths. Run "example.sh" and the output should match that found in the "correct-example-out" directory provided.


### Who do I talk to? ###

If you have any issues with lohhla, please send an email to lohhla@gmail.com

### How do I cite LOHHLA ? ###

If you use LOHHLA in your research, please cite the following paper:


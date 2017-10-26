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

LOHHLA also requires an HLA fasta file. This can be obtained from Polysolver (http://archive.broadinstitute.org/cancer/cga/polysolver)

### How do I install LOHHLA? ###

To install LOHHLA, simply clone the repository:

git clone https://bitbucket.org/mcgranahanlab/lohhla.git

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
		location of HLA FASTA [default=~/lohhla/data/hla_all.fasta]

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
		HLA exon boundaries for plotting [default=~/lohhla/data/hla.dat]

	-w CHARACTER, --ignoreWarnings=CHARACTER
		continue running with warnings [default= TRUE]

	-h, --help
		Show this help message and exit            
 

### What is the output of LOHHLA? ###

LOHHLA produces multiple different files (see correct-example-out for an example). To determine HLA LOH in a given sample, the most relevant output is the file which ends '.HLAlossPrediction CI.xls'. 
The most relavant columns are:

	HLA_A_type1  						 - the identity of allele 1
	HLA_A_type2  						 - the identity of allele 2
	Pval_unique  					     - this is a p-value relating to allelic imbalance 
	LossAllele      					 - this corresponds to the HLA allele that is subject to loss
	KeptAllele      					 - this corresponds to the HLA allele that is not subject to loss
	HLA_type1copyNum_withBAFBin          - the estimated raw copy number of HLA (allele 1)
	HLA_type2copyNum_withBAFBin          - the estimated raw copy number of HLA (allele 2)




### How can I test if LOHHLA is working? ###

Example data is included in the LOHHLA repository. To run LOHHLA on the example dataset, alter the "example.sh" script to match your local file structure and ensure the requisite dependencies are available / loaded.
The --HLAfastaLoc, --gatkDir, and --novoDir file paths should also be updated to the corresponding locations.
File paths must be full paths. Run "example.sh" and the output should match that found in the "correct-example-out" directory provided.
All BAM files (normal and tumour) should be found in or linked to the same directory.

### Who do I talk to? ###

If you have any issues with lohhla, please send an email to lohhla@gmail.com

### How do I cite LOHHLA ? ###

If you use LOHHLA in your research, please cite the following paper:

McGranahan et al., Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution, Cell (2017), https://doi.org/10.1016/j.cell.2017.10.001


# README #

Immune evasion is a hallmark of cancer. Losing the ability to present productive tumor neoantigens could facilitate evasion from immune predation. 
An integral part of neoantigen presentation is the HLA class I molecule, which presents epitopes to T-cells on the cell surface. Thus, loss of an 
HLA allele, resulting in HLA homozygosity, may be a mechanism of immune escape. However, the polymorphic nature of the HLA locus precludes accurate
copy number calling using conventional copy number tools.  

Here, we present **LOHHLA**, **L**oss **O**f **H**eterozygosity in **H**uman **L**eukocyte **A**ntigen, a computational tool to evaluate HLA loss 
using next-generation sequencing data. 

### LICENCE ###

THIS LOHHLA TOOL IS PROTECTED BY COPYRIGHT AND IS SUBJECT TO PATENT APPLICATION PCT/GB2018/052004.  THE TOOL IS PROVIDED “AS IS” FOR INTERNAL 
NON-COMMERCIAL ACADEMIC RESEARCH PURPOSES ONLY.  THIS LOHHLA TOOL IS PROVIDED TO YOU AS A PERSONAL COPY SUBJECT TO THESE TERMS AND YOU SHALL NOT 
FORWARD OR PROVIDE A COPY OF THE LOHHLA TOOL TO ANY OTHER PERSON OR PARTY.   NO RESPONSIBILITY IS ACCEPTED FOR ANY LIABILITY ARISING FROM USE OF 
THE LOHHLA TOOL OR ANY RESULTS ARISING THEREON BY YOU OR ANY THIRD PARTY.  COMMERCIAL USE OF THIS TOOL FOR ANY PURPOSE IS NOT PERMITTED.  
ALL COMMERCIAL USE OF THE TOOL OR ANY MODIFICATION OR DERIVATIVE OF THE TOOL, INCLUDING TRANSFER, SALE OR LICENCE TO A COMMERCIAL THIRD PARTY OR 
USE ON BEHALF OF A COMMERCIAL THIRD PARTY (INCLUDING BUT NOT LIMITED TO USE AS PART OF A SERVICE SUPPLIED TO ANY THIRD PARTY FOR FINANCIAL REWARD) 
REQUIRES A LICENCE.  FOR FURTHER INFORMATION PLEASE EMAIL Eileen Clark eileen.clark@crick.ac.uk. 
 
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


For a full definition of the columns, see below, in each case whether the column should be used [use], or can be ignored [legacy]is indicated:

	region								 - the region or tumor sample [use]
	HLA_A_type1							 - the identity of allele 1 [use]
	HLA_A_type2							 - the identity of allele 2 [use]
	HLAtype1Log2MedianCoverage	         - the median LogR coverage across allele 1 [use] 
	HLAtype2Log2MedianCoverage	         - the median LogR coverage across allele 2 [use]
	HLAtype1Log2MedianCoverageAtSites	 - the median LogR coverage across allele 1, restricted to mismatch sites [use]
	HLAtype2Log2MedianCoverageAtSites	 - the median LogR coverage across allele 2, restricted to mismatch sites [use]
	HLA_type1copyNum_withoutBAF	         - estimated copy number of allele 1, without using BAF [legacy] 
	HLA_type1copyNum_withoutBAF_lower	 - lower 95% confidence interval of estimated copy number of allele 1, without using BAF [legacy] 
	HLA_type1copyNum_withoutBAF_upper	 - upper 95% confidence interval of estimated copy number of allele 1, without using BAF [legacy] 
	HLA_type1copyNum_withBAF	         - estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type1copyNum_withBAF_lower	     - lower 95% confidence interval of estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type1copyNum_withBAF_upper	     - upper 95% confidence interval of estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type2copyNum_withoutBAF	         - estimated copy number of allele 2 without using BAF  [legacy] 
	HLA_type2copyNum_withoutBAF_lower	 - lower 95% confidence interval of estimated copy number of allele 2, without using BAF [legacy] 
	HLA_type2copyNum_withoutBAF_upper	 - upper 95% confidence interval of estimated copy number of allele 2, without using BAF [legacy] 
	HLA_type2copyNum_withBAF	         - estimated copy number of allele 2 using BAF, without binning sites [legacy] 
	HLA_type2copyNum_withBAF_lower	     - lower 95% confidence interval of estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type2copyNum_withBAF_upper	     - upper 95% confidence interval of estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type1copyNum_withoutBAFBin	     - estimated copy number of allele 1 using binning, but without BAF [legacy]  
	HLA_type1copyNum_withoutBAFBin_lower - lower 95% confidence interval of estimated copy number of allele 1 using binning, but without BAF [legacy]  	
	HLA_type1copyNum_withoutBAFBin_upper - upper 95% confidence interval of estimated copy number of allele 1 using binning, but without BAF [legacy] 	
	HLA_type1copyNum_withBAFBin	         - estimated copy number of allele 1 using binning and BAF [use] 
	HLA_type1copyNum_withBAFBin_lower	 - lower 95% confidence interval of estimated copy number of allele 1 using binning and BAF [use] 
	HLA_type1copyNum_withBAFBin_upper	 - upper 95% confidence interval of estimated copy number of allele 1 using binning and BAF [use]  
	HLA_type2copyNum_withoutBAFBin	     - estimated copy number of allele 2 using binning, but without BAF [legacy]  
	HLA_type2copyNum_withoutBAFBin_lower - lower 95% confidence interval of estimated copy number of allele 2 using binning, but without BAF [legacy]	
	HLA_type2copyNum_withoutBAFBin_upper - upper 95% confidence interval of estimated copy number of allele 2 using BAF, without binning sites [legacy] 	
	HLA_type2copyNum_withBAFBin	         - estimated copy number of allele 2 using binning and BAF [use] 
	HLA_type2copyNum_withBAFBin_lower	 - lower 95% confidence interval of estimated copy number of allele 2 using binning and BAF [use] 
	HLA_type2copyNum_withBAFBin_upper	 - upper 95% confidence interval of estimated copy number of allele 2 using binning and BAF [use
	PVal                                 - p-value relating to difference in logR between allele 1 and allele 2 (paired t-test)[legacy]
	UnPairedPval	                     - p-value relating to difference in logR between allele 1 and allele 2 (unpaired t-test)[legacy]
	PVal_unique	                         - p-value relating to difference in logR between allele 1 and allele 2, ensuring each read only contributes once (paired t-test) [use]
	UnPairedPval_unique                  - p-value relating to difference in logR between allele 1 and allele 2, ensuring each read only contributes once (unpaired t-test) [use]
	LossAllele	                         - HLA allele that is present at lower frequency (potentially subject to loss) [use]
	KeptAllele                           - HLA allele that is present at higher frequency (potentially not subject to loss) [use]
	numMisMatchSitesCov                  - number of mismatch sites with sufficient coverage [use]
	propSupportiveSites                  - proportion of missmatch sites that are consistent with loss or allelic imbalance [use]


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


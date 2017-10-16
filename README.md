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
* R (https://www.r-project.org/about.html)
* Novalign (http://www.novocraft.com/products/novoalign/)
* Picard (http://broadinstitute.github.io/picard/)


### How do I run LOHHLA? ###

LOHHLA is coded in R, and should be executed from the command line (Terminal, in Linux/UNIX/OSX, or Command Prompt in MS Windows). 
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

        -mm CHARACTER, --numMisMatch=CHARACTER
                number of mismatches allowed in read to map to HLA allele [default= 1]

        -m CHARACTER, --mappingStep=CHARACTER
                does mapping to HLA alleles need to be done [default= TRUE]

        -cu CHARACTER, --cleanUp=CHARACTER
                remove temporary files [default= TRUE]

        -no CHARACTER, --novoDir=CHARACTER
                path to novoalign executable [default= ]

        -ga CHARACTER, --gatkDIR=CHARACTER
                path to GATK executable [default= ]

        -ex CHARACTER, --HLAexonLoc=CHARACTER
                HLA exon boundaries for plotting [default= ]

        -h, --help
                Show this help message and exit


### How can I test if LOHHLA is working? ###


### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
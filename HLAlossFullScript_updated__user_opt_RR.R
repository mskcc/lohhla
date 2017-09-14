# before running
# ml BEDTools/2.26.0-foss-2016b
# ml SAMtools/1.3.1-foss-2016b
# ml R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3
# ml novoalign/3.07.00
# ml TracerX-Picard-GATK/0.1-Java-1.7.0_80

library(optparse)
option_list = list(
  make_option(c("-id", "--patientId"), type="character", default=NULL, 
              help="patient ID", metavar="character"),
  make_option(c("-o", "--outputDir"), type="character", default=NULL, 
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
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

print(opt)

##########################
# Command line arguments #
##########################

full.patient      <- opt$patientId
workDir           <- opt$outputDir
normalBAMfile     <- opt$normalBAMfile
BAMDir            <- opt$BAMDir
hlaPath           <- opt$hlaPath
HLAfastaLoc       <- opt$HLAfastaLoc
CopyNumLoc        <- opt$CopyNumLoc
overrideDir       <- opt$overrideDir
minCoverageFilter <- opt$minCoverageFilter
numMisMatch       <- opt$numMisMatch
mapping.step      <- opt$mappingStep
cleanUp           <- opt$cleanUp
NOVODir           <- opt$novoDir
GATKDir           <- opt$gatkDIR
HLAexonLoc        <- opt$HLAexonLoc

if (is.null(opt$BAMDir) | is.null(opt$outputDir) | is.null(opt$hlaPath) | is.null(opt$HLAfastaLoc)){
  print_help(opt_parser)
  stop("Missing arguments.\n", call.=FALSE)  
}




#############
# libraries #
#############

require(seqinr, quietly = TRUE)
require(Biostrings, quietly = TRUE)
require(beeswarm, quietly = TRUE)
require(zoo, quietly = TRUE)
require(Rsamtools, quietly = TRUE)



###########
# inputs #
###########

interactive           <- FALSE
if(interactive)
{
  full.patient          <- "L_LTX050"
  workDir               <- "/camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/tx100-noPoly/L_LTX050/exome/NeoAntigen/LOH/"
  hlaPath               <- '/camp/lab/swantonc/working/rosentr/projects/neoantigen/tx100/samples-20160818/LTX050/LTX050.polysolver/winners.hla.txt'
  normalBAMfile         <- '/farm/tracerx/lung/release_002.2/L_LTX050/exome/BAM/processed/L_LTX050_BS_GL.bam'
  BAMDir                <- '/farm/tracerx/lung/release_002.2/L_LTX050/exome/BAM/processed/'
  HLAfastaLoc           <- "/farm/home/lr-tct-lif/wilson52/installs/polysolver/data/abc_complete.fasta"
  mapping.step          <- TRUE
  cleanUp               <- FALSE
  CopyNumLoc            <- '/farm/tracerx/lung/release_002.2/L_LTX050/exome/ASCAT/solutions.txt'
  overrideDir           <- '/farm/tracerx/lung/release/L_LTX050/exome/QC/flagstat/'

  full.patient          <- "U_LTX156"
  workDir               <- "/camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/tx100-noPoly/U_LTX156/exome/NeoAntigen/LOH/"
  hlaPath               <- '/camp/lab/swantonc/working/rosentr/projects/neoantigen/tx100/samples-20160818/LTX156/LTX156.polysolver/winners.hla.txt'
  normalBAMfile         <- '/farm/tracerx/lung/release_002.2/U_LTX156/exome/BAM/processed/U_LTX156_BS_GL.bam'
  BAMDir                <- '/farm/tracerx/lung/release_002.2/U_LTX156/exome/BAM/processed/'
  HLAfastaLoc           <- "/farm/home/lr-tct-lif/wilson52/installs/polysolver/data/abc_complete.fasta"
  mapping.step          <- TRUE
  cleanUp               <- FALSE
  CopyNumLoc            <- '/farm/tracerx/lung/release_002.2/U_LTX156/exome/ASCAT/solutions.txt'
  overrideDir           <- '/farm/tracerx/lung/release/U_LTX156/exome/QC/flagstat/'

  full.patient          <- "D_LMS021"
  workDir               <- "/camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/Rizvi-noPoly/D_LMS021//exome/NeoAntigen/LOH/"
  hlaPath               <- '/camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/Rizvi/LMS021/exome/NeoAntigen/Polysolver/winners.hla.txt'
  normalBAMfile         <- '/farm/tracerx/lung/LMS/lms_280715/D_LMS021/exome/BAM/processed/D_LMS021_GL.bam'
  BAMDir                <- '/farm/tracerx/lung/LMS/lms_280715/D_LMS021/exome/BAM/processed/'
  HLAfastaLoc           <- "/farm/home/lr-tct-lif/wilson52/installs/polysolver/data/abc_complete.fasta"
  mapping.step          <- TRUE
  cleanUp               <- FALSE
  CopyNumLoc            <- '/farm/tracerx/lung/LMS/lms_280715/D_LMS021/exome/ASCAT/solutions.txt'
  overrideDir           <- '/farm/tracerx/lung/LMS/lms_280715/D_LMS021/exome/QC/flagstat/'

  full.patient          <- "A_LTX049"
  workDir               <- "/camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/tx100-noPoly/A_LTX049/rna/NeoAntigen/LOH/"
  hlaPath               <- '/camp/lab/swantonc/working/rosentr/projects/neoantigen/tx100/samples-20160818/LTX049/LTX049.polysolver/winners.hla.txt'
  normalBAMfile         <- 'FALSE'
  BAMDir                <- '/camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/tx100-noPoly/A_LTX049/rna/BAM/'
  HLAfastaLoc           <- "/camp/lab/swantonc/working/rosentr/data/IMGT/hla_abc_complete.rna.fasta"
  mapping.step          <- TRUE
  cleanUp               <- FALSE
  CopyNumLoc            <- 'FALSE'
}

print(full.patient)
system('echo ${SLURM_JOBID}')

figureDir <- paste(workDir,"/Figures/",sep="")

log.name  <- paste(workDir, "running.hla.loh.exome@", gsub(":", "-", gsub(" +", "_", date())), "_log.txt", sep="")


runWithNormal           <- TRUE
if(normalBAMfile == 'FALSE'){
  runWithNormal         <- FALSE
}

extractNONmismatchReads <- TRUE
extractUniqueReads      <- TRUE

performIntegerCopyNum   <- TRUE
useLogRbin              <- TRUE
if(CopyNumLoc == 'FALSE'){
  performIntegerCopyNum <- FALSE
  useLogRbin            <- FALSE

}

override <-  ifelse(overrideDir == FALSE, yes = FALSE, no = TRUE)

gamma                   <- 1
binSize                 <- 150


#############
# functions #
#############

document.params <- function(params, log.name){
  
  msg <- "\n#######################\ \n####### Inputs ######## \n#######################\n"
  
  for(i in seq_along(params))
  {
    tmp <- c()
    for(j in 1:length(params[[i]])){
      tmp <- paste(tmp, params[[i]][j], "\t", sep = "")
    }
    msg <- paste(msg, names(params)[i], "\t", tmp, "\n", sep="")
  }
  
  write.table(msg, file=log.name, row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
  
}

PasteVector <- function(v,sep=""){
  
  vt <- v[1];
  if(length(v) > 1){
    for(g in 2:length(v)){
      vt <- paste(vt,v[g],sep=sep)
      
    }
  }
  vt <- paste(vt," EnD",sep="");
  out.v <- sub(" EnD","",vt);
  out.v <- sub("NA , ","",out.v);
  out.v <- sub(" , NA","",out.v);
  out.v <- sub(" , NA , "," , ",out.v);
  return(out.v);
  
}

count.events <- function(BAMfile, n){
  x              <- scanBam(BAMfile, index = BAMfile, param=ScanBamParam(what = scanBamWhat(), tag = 'NM'))
  readIDs        <- x[[1]][['qname']]
  cigar          <- x[[1]][['cigar']]
  editDistance   <- unlist(x[[1]][['tag']])
  insertionCount <- sapply(cigar, FUN = function(boop) {return(length(grep(pattern = 'I', x = unlist(strsplit(boop, split = '')))))} )
  deletionCount  <- sapply(cigar, FUN = function(boop) {return(length(grep(pattern = 'D', x = unlist(strsplit(boop, split = '')))))} )
  indelTotals    <- sapply(cigar, FUN = function(boop) {
    tmp <- unlist(strsplit( gsub("([0-9]+)","~\\1~",boop), "~" ))
    Is  <- grep(pattern = 'I', x = tmp)
    Ds  <- grep(pattern = 'D', x = tmp)
    total <- sum(as.numeric(tmp[(Is-1)])) + sum(as.numeric(tmp[Ds-1]))
    return(total)
  })
  misMatchCount <- editDistance - indelTotals
  eventCount <- misMatchCount + insertionCount + deletionCount
  names(eventCount) <- 1:length(eventCount)
  passed     <- eventCount[which(eventCount <= n)]
  y <- readIDs[as.numeric(names(passed))]
  y <- names(table(y)[which(table(y) == 2)])
  return(y)
}

dont.count.twice <- function(BAMfile1, BAMfile2, normalBAMfile1, normalBAMfile2){

  x              <- scanBam(BAMfile1)
  readIDs_x      <- x[[1]][['qname']]
  y              <- scanBam(BAMfile2)
  readIDs_y      <- y[[1]][['qname']]
  
  table(readIDs_x %in% readIDs_y)
  
  xn             <- scanBam(normalBAMfile1)
  readIDs_xn     <- xn[[1]][['qname']]
  yn             <- scanBam(normalBAMfile2)
  readIDs_yn     <- yn[[1]][['qname']]

  remove.from.tumor  <- readIDs_x[which(readIDs_x %in% readIDs_y)]
  remove.from.normal <- readIDs_xn[which(readIDs_xn %in% readIDs_yn)] 
  
  out <- list(remove.from.tumor, remove.from.normal)
  names(out) <- c('remove.from.tumor', 'remove.from.normal')
  return(out)
  
}


getMisMatchPositionsPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE){
  
  seq1aln <- pattern(alignment) # Get the alignment for the first sequence
  seq2aln <- subject(alignment) # Get the alignment for the second sequence
  
  #let's check differences across the whole thing
  gapsSeq1        <- countPattern("-",as.character(seq1aln))
  seq1alnresidues <- length(unlist(strsplit(as.character(seq1aln),split="")))-gapsSeq1
  
  k <- 1
  seq1Positions   <- c()
  for (char in unlist(strsplit(as.character(seq1aln),split="")))
  {
    if(char%in%c('C','G','A','T'))
    {
      seq1Positions <- c(seq1Positions,k)
      k <- k+1
      next;
    }
    if(char%in%c('-'))
    {
      seq1Positions <- c(seq1Positions,k)
      next;
      #
    }
  }
  
  
  k <- 1
  seq2Positions   <- c()
  for (char in unlist(strsplit(as.character(seq2aln),split="")))
  {
    if(char%in%c('C','G','A','T'))
    {
      seq2Positions <- c(seq2Positions,k)
      k <- k+1
      next;
    }
    if(char%in%c('-'))
    {
      seq2Positions <- c(seq2Positions,k)
      next;
      #
    }
  }
  
  diffSeq1 <-  seq1Positions[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]
  diffSeq2 <-  seq2Positions[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]
  
  diffType1 <- rep(1,length(diffSeq1))
  diffType1[which(unlist(strsplit(as.character(seq1aln),split=""))[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]%in%'-')] <- 2
  
  diffType2 <- rep(1,length(diffSeq2))
  diffType2[which(unlist(strsplit(as.character(seq2aln),split=""))[unlist(strsplit(as.character(seq2aln),split=""))!=unlist(strsplit(as.character(seq1aln),split=""))]%in%'-')] <- 2
  
  
  out <- list()
  out$diffSeq1 <- diffSeq1
  out$diffSeq2 <- diffSeq2
  out$diffType1 <- diffType1
  out$diffType2 <- diffType2
  return(out)
  
}


getUniqMapReads <- function(workDir
                            ,BAMDir
                            ,override=FALSE
                            ,overrideDir = NULL
)
{

  if(!override){
    outDir     <- paste(workDir, '/flagstat/', sep = '')
    if( !file.exists(outDir)){
      if( !dir.create(outDir, recursive = TRUE) ){
        stop(paste("Unable to create directory: ",outDir, "!\n", sep = ''))
      }
    }

    BAMs       <- list.files(BAMDir, pattern = 'bam$', full.names = TRUE)

    for(BAM in BAMs){
      region <- unlist(strsplit(BAM, split = '/'))[length(unlist(strsplit(BAM, split = '/')))]
      cmd    <- paste('samtools flagstat ', BAM, ' > ', outDir, region, '.proc.flagstat', sep = '')
      system(cmd)
    }  
  }    

  if(override){

    outDir <- overrideDir

  }
  
  flagStatRegions  <- list.files(outDir,pattern=".proc.flagstat$")
  if(length(flagStatRegions) == 0){
    stop('Either run flagstat or do not override.')
  }  
  
  UniqMapReads <- list()
  
  for (flagStatRegion in flagStatRegions)
  {
    
    UniqMapReads[[unlist(strsplit(flagStatRegion,split="\\."))[1]]] <-as.numeric(read.table(paste(outDir,flagStatRegion,sep=""),stringsAsFactors=FALSE,header=FALSE,nrows =1)[,1])    
    
  }
  
  return(UniqMapReads)
  
}

funCalcN_withBAF <- function(logRSites,bafSites,tumorPloidy,tumorPurity,gamma)
{
  
  nA <- (tumorPurity-1+bafSites*2^(logRSites/gamma)*((1-tumorPurity)*2+tumorPurity*tumorPloidy))/tumorPurity
  nB <- (tumorPurity-1-(bafSites-1)*2^(logRSites/gamma)*((1-tumorPurity)*2+tumorPurity*tumorPloidy))/tumorPurity
  return(cbind(nA,nB))
  
}

funCalcN_withoutBAF <- function(rSites,tumorPloidy,tumorPurity,gamma)
{
  return(((((1-tumorPurity)+tumorPurity*tumorPloidy/2))*2^(rSites/gamma) - (1-tumorPurity))/tumorPurity )
}

t.test.NA <- function(x){
  if(length(x[!is.na(x)]) <= 1){
    stat <- NA
    ci   <- c(NA, NA)
    out <- list(stat, ci)
    names(out) <- c('stat', 'conf.int')
    return(out)
  }
  else{
    return(t.test(x))
  }
}

print('got here')
print(workDir)
#############################
# create output directories # 
#############################

if(!dir.exists(workDir))
{
  dir.create(workDir,recursive=TRUE)
}

if(!dir.exists(figureDir))
{
  dir.create(figureDir,recursive=TRUE)
}

params <- list(full.patient, workDir, hlaPath, normalBAMfile, BAMDir, HLAfastaLoc, CopyNumLoc, GATKDir, NOVODir)
names(params) <- c('full.patient', 'workDir', 'hlaPath', 'normalBAMfile', 'BAMDir', 'HLAfastaLoc', 'CopyNumLoc', 'GATKDir', 'NOVODir')
document.params(params, log.name)



#############################
# get HLA-specific mappings # 
#############################

BAMfiles  <- list.files(BAMDir, pattern = '.bam$')
regions   <- sapply(BAMfiles, FUN =function(x) {return(unlist(strsplit(x, split = '.bam'))[1])})

hlaAlleles <- read.table(hlaPath, sep = '\t', header = FALSE, as.is = TRUE)
hlaAlleles <- unique(sort(c(hlaAlleles$V2, hlaAlleles$V3)))

hlaFasta   <- read.fasta(HLAfastaLoc)

if(!all(hlaAlleles %in% names(hlaFasta))) {
  warning(paste('Missing HLAs from FASTA: ', paste(hlaAlleles[-which(hlaAlleles %in% names(hlaFasta))], collapse = ', '), '\nTrying to find alternative.', sep = ''))
  missing <- hlaAlleles[-which(hlaAlleles %in% names(hlaFasta))]
  for(i in missing){
    alt <- grep(pattern = i, x = names(hlaFasta), value = TRUE)[1]
    if(!is.na(alt)){
      warning(paste('Replacing: ', i, ' with ',alt,  '.', sep = ''))
      hlaAlleles[which(hlaAlleles == i)] <- alt
    }
  }
  hlaAlleles <- hlaAlleles[which(hlaAlleles %in% names(hlaFasta))]
}

# check for homozygous alleles here to save time on mapping step.
# also figure out if hla names will be uniformly 'hla_x'
if(length(grep('hla_a', x = hlaAlleles))== 1){
  write.table(paste('\nHomozygous for HLA-A -- not going to see any LOH here.', '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  hlaAlleles <- hlaAlleles[-grep('hla_a', x = hlaAlleles)]
}

if(length(grep('hla_b', x = hlaAlleles))== 1){
  write.table(paste('\nHomozygous for HLA-B -- not going to see any LOH here.', '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  hlaAlleles <- hlaAlleles[-grep('hla_b', x = hlaAlleles)]
}

if(length(grep('hla_c', x = hlaAlleles))== 1){
  write.table(paste('\nHomozygous for HLA-C -- not going to see any LOH here.', '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  hlaAlleles <- hlaAlleles[-grep('hla_c', x = hlaAlleles)]
}

if(length(hlaAlleles) == 0){
  stop('No suitable HLA alleles!')
}

if(mapping.step){

  # generate patient reference fasta
  write.table(paste('\ngenerate patient reference fasta at ', date(), '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  patient.hlaFasta <- hlaFasta[hlaAlleles]
  write.fasta(patient.hlaFasta, file = paste(workDir, full.patient, '.patient.hlaFasta.fa', sep = ''), names = names(patient.hlaFasta))

  # nix file for patient reference fasta -- novoalign
  write.table(paste('\nnix file for patient reference fasta at ', date(), '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  novoindexCMD <- paste('novoindex ', workDir, full.patient, '.patient.hlaFasta.nix', ' ', workDir, full.patient, '.patient.hlaFasta.fa', sep = '')
  write.table(novoindexCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  system(novoindexCMD)

  for(BAMfile in BAMfiles){
    
    BAMid <- unlist(strsplit(BAMfile, split = '.bam'))[1]
    
    if(paste(BAMDir, BAMfile, sep = '') == normalBAMfile){
      normalName <- BAMid
    }

    print(BAMid)
    
    regionDir <- paste(workDir, '/', BAMid, sep = '')
    if(!dir.exists(regionDir))
    {
      dir.create(regionDir,recursive=TRUE)
    }
    
    #extract HLA possible reads from BAM file
    write.table(paste('\nextract HLA possible reads from BAM file at ', date(), '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    samtoolsCMD <- paste("samtools view -H ", BAMDir, '/', BAMfile, " > " , regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6:29909037-29913661 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6:31321649-31324964 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6:31236526-31239869 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)

    # samtoolsCMD <- paste("samtools view -F 4 ",BAMDir, '/', BAMfile, ' >> ', regionDir, '/', BAMid, ".chr6region.sam ", sep = "")
    # write.table(paste(samtoolsCMD, '\n', sep = ''), file = log.name, row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    # system(samtoolsCMD)

    samtoolsCMD <- paste("samtools view -f 4 ",BAMDir, '/', BAMfile, ' >> ', regionDir, '/', BAMid, ".chr6region.sam ", sep = "")
    write.table(paste(samtoolsCMD, '\n', sep = ''), file = log.name, row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6_cox_hap2 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6_dbb_hap3 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6_mann_hap4 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6_mcf_hap5 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6_qbl_hap6 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    samtoolsCMD <- paste("samtools view ", BAMDir, '/', BAMfile, " chr6_ssto_hap7 >> ",regionDir,"/",BAMid,".chr6region.sam",sep="")
    write.table(samtoolsCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samtoolsCMD)
    
    # turn into fastq -- this step has an error with unpaired mates, but seems to work ok just the same (VALIDATION_STRINGENCY=SILENT)
    write.table(paste('\nturn into fastq at ', date(), '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    samToFastQ <- paste("java -jar ",GATKDir,"/SamToFastq.jar ","I=",regionDir,"/",BAMid,".chr6region.sam"," F=",regionDir,"/",BAMid,".chr6region.1.fastq"," F2=",regionDir,"/",BAMid,".chr6region.2.fastq"," VALIDATION_STRINGENCY=SILENT",sep="")
    write.table(samToFastQ, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(samToFastQ)
    
    # align to all HLA alleles
    write.table(paste('\nalign to all HLA alleles at ', date(), '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)  
    
    alignCMD <- paste(NOVODir, '/novoalign -d ', workDir, full.patient, '.patient.hlaFasta.nix', ' -f ', regionDir,"/",BAMid,".chr6region.1.fastq", ' ', regionDir,"/",BAMid,".chr6region.2.fastq", ' -F STDFQ -R 0 -r All 9999 -o SAM -o FullNW 1> ', regionDir, '/', BAMid, '.chr6region.patient.reference.hlas.sam ', '2> ', regionDir, '/', BAMid, '_BS_GL.chr6region.patient.reference.hlas.metrics', sep = '')
    write.table(alignCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(alignCMD)

    convertToBam <- paste("samtools view -bS -o ",regionDir, '/', BAMid, '.chr6region.patient.reference.hlas.bam'," ",regionDir, '/', BAMid, '.chr6region.patient.reference.hlas.sam' , sep="")
    write.table(convertToBam, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(convertToBam)
    
    # sort
    sortBAM <- paste("java -jar ",GATKDir,"/SortSam.jar"," I=",regionDir, '/', BAMid, '.chr6region.patient.reference.hlas.bam'," ","O=",regionDir,"/",BAMid,".chr6region.patient.reference.hlas.csorted.bam", " SORT_ORDER=coordinate",sep="")
    write.table(sortBAM, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(sortBAM)
    
    # remove duplicates
    removeDup <- paste('samtools rmdup ', regionDir, '/', BAMid, '.chr6region.patient.reference.hlas.csorted.bam ', regionDir, '/', BAMid, '.chr6region.patient.reference.hlas.csorted.noduplicates.bam', sep = '')
    write.table(removeDup, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(removeDup)
    
    # only take reads that are in proper pair
    readPairs <- paste("samtools view -f 2 -b -o ",regionDir,"/",BAMid,".chr6region.patient.reference.hlas.csorted.noduplicates.filtered.bam"," ",regionDir,"/",BAMid,".chr6region.patient.reference.hlas.csorted.noduplicates.bam",sep="")
    write.table(readPairs, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(readPairs)
    
    # let's index the aligned bam
    indexBAM <- paste("samtools index ",regionDir,"/",BAMid,".chr6region.patient.reference.hlas.csorted.noduplicates.filtered.bam",sep="")
    write.table(indexBAM, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(indexBAM)
    
    hlaBAMfile <- paste(regionDir, '/', BAMid, '.chr6region.patient.reference.hlas.csorted.noduplicates.filtered.bam', sep = '')
    
    for (allele in hlaAlleles)
    {
      
      write.table(paste('\nget HLA specific SAM for allele: ', allele,  ' at ', date(), '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      
      getReads <- paste("samtools view -b -o ", regionDir,"/",BAMid,".temp.",allele, ".bam ", hlaBAMfile, " ", allele, sep = '')
      write.table(getReads, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      system(getReads)
      
      samtoolsSort <- paste("samtools sort ",regionDir,"/",BAMid,".temp.",allele, ".bam"," -o ",regionDir,"/",BAMid,".type.",allele, ".bam",sep="")
      write.table(samtoolsSort, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      system(samtoolsSort)  
      
      samtoolsIndex <- paste("samtools index ",regionDir,"/",BAMid,".type.",allele, ".bam",sep="")
      write.table(samtoolsIndex, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      system(samtoolsIndex)  
      
      # and filter out reads that have too many events -- has to be done here because some reads map to multiple alleles
      passed.reads <- count.events(paste(regionDir, '/', BAMid, '.type.', allele, '.bam', sep = ''), n = numMisMatch)
      write.table(passed.reads, file = paste(regionDir, '/', BAMid, '.', allele, '.passed.reads.txt', sep = ''), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
      extractCMD <- paste("java -jar ",GATKDir,"/FilterSamReads.jar ", "I=",  regionDir, '/', BAMid, '.type.', allele, '.bam' , " FILTER=includeReadList READ_LIST_FILE=", regionDir, "/", BAMid, '.', allele, ".passed.reads.txt", ' OUTPUT=', regionDir, '/', BAMid, '.type.', allele, '.filtered.bam', sep = '')
      write.table(extractCMD, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      system(extractCMD)

      samtoolsIndex <- paste("samtools index ",regionDir,"/",BAMid,".type.",allele, ".filtered.bam",sep="")
      write.table(samtoolsIndex, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      system(samtoolsIndex)
      
    }
    
  }

}



############################
# get coverage for regions # 
############################

for (region in regions)
{
  
  write.table(paste('\nget coverage of HLA alleles for region: ', region, ' at ', date(), '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
  regionDir <- paste(workDir, '/', region, sep = '')
  BAMfiles  <- grep('filtered.bam$', grep('type',list.files(regionDir),value=TRUE),value=TRUE)

  if(paste(BAMDir, region, '.bam', sep = '') == normalBAMfile){
    type <- 'normal'
  } else{
    type <- 'tumor'
  }
  
  #let's get pileup files for each bam
  for (BAMfile in c(BAMfiles))
  {      
    
    hlaAllele <- grep(pattern = 'hla', x = unlist(strsplit(BAMfile, split = '\\.')), value = TRUE)

    mpileupFile <- paste(workDir, "/",region,".",hlaAllele,".",type,".mpileup",sep="")
    
    #cmd         <- paste("samtools mpileup -q 20 -Q 20 ", regionDir, "/", BAMfile," -f ", HLAfastaLoc, " > ",mpileupFile,sep="")
    cmd         <- paste("samtools mpileup ", regionDir, "/", BAMfile," -f ", HLAfastaLoc, " > ",mpileupFile,sep="")

    write.table(cmd, file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    system(cmd)
    
  }

  
}

# also extract number of unique reads sequenced in tumor and normal
if(runWithNormal){
  
  if(!override){
    regionUniqMappedRegions <- getUniqMapReads(workDir = workDir, BAMDir = BAMDir, override = FALSE)
  }
  if(override){
   #regionUniqMappedRegions <- getUniqMapReads(workDir = workDir, BAMDir = BAMDir, override = TRUE, overrideDir = paste('/farm/tracerx/lung/release/', full.patient, '/exome/QC/flagstat/', sep = '')) 
    regionUniqMappedRegions <- getUniqMapReads(workDir = workDir, BAMDir = BAMDir, override = TRUE, overrideDir = overrideDir) 
  }

  GermLineUniqMappedReads <- regionUniqMappedRegions[[grep("GL",names(regionUniqMappedRegions),value=TRUE)]]

}


print(regionUniqMappedRegions)


####################################
# compare coverage between alleles # 
####################################

normalName <- regions[which(paste(BAMDir, regions, '.bam', sep = '') == normalBAMfile)]

# let's load the winners
# next, we can look at each mpileupFile, and assess whether we see differences in coverage between the two. 
# let's look at a region of interest. 
PatientOutPut <- c()
for (region in regions)
{
 
  if(paste(BAMDir, region, '.bam', sep = '') == normalBAMfile){
    next
  }
  
  if(runWithNormal){
    UniqMappedReads <- regionUniqMappedRegions[[region]]
    MultFactor      <- as.numeric(GermLineUniqMappedReads/UniqMappedReads)
  }

  if(!runWithNormal){
    MultFactor      <- 1
  }
  
  write.table(paste('\nanalyzing coverage differences in region: ', region, ' at ', date(), '\n', sep = ''), file = log.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
  print(region)
  HLAoutPut <- c()
  
  for (HLA_gene in c('hla_a','hla_b','hla_c'))
  {
    
    print(HLA_gene)
    
    HLA_As            <- grep(HLA_gene,hlaAlleles,value=TRUE)
    if(length(HLA_As)<=1)
    {
      next
    }
    HLA_A_type1 <- HLA_As[1]
    HLA_A_type2 <- HLA_As[2]
    
    # the reference sequence for the patient's two HLA alleles
    HLA_type1Fasta <- hlaFasta[[HLA_A_type1]]
    HLA_type2Fasta <- hlaFasta[[HLA_A_type2]]
    
    #perform local pairwise alignement#
    {

      seq1 <- PasteVector(toupper(HLA_type1Fasta),sep="")
      seq2 <- PasteVector(toupper(HLA_type2Fasta),sep="")
      sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
      
      
      tmp <- pairwiseAlignment(seq1, seq2, substitutionMatrix = sigma, gapOpening = -2,gapExtension = -4, scoreOnly = FALSE,type='local')
      
      missMatchPositions <- getMisMatchPositionsPairwiseAlignment(tmp,returnlist=TRUE)
      if(length(missMatchPositions$diffSeq1) == 0){
        next
      }
      
    }
    
    # getting exons for both alleles
    HLA_A_type1DatFormat <- unlist(strsplit(HLA_A_type1,split="_"))[2:length(unlist(strsplit(HLA_A_type1,split="_")))]
    HLA_A_type1Formatted <- toupper(paste('HLA-', HLA_A_type1DatFormat[1], '\\*', paste(HLA_A_type1DatFormat[2:length(HLA_A_type1DatFormat)], collapse = ':'), sep = ''))
    cmd        <- paste("awk \'/^DE   ", HLA_A_type1Formatted, ",/ {p=1}; p; /Sequence/ {p=0}\' ", HLAexonLoc, sep  = '')
    awk.result <- system(cmd, intern = TRUE)
    HLAtype1exons <- grep('^FT   exon',awk.result,value=TRUE)
    HLAtype1exons <- do.call(rbind,strsplit(HLAtype1exons,split=" "))[,ncol(do.call(rbind,strsplit(HLAtype1exons,split=" ")))]
    
    if(length(HLAtype1exons)!=0)
    {
      HLAtype1exonTable <- c()
      for (i in 1:length(HLAtype1exons))
      {
        HLAtype1exonTable <- rbind(HLAtype1exonTable,(unlist(strsplit(HLAtype1exons[i],split="\\.\\."))))
      }
    }
    
    HLA_A_type2DatFormat <- unlist(strsplit(HLA_A_type2,split="_"))[2:length(unlist(strsplit(HLA_A_type2,split="_")))]
    HLA_A_type2Formatted <- toupper(paste('HLA-', HLA_A_type2DatFormat[1], '\\*', paste(HLA_A_type2DatFormat[2:length(HLA_A_type2DatFormat)], collapse = ':'), sep = ''))
    cmd        <- paste("awk \'/^DE   ", HLA_A_type2Formatted, ",/ {p=1}; p; /Sequence/ {p=0}\' ", HLAexonLoc, sep  = '')
    awk.result <- system(cmd, intern = TRUE)
    HLAtype2exons        <- grep('^FT   exon',awk.result,value=TRUE)
    HLAtype2exons        <- do.call(rbind,strsplit(HLAtype2exons,split=" "))[,ncol(do.call(rbind,strsplit(HLAtype2exons,split=" ")))]
        
    if(length(HLAtype2exons)!=0)
    {
      HLAtype2exonTable <- c()
      for (i in 1:length(HLAtype2exons))
      {
        HLAtype2exonTable <- rbind(HLAtype2exonTable,(unlist(strsplit(HLAtype2exons[i],split="\\.\\."))))
      }
    }
    
        
    #load the normal and tumour for each type
    #HLAtype 1 coverage
    {
      
      HLA_A_type1normalLoc <- grep(pattern = HLA_A_type1, x = list.files(workDir, pattern = "normal.mpileup", full.names = TRUE), value = TRUE)
      if(runWithNormal){
        HLA_A_type1normal <- read.table(HLA_A_type1normalLoc ,sep="\t",stringsAsFactors=FALSE,quote="",fill=TRUE)
      }
      if(!runWithNormal){
        HLA_A_type1normal <- data.frame(cbind(HLA_A_type1, 1:length(HLA_type1Fasta), toupper(HLA_type1Fasta), minCoverageFilter+1), stringsAsFactors=FALSE)
        colnames(HLA_A_type1normal) <- paste('V',1:ncol(HLA_A_type1normal), sep = '')
        HLA_A_type1normal$V4 <- as.numeric(HLA_A_type1normal$V4)
      }
      rownames(HLA_A_type1normal) <- HLA_A_type1normal$V2
      HLA_A_type1tumor  <- read.table(paste(workDir, "/",region,".",HLA_A_type1,".","tumor.mpileup",sep=""),sep="\t",stringsAsFactors=FALSE,quote="",fill=TRUE)
      rownames(HLA_A_type1tumor) <- HLA_A_type1tumor$V2
      
      #apply minimum coverage thresholds (we only apply this to the normal for now)
      HLA_A_type1normal <- HLA_A_type1normal[HLA_A_type1normal$V4>minCoverageFilter,,drop=FALSE]
      if(nrow(HLA_A_type1normal) == 0){
        print('No position has greater than minimum coverage filter')
      }
      
      tmp <- intersect(rownames(HLA_A_type1tumor),rownames(HLA_A_type1normal))
      HLA_A_type1tumor  <- HLA_A_type1tumor[tmp,,drop=FALSE]
      # HLA_A_type1normal <- HLA_A_type1normal[tmp,,drop=FALSE]
      
      HLA_A_type1normalCov <- HLA_A_type1normal$V4#rep(0,max(HLA_A_type1normal$V2,HLA_A_type1tumor$V2))
      names(HLA_A_type1normalCov) <- HLA_A_type1normal$V2 #1:length(HLA_A_type1tumorCov)
      
      HLA_A_type1tumorCov <- rep(0,length(HLA_A_type1normalCov))      
      names(HLA_A_type1tumorCov) <- names(HLA_A_type1normalCov)
      HLA_A_type1tumorCov[rownames(HLA_A_type1tumor)] <- HLA_A_type1tumor$V4
      
    }
    
    #HLAtype2 coverage
    {
      
      HLA_A_type2normalLoc <- grep(pattern = HLA_A_type2, x = list.files(workDir, pattern = "normal.mpileup", full.names = TRUE), value = TRUE)
      if(runWithNormal){
        HLA_A_type2normal <- read.table(HLA_A_type2normalLoc ,sep="\t",stringsAsFactors=FALSE,quote="",fill=TRUE)
      }
      if(!runWithNormal){
        HLA_A_type2normal <- data.frame(cbind(HLA_A_type2, 1:length(HLA_type2Fasta), toupper(HLA_type2Fasta), minCoverageFilter+1), stringsAsFactors=FALSE)
        colnames(HLA_A_type2normal) <- paste('V',1:ncol(HLA_A_type2normal), sep = '')
        HLA_A_type2normal$V4 <- as.numeric(HLA_A_type2normal$V4)
      }
      rownames(HLA_A_type2normal) <- HLA_A_type2normal$V2
      HLA_A_type2tumor  <- read.table(paste(workDir, "/",region,".",HLA_A_type2,".","tumor.mpileup",sep=""),sep="\t",stringsAsFactors=FALSE,quote="",fill=TRUE)
      rownames(HLA_A_type2tumor) <- HLA_A_type2tumor$V2
      
      #apply minimum coverage thresholds (we only apply this to the normal for now)
      HLA_A_type2normal <- HLA_A_type2normal[HLA_A_type2normal$V4>minCoverageFilter,,drop=FALSE]
      
      tmp <- intersect(rownames(HLA_A_type2tumor),rownames(HLA_A_type2normal))
      HLA_A_type2tumor  <- HLA_A_type2tumor[tmp,,drop=FALSE]
      # HLA_A_type2normal <- HLA_A_type2normal[tmp,,drop=FALSE]
      
      HLA_A_type2normalCov <- HLA_A_type2normal$V4#rep(0,max(HLA_A_type2normal$V2,HLA_A_type2tumor$V2))
      names(HLA_A_type2normalCov) <- HLA_A_type2normal$V2 #1:length(HLA_A_type2tumorCov)
      
      HLA_A_type2tumorCov <- rep(0,length(HLA_A_type2normalCov))      
      names(HLA_A_type2tumorCov) <- names(HLA_A_type2normalCov)
      HLA_A_type2tumorCov[rownames(HLA_A_type2tumor)] <- HLA_A_type2tumor$V4
      
    }


    
    if(extractNONmismatchReads%in%TRUE)
    {
      #next, let's extract the reads that do not cover any mis-matches, and then mpileup them
      mismatchPosSeq1 <- HLA_A_type1tumor[HLA_A_type1tumor$V2%in%missMatchPositions$diffSeq1,,drop=FALSE]
      missMatchBed1    <- mismatchPosSeq1[,c(1,2,2)]
      missMatchBed1$V2 <- as.numeric(missMatchBed1$V2)-1
      missMatchBed1$V3 <- as.numeric(missMatchBed1$V2.1)+1
      
      mismatchPosSeq2 <- HLA_A_type2tumor[HLA_A_type2tumor$V2%in%missMatchPositions$diffSeq2,,drop=FALSE]
      missMatchBed2    <- mismatchPosSeq2[,c(1,2,2)]
      missMatchBed2$V2 <- as.numeric(missMatchBed2$V2)-1
      missMatchBed2$V3 <- as.numeric(missMatchBed2$V2.1)+1
      
      write.table(missMatchBed1,file=paste(workDir,region,".",HLA_A_type1,".bed",sep=""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
      write.table(missMatchBed2,file=paste(workDir,region,".",HLA_A_type2,".bed",sep=""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
      
      #create the command (assume bedtools has been loaded)
      Type1TumorCmd <- paste("bedtools intersect -v -a ",workDir,region,"/",region,".type.",HLA_A_type1,".filtered.bam"," -b ",workDir,region,".",HLA_A_type1,".bed"," > ",workDir,region,".",HLA_A_type1,".tumor.NoMissMatch.bam",sep="")
      system(Type1TumorCmd)
      if(runWithNormal){
        Type1NormalCmd <- paste("bedtools intersect -v -a ",workDir,normalName,"/",normalName,".type.",HLA_A_type1,".filtered.bam"," -b ",workDir,region,".",HLA_A_type1,".bed"," > ",workDir,region,".",HLA_A_type1,".normal.NoMissMatch.bam",sep="")
        system(Type1NormalCmd)
      }

      Type2TumorCmd <- paste("bedtools intersect -v -a ",workDir,region,"/",region, ".type.",HLA_A_type2,".filtered.bam"," -b ",workDir,region,".",HLA_A_type2,".bed"," > ",workDir,region,".",HLA_A_type2,".tumor.NoMissMatch.bam",sep="")
      system(Type2TumorCmd)
      if(runWithNormal){
        Type2NormalCmd <- paste("bedtools intersect -v -a ",workDir,normalName,"/",normalName,".type.",HLA_A_type2,".filtered.bam"," -b ",workDir,region,".",HLA_A_type2,".bed"," > ",workDir,region,".",HLA_A_type2,".normal.NoMissMatch.bam",sep="")
        system(Type2NormalCmd)
      }

      #next, let's do the mpileup step. 
      MpilupType1TumorCmd <- paste("samtools mpileup ",workDir,region,".",HLA_A_type1,".tumor.NoMissMatch.bam"," -f ", HLAfastaLoc, " > ",workDir,region,".",HLA_A_type1,".tumor.NoMissMatch.pileup",sep="")
      system(MpilupType1TumorCmd)
      if(runWithNormal){
        MpilupType1NormalCmd <- paste("samtools mpileup ",workDir,region,".",HLA_A_type1,".normal.NoMissMatch.bam"," -f ", HLAfastaLoc, " > ",workDir,region,".",HLA_A_type1,".normal.NoMissMatch.pileup",sep="")
        system(MpilupType1NormalCmd)
      }

      MpilupType2TumorCmd <- paste("samtools mpileup ",workDir,region,".",HLA_A_type2,".tumor.NoMissMatch.bam"," -f ", HLAfastaLoc, " > ",workDir,region,".",HLA_A_type2,".tumor.NoMissMatch.pileup",sep="")
      system(MpilupType2TumorCmd)
      if(runWithNormal){
        MpilupType2NormalCmd <- paste("samtools mpileup ",workDir,region,".",HLA_A_type2,".normal.NoMissMatch.bam"," -f ", HLAfastaLoc, " > ",workDir,region,".",HLA_A_type2,".normal.NoMissMatch.pileup",sep="")
        system(MpilupType2NormalCmd)
      }

      #get the coverage for the sites
      HLA_type1tumor_nomissmatch <- read.table(paste(workDir,region,".",HLA_A_type1,".tumor.NoMissMatch.pileup",sep=""),stringsAsFactors=FALSE,fill=TRUE,quote="", sep = '\t')
      HLA_type2tumor_nomissmatch <- read.table(paste(workDir,region,".",HLA_A_type2,".tumor.NoMissMatch.pileup",sep=""),stringsAsFactors=FALSE,fill=TRUE,quote="", sep = '\t')

      if(runWithNormal){
        HLA_type1normal_nomissmatch <- read.table(paste(workDir,region,".",HLA_A_type1,".normal.NoMissMatch.pileup",sep=""),stringsAsFactors=FALSE,fill=TRUE,quote="", sep = '\t')
        HLA_type2normal_nomissmatch <- read.table(paste(workDir,region,".",HLA_A_type2,".normal.NoMissMatch.pileup",sep=""),stringsAsFactors=FALSE,fill=TRUE,quote="", sep = '\t')
      }

      if(!runWithNormal){
        HLA_type1normal_nomissmatch <-  data.frame(cbind(HLA_A_type1, 1:length(HLA_type1Fasta), toupper(HLA_type1Fasta), minCoverageFilter+1), stringsAsFactors=FALSE)
        HLA_type1normal_nomissmatch$V4 <- as.numeric(HLA_type1normal_nomissmatch$V4)
        HLA_type2normal_nomissmatch <-  data.frame(cbind(HLA_A_type2, 1:length(HLA_type2Fasta), toupper(HLA_type2Fasta), minCoverageFilter+1), stringsAsFactors=FALSE)
        HLA_type2normal_nomissmatch$V4 <- as.numeric(HLA_type2normal_nomissmatch$V4)
      }

      rownames(HLA_type1tumor_nomissmatch) <- HLA_type1tumor_nomissmatch$V2
      rownames(HLA_type1normal_nomissmatch) <- HLA_type1normal_nomissmatch$V2
      rownames(HLA_type2tumor_nomissmatch) <- HLA_type2tumor_nomissmatch$V2
      rownames(HLA_type2normal_nomissmatch) <- HLA_type2normal_nomissmatch$V2
      
      #apply minimum coverage thresholds (we only apply this to the normal for now)
      HLA_type1normal_nomissmatch <- HLA_type1normal_nomissmatch[HLA_type1normal_nomissmatch$V4>minCoverageFilter,,drop=FALSE]
      HLA_type2normal_nomissmatch <- HLA_type2normal_nomissmatch[HLA_type2normal_nomissmatch$V4>minCoverageFilter,,drop=FALSE]
      
      tmp <- intersect(rownames(HLA_type1tumor_nomissmatch),rownames(HLA_type1normal_nomissmatch))
      HLA_type1tumor_nomissmatch  <- HLA_type1tumor_nomissmatch[tmp,,drop=FALSE]
      
      tmp <- intersect(rownames(HLA_type2tumor_nomissmatch),rownames(HLA_type2normal_nomissmatch))
      HLA_type2tumor_nomissmatch  <- HLA_type2tumor_nomissmatch[tmp,,drop=FALSE]
      
      HLA_type2normal_nomissmatchCov <- HLA_type2normal_nomissmatch$V4
      names(HLA_type2normal_nomissmatchCov) <- HLA_type2normal_nomissmatch$V2
      
      HLA_type1normal_nomissmatchCov <- HLA_type1normal_nomissmatch$V4
      names(HLA_type1normal_nomissmatchCov) <- HLA_type1normal_nomissmatch$V2 
      
      HLA_type2tumor_nomissmatchCov <- rep(0,length(HLA_type2normal_nomissmatchCov))      
      names(HLA_type2tumor_nomissmatchCov) <- names(HLA_type2normal_nomissmatchCov)
      HLA_type2tumor_nomissmatchCov[rownames(HLA_type2tumor_nomissmatch)] <- HLA_type2tumor_nomissmatch$V4
      
      HLA_type1tumor_nomissmatchCov <- rep(0,length(HLA_type1normal_nomissmatchCov))      
      names(HLA_type1tumor_nomissmatchCov) <- names(HLA_type1normal_nomissmatchCov)
      HLA_type1tumor_nomissmatchCov[rownames(HLA_type1tumor_nomissmatch)] <- HLA_type1tumor_nomissmatch$V4
      
    }
    
    if(extractUniqueReads == TRUE) {
      
      # use the same mismatch positions, already have the bed to get reads that do overlap a mismatch
      Type1TumorCmd <- paste("bedtools intersect -loj -bed -b ",workDir,region,"/",region, ".type.",HLA_A_type1,".filtered.bam"," -a ",workDir,region,".",HLA_A_type1,".bed"," > ",workDir,region,".",HLA_A_type1,".tumor.mismatch.reads.bed",sep="")
      system(Type1TumorCmd)
      if(runWithNormal){
        Type1NormalCmd <- paste("bedtools intersect -loj -bed -b ",workDir,normalName,"/",normalName, ".type.",HLA_A_type1,".filtered.bam"," -a ",workDir,region,".",HLA_A_type1,".bed"," > ",workDir,region,".",HLA_A_type1,".normal.mismatch.reads.bed",sep="")
        system(Type1NormalCmd)
      }

      Type2TumorCmd <- paste("bedtools intersect -loj -bed -b ",workDir,region,"/",region, ".type.",HLA_A_type2,".filtered.bam"," -a ",workDir,region,".",HLA_A_type2,".bed"," > ",workDir,region,".",HLA_A_type2,".tumor.mismatch.reads.bed",sep="")
      system(Type2TumorCmd)
      if(runWithNormal){
        Type2NormalCmd <- paste("bedtools intersect -loj -bed -b ",workDir,normalName,"/",normalName, ".type.",HLA_A_type2,".filtered.bam"," -a ",workDir,region,".",HLA_A_type2,".bed"," > ",workDir,region,".",HLA_A_type2,".normal.mismatch.reads.bed",sep="")
        system(Type2NormalCmd)
      }

      # only take unique reads 
      for(i in grep(pattern = HLA_gene, x = grep(pattern = region, x = list.files(workDir, pattern = 'mismatch.reads.bed', full.names = TRUE), value = TRUE), value = TRUE)){
        #print(i)
        if(file.size(i) > 0){
          x <- read.csv(i, sep = '\t', as.is = TRUE, header = FALSE)
          x <- x[!is.na(x$V2),]
          x <- x[!duplicated(x$V8),]
        } else{
          x <- t(c(paste('V', seq(1, 10), sep = '')))
        }
        write.table(x, file = gsub(pattern = 'mismatch.reads.bed', replacement = 'mismatch.unique.reads.bed', x = i), sep = '\t', quote = FALSE, row.names = FALSE)
      }
      
      HLA_A_type1normalCov_mismatch        <- HLA_A_type1normalCov[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1]
      HLA_A_type1normalCov_mismatch_unique <- HLA_A_type1normalCov_mismatch
      if(runWithNormal){
        HLA_A_type1normal_unique             <- read.table(paste(workDir, region, '.', HLA_A_type1, '.normal.mismatch.unique.reads.bed', sep = ''), sep= '\t', header = TRUE, as.is = TRUE)
        HLA_A_type1normalCov_mismatch_unique <- sapply(X = names(HLA_A_type1normalCov_mismatch_unique), FUN = function(x) {return(table(HLA_A_type1normal_unique$V3)[x])} )  
        names(HLA_A_type1normalCov_mismatch_unique) <- names(HLA_A_type1normalCov_mismatch)
        HLA_A_type1normalCov_mismatch_unique[is.na(HLA_A_type1normalCov_mismatch_unique)] <- 0
        if(length(HLA_A_type1normalCov_mismatch_unique) != 0){
          HLA_A_type1normalCov_mismatch_unique <- HLA_A_type1normalCov_mismatch_unique + 1
        }
      }

      HLA_A_type1tumorCov_mismatch        <- HLA_A_type1tumorCov[names(HLA_A_type1tumorCov)%in%missMatchPositions$diffSeq1]
      HLA_A_type1tumor_unique             <- read.table(paste(workDir, region, '.', HLA_A_type1, '.tumor.mismatch.unique.reads.bed', sep = ''), sep= '\t', header = TRUE, as.is = TRUE)
      HLA_A_type1tumorCov_mismatch_unique <- HLA_A_type1tumorCov_mismatch
      HLA_A_type1tumorCov_mismatch_unique <- sapply(X = names(HLA_A_type1tumorCov_mismatch_unique), FUN = function(x) {return(table(HLA_A_type1tumor_unique$V3)[x])} )  
      names(HLA_A_type1tumorCov_mismatch_unique) <- names(HLA_A_type1tumorCov_mismatch)
      HLA_A_type1tumorCov_mismatch_unique[is.na(HLA_A_type1tumorCov_mismatch_unique)] <- 0
      if(length(HLA_A_type1tumorCov_mismatch_unique) != 0){
        HLA_A_type1tumorCov_mismatch_unique <- HLA_A_type1tumorCov_mismatch_unique + 1
      }
      
      HLA_A_type2normalCov_mismatch        <- HLA_A_type2normalCov[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2]
      HLA_A_type2normalCov_mismatch_unique <- HLA_A_type2normalCov_mismatch
      if(runWithNormal){
        HLA_A_type2normal_unique             <- read.table(paste(workDir, region, '.', HLA_A_type2, '.normal.mismatch.unique.reads.bed', sep = ''), sep= '\t', header = TRUE, as.is = TRUE)
        HLA_A_type2normalCov_mismatch_unique <- sapply(X = names(HLA_A_type2normalCov_mismatch_unique), FUN = function(x) {return(table(HLA_A_type2normal_unique$V3)[x])} )  
        names(HLA_A_type2normalCov_mismatch_unique) <- names(HLA_A_type2normalCov_mismatch)
        HLA_A_type2normalCov_mismatch_unique[is.na(HLA_A_type2normalCov_mismatch_unique)] <- 0
        if(length(HLA_A_type2normalCov_mismatch_unique) != 0){
          HLA_A_type2normalCov_mismatch_unique <- HLA_A_type2normalCov_mismatch_unique + 1
        }
      }

      HLA_A_type2tumorCov_mismatch        <- HLA_A_type2tumorCov[names(HLA_A_type2tumorCov)%in%missMatchPositions$diffSeq2]
      HLA_A_type2tumor_unique             <- read.table(paste(workDir, region, '.', HLA_A_type2, '.tumor.mismatch.unique.reads.bed', sep = ''), sep= '\t', header = TRUE, as.is = TRUE)
      HLA_A_type2tumorCov_mismatch_unique <- HLA_A_type2tumorCov_mismatch
      HLA_A_type2tumorCov_mismatch_unique <- sapply(X = names(HLA_A_type2tumorCov_mismatch_unique), FUN = function(x) {return(table(HLA_A_type2tumor_unique$V3)[x])} )  
      names(HLA_A_type2tumorCov_mismatch_unique) <- names(HLA_A_type2tumorCov_mismatch)
      HLA_A_type2tumorCov_mismatch_unique[is.na(HLA_A_type2tumorCov_mismatch_unique)] <- 0
      if(length(HLA_A_type2tumorCov_mismatch_unique) != 0){
        HLA_A_type2tumorCov_mismatch_unique <- HLA_A_type2tumorCov_mismatch_unique + 1
      }
    
    }
    
    
    #plotting
    {
      pdf(paste(figureDir,region,".minCoverage_",minCoverageFilter,".HLA.pdf",sep=""),width=10,height=6)

      if(runWithNormal){
        # rolling mean
        par(mfrow=c(2,1))
        par(mar=c(2,5,2,2))
        barplot(c(rollmean(HLA_A_type2tumorCov*MultFactor,min(500, length(HLA_A_type2tumorCov)))/rollmean(HLA_A_type2normalCov,min(500, length(HLA_A_type2normalCov)))),ylim=c(0,3),xaxt='n',main=HLA_A_type2,las=1, ylab = 'Tumor/Normal Coverage')
        abline(h=1,lty='dashed',col='blue',lwd=1.5)
        barplot(c(rollmean(HLA_A_type1tumorCov*MultFactor,min(500, length(HLA_A_type1tumorCov)))/rollmean(HLA_A_type1normalCov,min(500, length(HLA_A_type1normalCov)))),ylim=c(0,3),xaxt='n',main=HLA_A_type1,las=1, ylab = 'Tumor/Normal Coverage')
        abline(h=1,lty='dashed',col='blue',lwd=1.5)

        
          
        # log ratio and density of mismatches
        par(mfrow=c(1,1))
        par(mar=c(5,5,5,2))
        plot(c(1:max(HLA_A_type1normal$V2,HLA_A_type2normal$V2,HLA_A_type2tumor$V2,HLA_A_type1tumor$V2))
             ,ylim=c(-3,3),col='#3182bd99',pch=16
             ,xlab='HLA genomic position'
             ,ylab='Log Ratio'
             ,main=c(paste("HLA raw balance",region))
             ,cex=0.75
             ,type='n'
             ,las=1)
        # add the exonic positions
        if(length(HLAtype1exons)!=0)
        {
          for (i in 1:nrow(HLAtype1exonTable))
          {
            rect(xleft = as.numeric(HLAtype1exonTable[i,1]),xright = as.numeric(HLAtype1exonTable[i,2]),ybottom = -2,ytop=2,col='#bdbdbd25',border=FALSE)
          }
        }
        
        if(length(HLAtype2exons)!=0)
        {
          for (i in 1:nrow(HLAtype2exonTable))
          {
            rect(xleft = as.numeric(HLAtype2exonTable[i,1]),xright = as.numeric(HLAtype2exonTable[i,2]),ybottom = -2,ytop=2,col='#bdbdbd25',border=FALSE)
          }
        }
        
        points(c(names(HLA_A_type2normalCov)),log2(c((HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)),col='#3182bd99',pch=16)
        points(names(HLA_A_type2normalCov)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2],log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2],col='black',bg='#3182bd99',pch=21,cex=1)
        points(names(HLA_A_type1normalCov),log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor),col='#de2d2699',pch=16,cex=0.75)
        points(names(HLA_A_type1normalCov)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1],log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1],col='black',bg='#de2d2699',pch=21,cex=1)
        
        points(missMatchPositions$diffSeq1,rep(-3,length(missMatchPositions$diffSeq1)),col='#63636399',bg='#63636399',pch=21,cex=1)
        d <- density(missMatchPositions$diffSeq1,bw=40)
        d$y <- (d$y/max(d$y))
        d$y <- d$y -3
        lines(d)
        abline(h=0,lty='dashed')
        
        legend('topright', legend = c(HLA_A_type2,HLA_A_type1) , 
               lty=1, col=c('#3182bd99', '#de2d2699'), bty='n', cex=1,lwd=3)



        
        # normal coverage comparison
        plot(c(1:max(HLA_A_type1normal$V2,HLA_A_type2normal$V2))
             ,ylim=c(0, max(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#3182bd99',pch=16
             ,xlab='HLA genomic position'
             ,ylab='Coverage'
             ,main=c(paste("HLA normal coverage",region))
             ,cex=0.75
             ,type='n'
             ,las=1)
        # add the exonic positions
        if(length(HLAtype1exons)!=0)
        {
          for (i in 1:nrow(HLAtype1exonTable))
          {
            rect(xleft = as.numeric(HLAtype1exonTable[i,1]),xright = as.numeric(HLAtype1exonTable[i,2]),ybottom = -2,ytop=max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#bdbdbd25',border=FALSE)
          }
        }
        
        if(length(HLAtype2exons)!=0)
        {
          for (i in 1:nrow(HLAtype2exonTable))
          {
            rect(xleft = as.numeric(HLAtype2exonTable[i,1]),xright = as.numeric(HLAtype2exonTable[i,2]),ybottom = -2,ytop=max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#bdbdbd25',border=FALSE)
          }
        }
        
        lines(c(HLA_A_type2normal$V2),c(HLA_A_type2normalCov),col='#3182bd')
        lines(c(HLA_A_type1normal$V2),c(HLA_A_type1normalCov),col='#de2d26')
        legend('topleft', legend = c(HLA_A_type2,HLA_A_type1) , 
               lty=1, col=c('#3182bd', '#de2d26'), bty='n', cex=1,lwd=3)
        
        points(c(HLA_A_type1normal$V2)[HLA_A_type1normal$V2%in%missMatchPositions$diffSeq1],c(HLA_A_type1normalCov)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1],col='#de2d26',pch=16)
        points(c(HLA_A_type2normal$V2)[HLA_A_type2normal$V2%in%missMatchPositions$diffSeq2],c(HLA_A_type2normalCov)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2],col='#3182bd',pch=16)
        
        

        # tumor and normal coverage for allele 1
        plot(names(HLA_A_type1normalCov),apply(cbind(HLA_A_type1tumorCov*MultFactor,HLA_A_type1normalCov),1,max),type='n',xlab='HLA genomic position',ylab='Coverage',las=1,main=HLA_A_type1,ylim=c(0,max(c(HLA_A_type1tumorCov,HLA_A_type1normalCov))))
        if(length(HLAtype1exons)!=0)
        {
          for (i in 1:nrow(HLAtype1exonTable))
          {
            rect(xleft = as.numeric(HLAtype1exonTable[i,1]),xright = as.numeric(HLAtype1exonTable[i,2]),ybottom = 0,ytop=max(c(HLA_A_type1tumorCov,HLA_A_type1normalCov)),col='#bdbdbd50',border=FALSE)
          }
        }
        
        lines(c(HLA_A_type1normal$V2),c(HLA_A_type1tumorCov*MultFactor),col='#3182bd')
        lines(c(HLA_A_type1normal$V2),c(HLA_A_type1normalCov),col='#9ecae1')
        legend('topleft', legend = c('tumour','normal') , 
               lty=1, col=c('#3182bd', '#9ecae1'), bty='n', cex=1,lwd=3)
        
        points(c(HLA_A_type1normal$V2)[HLA_A_type1normal$V2%in%missMatchPositions$diffSeq1],c(HLA_A_type1tumorCov*MultFactor)[names(HLA_A_type1tumorCov)%in%missMatchPositions$diffSeq1],col='#3182bd',pch=16)
        
        
        
        # tumor and normal coverage for allele 2
        plot(c(HLA_A_type2normal$V2),apply(cbind(HLA_A_type2tumorCov*MultFactor,HLA_A_type2normalCov),1,max),type='n',xlab='HLA genomic position',ylab='Coverage',las=1,main=HLA_A_type2,ylim=c(0,max(c(HLA_A_type2tumorCov,HLA_A_type2normalCov))))
        if(length(HLAtype2exons)!=0)
        {
          for (i in 1:nrow(HLAtype2exonTable))
          {
            rect(xleft = as.numeric(HLAtype2exonTable[i,1]),xright = as.numeric(HLAtype2exonTable[i,2]),ybottom = 0,ytop=max(c(HLA_A_type2tumorCov,HLA_A_type2normalCov)),col='#bdbdbd50',border=FALSE)
          }
        }
        
        lines(c(HLA_A_type2normal$V2),c(HLA_A_type2tumorCov*MultFactor),col='#3182bd')
        lines(c(HLA_A_type2normal$V2),c(HLA_A_type2normalCov),col='#9ecae1')
        legend('topleft', legend = c('tumour','normal') , 
               lty=1, col=c('#3182bd', '#9ecae1'), bty='n', cex=1,lwd=3)
        # points(c(HLA_A_type2normal$V2)[HLA_A_type2tumor$V2%in%missMatchPositions$diffSeq2],c(HLA_A_type2tumorCov)[names(HLA_A_type2tumorCov)%in%missMatchPositions$diffSeq2],col='#3182bd',pch=16)
        points(c(HLA_A_type2normal$V2)[HLA_A_type2normal$V2%in%missMatchPositions$diffSeq2],c(HLA_A_type2tumorCov*MultFactor)[names(HLA_A_type2tumorCov*MultFactor)%in%missMatchPositions$diffSeq2],col='#3182bd',pch=16)
        
        

        #let's now just plot the mismatch positions
        plot(c(1:max(HLA_A_type1normal$V2,HLA_A_type2normal$V2,HLA_A_type2tumor$V2,HLA_A_type1tumor$V2)),ylim=c(-2,2),col='#3182bd99',pch=16
             ,xlab='HLA genomic position'
             ,ylab='Log Ratio'
             ,main=c(paste("HLA raw balance",region))
             ,cex=0.75
             ,type='n')
        points(c(HLA_A_type2normal$V2)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2],log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2],col='#3182bd99',pch=16,cex=1)
        points(c(HLA_A_type1normal$V2)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1],log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1],col='#de2d2699',pch=16,cex=1)
        abline(h=0,lty='dashed')
      
      }

      if(!runWithNormal){
        par(mfrow=c(2,1))
        par(mar=c(2,5,2,2))
        Ymax <- max(c(HLA_A_type1tumorCov,HLA_A_type2tumorCov))+100
        
        barplot(c(rollmean(HLA_A_type1tumorCov,150)),ylim=c(0,Ymax),xaxt='n',main=HLA_A_type1,las=1,col='#de2d2699',border='#de2d2650')
        abline(h=median(HLA_A_type1tumorCov),lty='dashed',col='#de2d26',lwd=1.5, na.rm = TRUE)
        abline(h=median(HLA_A_type2tumorCov),lty='dashed',col='#3182bd',lwd=1.5, na.rm = TRUE)
        barplot(c(rollmean(HLA_A_type2tumorCov,150)),ylim=c(0,Ymax),xaxt='n',main=HLA_A_type2,las=1,col='#3182bd',border='#3182bd50')
        abline(h=median(HLA_A_type1tumorCov),lty='dashed',col='#de2d26',lwd=1.5, na.rm = TRUE)
        abline(h=median(HLA_A_type2tumorCov),lty='dashed',col='#3182bd',lwd=1.5, na.rm = TRUE)
        

        par(mfrow=c(1,1))
        par(mar=c(5,5,5,2))
        plot(c(1:max(c(HLA_A_type1tumor$V2,HLA_A_type2tumor$V2)))
             ,ylim=c(0,Ymax),col='#3182bd99',pch=16
             ,xlab='HLA genomic position'
             ,ylab='Coverage'
             ,main=c(paste("HLA raw coverage",region))
             ,cex=0.75
             ,type='n'
             ,las=1)
        
        points(HLA_A_type1tumor$V2,HLA_A_type1tumor$V4,col='#de2d2699',pch=16)
        points(HLA_A_type1tumor$V2[HLA_A_type1tumor$V2%in%missMatchPositions$diffSeq1],HLA_A_type1tumor$V4[HLA_A_type1tumor$V2%in%missMatchPositions$diffSeq1],col='black',bg='#de2d2699',pch=21,cex=1)
        
        points(HLA_A_type2tumor$V2,HLA_A_type2tumor$V4,col='#3182bd99',pch=16)
        points(HLA_A_type2tumor$V2[HLA_A_type2tumor$V2%in%missMatchPositions$diffSeq2],HLA_A_type2tumor$V4[HLA_A_type2tumor$V2%in%missMatchPositions$diffSeq2],col='black',bg='#3182bd99',pch=21,cex=1)
        
        legend('topleft', legend = c(HLA_A_type2,HLA_A_type1) , 
               lty=1, col=c('#3182bd99', '#de2d2699'), bty='n', cex=1,lwd=3)
        
        
        plot(c(1:max(c(HLA_A_type2tumor$V2,HLA_A_type1tumor$V2)))
             ,ylim=c(0,Ymax),col='#3182bd99',pch=16
             ,xlab='HLA genomic position'
             ,ylab='Coverage'
             ,main=c(paste("HLA raw balance",region))
             ,cex=0.75
             ,type='n'
             ,las=1)
        points(HLA_A_type1tumor$V2[HLA_A_type1tumor$V2%in%missMatchPositions$diffSeq1],HLA_A_type1tumor$V4[HLA_A_type1tumor$V2%in%missMatchPositions$diffSeq1],col='#de2d2699',pch=16,cex=1)
        points(HLA_A_type2tumor$V2[HLA_A_type2tumor$V2%in%missMatchPositions$diffSeq2],HLA_A_type2tumor$V4[HLA_A_type2tumor$V2%in%missMatchPositions$diffSeq2],col='#3182bd99',pch=16,cex=1)
        
        legend('topleft', legend = c(HLA_A_type2,HLA_A_type1) , 
               lty=1, col=c('#3182bd99', '#de2d2699'), bty='n', cex=1,lwd=3)

      }
      
      #next, let's plot the predicted integer copy numbers 
      
      if(performIntegerCopyNum)
      {
        copyNumSolutions           <- read.table(CopyNumLoc,sep="\t",header=TRUE,stringsAsFactors=FALSE)
        tumorPloidy                <- copyNumSolutions[region,'psi']
        tumorPurity                <- copyNumSolutions[region,'Cellularity']
        
           
        # infer copy number using combined BAF and logR
        
        misMatchCoveredInBoth <- cbind(ifelse(missMatchPositions$diffSeq1%in%names(HLA_A_type1normalCov),1,0),ifelse(missMatchPositions$diffSeq2%in%names(HLA_A_type2normalCov),1,0))
        missMatchseq1 <- missMatchPositions$diffSeq1[rowSums(misMatchCoveredInBoth)==2]
        missMatchseq2 <- missMatchPositions$diffSeq2[rowSums(misMatchCoveredInBoth)==2]

        #let's get bins to collect for the coverage estimates
        startChar <- min(c(as.numeric(names(HLA_A_type1tumorCov))),as.numeric(names(HLA_A_type2tumorCov)))
        endChar   <- max(c(as.numeric(names(HLA_A_type1tumorCov))),as.numeric(names(HLA_A_type2tumorCov)))
        seqToConsider <- seq(startChar,endChar,by=binSize)
        seqToConsider <- c(seqToConsider[-length(seqToConsider)],endChar)
        
        binLogR       <- c()
        for (i in 1:(length(seqToConsider)-1))
        {
          
          PotentialSites   <- as.character(seqToConsider[i]:seqToConsider[i+1])
          combinedBinTumor  <- median(as.numeric(as.numeric(HLA_A_type1tumorCov[names(HLA_A_type1tumorCov)%in%PotentialSites])), na.rm = TRUE)+median(as.numeric(as.numeric(HLA_A_type2tumorCov[names(HLA_A_type2tumorCov)%in%PotentialSites])), na.rm = TRUE)
          combinedBinNormal <- median(as.numeric(as.numeric(HLA_A_type1normalCov[names(HLA_A_type1normalCov)%in%PotentialSites])), na.rm= TRUE)+median(as.numeric(as.numeric(HLA_A_type2normalCov[names(HLA_A_type2normalCov)%in%PotentialSites])), na.rm = TRUE)
          combinedBinlogR   <- log2(combinedBinTumor/combinedBinNormal*MultFactor)
          type1BinlogR     <- median(log2(as.numeric(as.numeric(HLA_A_type1tumorCov[names(HLA_A_type1tumorCov)%in%PotentialSites])/as.numeric(HLA_A_type1normalCov[names(HLA_A_type1normalCov)%in%PotentialSites])*MultFactor)), na.rm = TRUE)
          type2BinlogR     <- median(log2(as.numeric(as.numeric(HLA_A_type2tumorCov[names(HLA_A_type2tumorCov)%in%PotentialSites])/as.numeric(HLA_A_type2normalCov[names(HLA_A_type2normalCov)%in%PotentialSites])*MultFactor)), na.rm = TRUE)
          binLogR <- rbind(binLogR,cbind(seqToConsider[i],seqToConsider[i+1],combinedBinlogR,type1BinlogR,type2BinlogR))
        }


        tmpOut <- cbind(missMatchseq1
                        ,log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[as.character(missMatchseq1)]
                        ,HLA_A_type1tumorCov[as.character(missMatchseq1)]
                        ,missMatchseq2
                       ,log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[as.character(missMatchseq2)]
                       ,HLA_A_type2tumorCov[as.character(missMatchseq2)]
                       ,HLA_A_type1normalCov[as.character(missMatchseq1)]
                       ,HLA_A_type2normalCov[as.character(missMatchseq2)]
                       )
        colnames(tmpOut) <- c('missMatchseq1','logR_type1','TumorCov_type1','missMatchseq2','logR_type2','TumorCov_type2','NormalCov_type1','NormalCov_type2')
      
        
        dup1   <- unique(tmpOut[duplicated(tmpOut[,1]),1])
        dup2   <- unique(tmpOut[duplicated(tmpOut[,4]),4])
        
        for(duplicationIn1 in dup1)
        {
          tmpOut[tmpOut[,1]==duplicationIn1,'TumorCov_type2']  <- mean(tmpOut[tmpOut[,1]==duplicationIn1,'TumorCov_type2'])
          tmpOut[tmpOut[,1]==duplicationIn1,'NormalCov_type2'] <- mean(tmpOut[tmpOut[,1]==duplicationIn1,'NormalCov_type2'])  
          tmpOut[tmpOut[,1]==duplicationIn1,'logR_type2']      <- mean(tmpOut[tmpOut[,1]==duplicationIn1,'logR_type2'])  
        }
        
        for(duplicationIn2 in dup2)
        {
          
          tmpOut[tmpOut[,4]==duplicationIn2,'TumorCov_type1']  <- mean(tmpOut[tmpOut[,4]==duplicationIn2,'TumorCov_type1'])
          tmpOut[tmpOut[,4]==duplicationIn2,'NormalCov_type1'] <- mean(tmpOut[tmpOut[,4]==duplicationIn2,'NormalCov_type1'])  
          tmpOut[tmpOut[,4]==duplicationIn2,'logR_type1']      <- mean(tmpOut[tmpOut[,4]==duplicationIn2,'logR_type1'])  
          
        }
        
        tmpOut  <- tmpOut[!duplicated(tmpOut[,1]),,drop=FALSE]
        tmpOut  <- tmpOut[!duplicated(tmpOut[,4]),,drop=FALSE]
        colnames(tmpOut) <- c('missMatchseq1','logR_type1','TumorCov_type1','missMatchseq2','logR_type2','TumorCov_type2','NormalCov_type1','NormalCov_type2')


        combinedTable <- data.frame(tmpOut,stringsAsFactors=FALSE)
        combinedTable$logRcombined <- log2(((combinedTable$TumorCov_type1+combinedTable$TumorCov_type2)/(combinedTable$NormalCov_type1+combinedTable$NormalCov_type2))*MultFactor)
        combinedTable$BAFcombined  <- combinedTable$TumorCov_type1/(combinedTable$TumorCov_type1+combinedTable$TumorCov_type2)
        

        if(nrow(combinedTable) != 0){
          combinedTable$binlogRCombined <- NA
          combinedTable$binlogRtype1    <- NA
          combinedTable$binlogRtype2    <- NA
          combinedTable$binNum          <- NA

          # next, add the binLogR to this table
          for (i in 1:nrow(combinedTable))
          {
            combinedTable[i,]$binlogRCombined <- binLogR[which(binLogR[,1]<=as.numeric(combinedTable$missMatchseq1[i])&binLogR[,2]>as.numeric(combinedTable$missMatchseq1[i])),3]
            combinedTable[i,]$binlogRtype1   <- binLogR[which(binLogR[,1]<=as.numeric(combinedTable$missMatchseq1[i])&binLogR[,2]>as.numeric(combinedTable$missMatchseq1[i])),4]
            combinedTable[i,]$binlogRtype2 <- binLogR[which(binLogR[,1]<=as.numeric(combinedTable$missMatchseq1[i])&binLogR[,2]>as.numeric(combinedTable$missMatchseq1[i])),5]
            combinedTable[i,]$binNum       <- binLogR[which(binLogR[,1]<=as.numeric(combinedTable$missMatchseq1[i])&binLogR[,2]>as.numeric(combinedTable$missMatchseq1[i])),1]
          }
        }

        if(nrow(combinedTable) == 0){
          combinedTable$binlogRCombined <- numeric(0)
          combinedTable$binlogRtype1    <- numeric(0)
          combinedTable$binlogRtype2    <- numeric(0)
          combinedTable$binNum          <- numeric(0)
        }
        
        rawVals <- funCalcN_withBAF(combinedTable$logRcombined,combinedTable$BAFcombined,tumorPloidy,tumorPurity,gamma)
        combinedTable$nAcombined <- rawVals[,1]
        combinedTable$nBcombined <- rawVals[,2]

        nB_rawVal_withBAF        <- median(combinedTable$nBcombined, na.rm = TRUE)
        nB_rawVal_withBAF_conf   <- t.test.NA(combinedTable$nBcombined)
        nB_rawVal_withBAF_lower  <- nB_rawVal_withBAF_conf$conf.int[1]
        nB_rawVal_withBAF_upper  <- nB_rawVal_withBAF_conf$conf.int[2]
        
        nA_rawVal_withBAF        <- median(combinedTable$nAcombined, na.rm = TRUE)
        nA_rawVal_withBAF_conf   <- t.test.NA(combinedTable$nAcombined)
        nA_rawVal_withBAF_lower  <- nA_rawVal_withBAF_conf$conf.int[1]
        nA_rawVal_withBAF_upper  <- nA_rawVal_withBAF_conf$conf.int[2]
        
        
        
        rawValsBin      <- funCalcN_withBAF(combinedTable$binlogRCombined,combinedTable$BAFcombined,tumorPloidy,tumorPurity,gamma)
        combinedTable$nAcombinedBin <- rawValsBin[,1]
        combinedTable$nBcombinedBin <- rawValsBin[,2]
        
        nB_rawVal_withBAF           <- median(combinedTable$nBcombined, na.rm = TRUE)
        nB_rawVal_withBAF_conf      <- t.test.NA(combinedTable$nBcombined)
        nB_rawVal_withBAF_lower     <- nB_rawVal_withBAF_conf$conf.int[1]
        nB_rawVal_withBAF_upper     <- nB_rawVal_withBAF_conf$conf.int[2]
        
        nA_rawVal_withBAF           <- median(combinedTable$nAcombined, na.rm = TRUE)
        nA_rawVal_withBAF_conf      <- t.test.NA(combinedTable$nAcombined)
        nA_rawVal_withBAF_lower     <- nA_rawVal_withBAF_conf$conf.int[1]
        nA_rawVal_withBAF_upper     <- nA_rawVal_withBAF_conf$conf.int[2]
        
        #let's only count non duplicates
        nB_rawVal_withBAF_bin <- median(combinedTable[!duplicated(combinedTable$binNum),]$nBcombinedBin, na.rm = TRUE)
        nB_rawVal_withBAF_bin_conf <- t.test.NA(combinedTable[!duplicated(combinedTable$binNum),]$nBcombinedBin)
        nB_rawVal_withBAF_bin_lower <-nB_rawVal_withBAF_bin_conf$conf.int[1]
        nB_rawVal_withBAF_bin_upper <-nB_rawVal_withBAF_bin_conf$conf.int[2]
        
        nA_rawVal_withBAF_bin <- median(combinedTable[!duplicated(combinedTable$binNum),]$nAcombinedBin, na.rm = TRUE)
        nA_rawVal_withBAF_bin_conf <- t.test.NA(combinedTable[!duplicated(combinedTable$binNum),]$nAcombinedBin)
        nA_rawVal_withBAF_bin_lower <- nA_rawVal_withBAF_bin_conf$conf.int[1]
        nA_rawVal_withBAF_bin_upper <- nA_rawVal_withBAF_bin_conf$conf.int[2]

        if(!is.na(tumorPurity) & nrow(combinedTable) != 0){
          
          plot(c(1:max(HLA_A_type1normal$V2,HLA_A_type2normal$V2,HLA_A_type2tumor$V2,HLA_A_type1tumor$V2)),ylim=c(-0.5,max(round(c(combinedTable$nBcombined,combinedTable$nAcombined)), na.rm = TRUE)+1),col='#3182bd99',pch=16
               ,xlab='HLA genomic position'
               ,ylab='Copy Number'
               ,main=c(paste("HLA copyNum balance",region))
               ,cex=0.75
               ,type='n'
               ,yaxt='n')
          
          if(length(HLAtype1exons)!=0)
          {
            for (i in 1:nrow(HLAtype1exonTable))
            {
              rect(xleft = as.numeric(HLAtype1exonTable[i,1]),xright = as.numeric(HLAtype1exonTable[i,2]),ybottom = -2,ytop=max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#bdbdbd25',border=FALSE)
            }
          }
          
          if(length(HLAtype2exons)!=0)
          {
            for (i in 1:nrow(HLAtype2exonTable))
            {
              rect(xleft = as.numeric(HLAtype2exonTable[i,1]),xright = as.numeric(HLAtype2exonTable[i,2]),ybottom = -2,ytop=max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#bdbdbd25',border=FALSE)
            }
          }
          
          if(useLogRbin)
          {
            axis(side = 2,at = 0:max(round(c(combinedTable$nAcombinedBin,combinedTable$nBcombinedBin)), na.rm = TRUE),las=1)
            points(combinedTable$missMatchseq1,combinedTable$nAcombinedBin,col='#de2d2699',pch=16,cex=1)
            points(combinedTable$missMatchseq2,combinedTable$nBcombinedBin,col='#3182bd99',pch=16,cex=1)
            
            abline(h=nA_rawVal_withBAF_bin,lty='dashed',col='#de2d2699',lwd=3)
            abline(h=nB_rawVal_withBAF_bin,lty='dashed',col='#3182bd99',lwd=3)
            
            legend('topleft', legend = c(HLA_A_type2,HLA_A_type1) , 
                   lty=1, col=c('#3182bd', '#de2d26'), bty='n', cex=1,lwd=3)
          }
          
          if(!useLogRbin)
          {
            axis(side = 2,at = 0:max(round(c(combinedTable$nAcombined,combinedTable$nBcombined)), na.rm = TRUE),las=1)
            points(combinedTable$missMatchseq1,combinedTable$nAcombined,col='#de2d2699',pch=16,cex=1)
            points(combinedTable$missMatchseq2,combinedTable$nBcombined,col='#3182bd99',pch=16,cex=1)
            
            abline(h=nA_rawVal_withBAF,lty='dashed',col='#de2d2699',lwd=3)
            abline(h=nB_rawVal_withBAF,lty='dashed',col='#3182bd99',lwd=3)
            
            legend('topleft', legend = c(HLA_A_type2,HLA_A_type1) , 
                   lty=1, col=c('#3182bd', '#de2d26'), bty='n', cex=1,lwd=3)
          }

        }
        
        combinedTable$nAsep    <- funCalcN_withoutBAF(combinedTable$logR_type1,tumorPloidy,tumorPurity,gamma)
        combinedTable$nAsepBin <- funCalcN_withoutBAF(combinedTable$binlogRtype1,tumorPloidy,tumorPurity,gamma)
        combinedTable$nBsep    <- funCalcN_withoutBAF(combinedTable$logR_type2,tumorPloidy,tumorPurity,gamma)
        combinedTable$nBsepBin <- funCalcN_withoutBAF(combinedTable$binlogRtype2,tumorPloidy,tumorPurity,gamma)


        nB_rawVal_withoutBAF <- median(combinedTable$nBsep, na.rm = TRUE)
        nB_rawVal_withoutBAF_conf <- t.test.NA(combinedTable$nBsep)
        nB_rawVal_withoutBAF_lower <- nB_rawVal_withoutBAF_conf$conf.int[1]
        nB_rawVal_withoutBAF_upper <- nB_rawVal_withoutBAF_conf$conf.int[2]
    
        nA_rawVal_withoutBAF <- median(combinedTable$nAsep, na.rm = TRUE)

        nA_rawVal_withoutBAF_conf <- t.test.NA(combinedTable$nAsep)
        nA_rawVal_withoutBAF_lower <- nA_rawVal_withoutBAF_conf$conf.int[1]
        nA_rawVal_withoutBAF_upper <- nA_rawVal_withoutBAF_conf$conf.int[2]

        #median(combinedTable[!duplicated(combinedTable$binNum),]$nAcombinedBin)
        
        nB_rawVal_withoutBAFBin <-  median(combinedTable[!duplicated(combinedTable$binNum),]$nBsepBin, na.rm = TRUE)
        nB_rawVal_withoutBAFBin_conf <-  t.test.NA(combinedTable[!duplicated(combinedTable$binNum),]$nBsepBin)
        nB_rawVal_withoutBAFBin_lower <- nB_rawVal_withoutBAFBin_conf$conf.int[1]
        nB_rawVal_withoutBAFBin_upper <- nB_rawVal_withoutBAFBin_conf$conf.int[2]
        
        
        nA_rawVal_withoutBAFBin <-  median(combinedTable[!duplicated(combinedTable$binNum),]$nAsepBin, na.rm = TRUE)
        nA_rawVal_withoutBAFBin_conf <-  t.test.NA(combinedTable[!duplicated(combinedTable$binNum),]$nAsepBin)
        nA_rawVal_withoutBAFBin_lower <- nA_rawVal_withoutBAFBin_conf$conf.int[1]
        nA_rawVal_withoutBAFBin_upper <- nA_rawVal_withoutBAFBin_conf$conf.int[2]
        
        if(!is.na(tumorPurity) & nrow(combinedTable) != 0){

          if(!useLogRbin)
          {
            plot(c(1:max(HLA_A_type1normal$V2,HLA_A_type2normal$V2,HLA_A_type2tumor$V2,HLA_A_type1tumor$V2)),ylim=c(-0.5,max(round(c(combinedTable$nAsep,combinedTable$nBsep)), na.rm = TRUE)+1),col='#3182bd99',pch=16
                 ,xlab='HLA genomic position'
                 ,ylab='Copy Number'
                 ,main=c(paste("HLA copyNum balance",region))
                 ,cex=0.75
                 ,type='n'
                 ,yaxt='n')
            axis(side = 2,at = 0:max(round(c(combinedTable$nAsep,combinedTable$nBsep))),las=1)
            
            if(length(HLAtype1exons)!=0)
            {
              for (i in 1:nrow(HLAtype1exonTable))
              {
                rect(xleft = as.numeric(HLAtype1exonTable[i,1]),xright = as.numeric(HLAtype1exonTable[i,2]),ybottom = -2,ytop=max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#bdbdbd25',border=FALSE)
              }
            }
            
            if(length(HLAtype2exons)!=0)
            {
              for (i in 1:nrow(HLAtype2exonTable))
              {
                rect(xleft = as.numeric(HLAtype2exonTable[i,1]),xright = as.numeric(HLAtype2exonTable[i,2]),ybottom = -2,ytop=max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#bdbdbd25',border=FALSE)
              }
            }
            
            points(combinedTable$missMatchseq1,combinedTable$nAsep,col='#de2d2699',pch=16,cex=1)
            points(combinedTable$missMatchseq2,combinedTable$nBsep,col='#3182bd99',pch=16,cex=1)
            
            abline(h=nA_rawVal_withoutBAF,lty='dashed',col='#de2d2699',lwd=3)
            abline(h=nB_rawVal_withoutBAF,lty='dashed',col='#3182bd99',lwd=3)
            
            legend('topleft', legend = c(HLA_A_type2,HLA_A_type1) , 
                   lty=1, col=c('#3182bd', '#de2d26'), bty='n', cex=1,lwd=3)
            
          }
          
          if(useLogRbin)
          {
            plot(c(1:max(HLA_A_type1normal$V2,HLA_A_type2normal$V2,HLA_A_type2tumor$V2,HLA_A_type1tumor$V2)),ylim=c(-0.5,max(round(c(combinedTable$nAsep,combinedTable$nBsep)), na.rm = TRUE)+1),col='#3182bd99',pch=16
                 ,xlab='HLA genomic position'
                 ,ylab='Copy Number'
                 ,main=c(paste("HLA copyNum balance",region))
                 ,cex=0.75
                 ,type='n'
                 ,yaxt='n')
            axis(side = 2,at = 0:max(round(c(combinedTable$nAsep,combinedTable$nBsep))),las=1)
            
            if(length(HLAtype1exons)!=0)
            {
              for (i in 1:nrow(HLAtype1exonTable))
              {
                rect(xleft = as.numeric(HLAtype1exonTable[i,1]),xright = as.numeric(HLAtype1exonTable[i,2]),ybottom = -2,ytop=max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#bdbdbd25',border=FALSE)
              }
            }
            
            if(length(HLAtype2exons)!=0)
            {
              for (i in 1:nrow(HLAtype2exonTable))
              {
                rect(xleft = as.numeric(HLAtype2exonTable[i,1]),xright = as.numeric(HLAtype2exonTable[i,2]),ybottom = -2,ytop=max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),col='#bdbdbd25',border=FALSE)
              }
            }
            
            points(combinedTable$missMatchseq1,combinedTable$nAsepBin,col='#de2d2699',pch=16,cex=1)
            points(combinedTable$missMatchseq2,combinedTable$nBsepBin,col='#3182bd99',pch=16,cex=1)
            
            abline(h=nA_rawVal_withoutBAFBin,lty='dashed',col='#de2d2699',lwd=3)
            abline(h=nB_rawVal_withoutBAFBin,lty='dashed',col='#3182bd99',lwd=3)
            
            legend('topleft', legend = c(HLA_A_type2,HLA_A_type1) , 
                   lty=1, col=c('#3182bd', '#de2d26'), bty='n', cex=1,lwd=3)
            
          }

        }
        
        # we can also predict the BAF from our logR
        combinedTable$expectedBAF <- (1-tumorPurity+tumorPurity*combinedTable$nAsep)/(2-2*tumorPurity+tumorPurity*(combinedTable$nAsep+combinedTable$nBsep))  
  
      }
      
      
      
      # t-test plot
      #let's also perform a t-test looking at the mismatch positions
      par(mfrow=c(1,2))
      par(mar=c(5,5,5,2))
      
      misMatchCoveredInBoth <- cbind(ifelse(missMatchPositions$diffSeq1%in%names(HLA_A_type1normalCov),1,0),ifelse(missMatchPositions$diffSeq2%in%names(HLA_A_type2normalCov),1,0))
      missMatchseq1 <- missMatchPositions$diffSeq1[rowSums(misMatchCoveredInBoth)==2]
      missMatchseq2 <- missMatchPositions$diffSeq2[rowSums(misMatchCoveredInBoth)==2]
      
      # we don't want to count the same miss-match multiple times - let's take an average where this is the case
      tmpOut <- cbind(missMatchseq1,log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov+0.0001)*MultFactor)[as.character(missMatchseq1)],missMatchseq2,log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov+0.0001)*MultFactor)[as.character(missMatchseq2)])
      dup1   <- unique(tmpOut[duplicated(tmpOut[,1]),1])
      dup2   <- unique(tmpOut[duplicated(tmpOut[,3]),3])
      
      for(duplicationIn1 in dup1)
      {
        tmpOut[tmpOut[,1]==duplicationIn1,4] <- median(tmpOut[tmpOut[,1]==duplicationIn1,4])
      }
      
      for(duplicationIn2 in dup2)
      {
        tmpOut[tmpOut[,3]==duplicationIn2,2] <- median(tmpOut[tmpOut[,3]==duplicationIn2,2])  
      }
      
      tmpOut  <- tmpOut[!duplicated(tmpOut[,1]),,drop=FALSE]
      tmpOut  <- tmpOut[!duplicated(tmpOut[,3]),,drop=FALSE]
      
      if(nrow(tmpOut) > 1){
        PairedTtest <- t.test(tmpOut[,2],tmpOut[,4],paired=TRUE)
        UnPairedTtest <- t.test(tmpOut[,2],tmpOut[,4],paired=FALSE)
      } else{
        PairedTtest <- list(p.value = NA)
        UnPairedTtest <- list(p.value = NA)
      }
      
      if(nrow(tmpOut) > 0){
        if(runWithNormal){
          boxplot(tmpOut[,2],tmpOut[,4],col=c('#de2d2699','#3182bd99'),boxwex=0.2,ylim=c(-2,2),names=c(HLA_A_type1,HLA_A_type2),las=1,main=paste('Paired t.test p.val=',signif(PairedTtest$p.value,3)),ylab=('logR ratio'))
        }
        if(!runWithNormal){
          boxplot(tmpOut[,2],tmpOut[,4],col=c('#de2d2699','#3182bd99'),boxwex=0.2,names=c(HLA_A_type1,HLA_A_type2),las=1,main=paste('Paired t.test p.val=',signif(PairedTtest$p.value,3)),ylab=('Coverage'))
        }
        beeswarm(tmpOut[,2],col=c('#de2d2699'),add=TRUE,corral = 'wrap',method='swarm',corralWidth=0.25,pch=16,at=1,cex=1.75)
        beeswarm(tmpOut[,4],col=c('#3182bd99'),add=TRUE,corral = 'wrap',method='swarm',corralWidth=0.25,pch=16,at=2,cex=1.75)
        barplot(sort(c(tmpOut[,2]-tmpOut[,4])),horiz = TRUE,yaxt='n',col=ifelse(sort(c(tmpOut[,2]-tmpOut[,4]))>0,'#de2d26','#3182bd'),border=FALSE,main='Paired Differences in logR')
      }
      
      
      
      # t-test of mismatch sites without counting the same read twice
      # re-create tmpOut
      
      if(!any(c(length(HLA_A_type1tumorCov_mismatch_unique), length(HLA_A_type2tumorCov_mismatch_unique)) == 0)){
        
        tmpOut_unique <- cbind(missMatchseq1,log2(c(HLA_A_type1tumorCov_mismatch_unique/HLA_A_type1normalCov_mismatch_unique + 0.0001)*MultFactor)[as.character(missMatchseq1)],missMatchseq2,log2(c(HLA_A_type2tumorCov_mismatch_unique/HLA_A_type2normalCov_mismatch_unique+0.0001)*MultFactor)[as.character(missMatchseq2)])
        dup1_unique   <- unique(tmpOut_unique[duplicated(tmpOut_unique[,1]),1])
        dup2_unique   <- unique(tmpOut_unique[duplicated(tmpOut_unique[,3]),3])
        
        for(duplicationIn1 in dup1_unique)
        {
          tmpOut_unique[tmpOut_unique[,1]==duplicationIn1,4] <- median(tmpOut_unique[tmpOut_unique[,1]==duplicationIn1,4])
          
        }
        
        for(duplicationIn2 in dup2_unique)
        {
          tmpOut_unique[tmpOut_unique[,3]==duplicationIn2,2] <- median(tmpOut_unique[tmpOut_unique[,3]==duplicationIn2,2])
          
        }
        
        tmpOut_unique  <- tmpOut_unique[!duplicated(tmpOut_unique[,1]),,drop=FALSE]
        tmpOut_unique  <- tmpOut_unique[!duplicated(tmpOut_unique[,3]),,drop=FALSE]
        
        if(nrow(tmpOut_unique) > 1){
          PairedTtest_unique <- t.test(tmpOut_unique[,2],tmpOut_unique[,4],paired=TRUE)
          UnPairedTtest_unique <- t.test(tmpOut_unique[,2],tmpOut_unique[,4],paired=FALSE)
        } else{
          PairedTtest_unique <- list(p.value = NA)
          UnPairedTtest_unique <- list(p.value = NA)
        }
        
        if(nrow(tmpOut_unique) > 0){
          if(runWithNormal){
            boxplot(tmpOut_unique[,2],tmpOut_unique[,4],col=c('#de2d2699','#3182bd99'),boxwex=0.2,ylim=c(-10,10),names=c(HLA_A_type1,HLA_A_type2),las=1,main=paste('Paired t.test p.val=',signif(PairedTtest_unique$p.value,3)),ylab=('"logR ratio"'))
          }
          if(!runWithNormal){
            boxplot(tmpOut_unique[,2],tmpOut_unique[,4],col=c('#de2d2699','#3182bd99'),boxwex=0.2,names=c(HLA_A_type1,HLA_A_type2),las=1,main=paste('Paired t.test p.val=',signif(PairedTtest_unique$p.value,3)),ylab=('"Coverage"'))       
          }
          beeswarm(tmpOut_unique[,2],col=c('#de2d2699'),add=TRUE,corral = 'wrap',method='swarm',corralWidth=0.25,pch=16,at=1,cex=1.75)
          beeswarm(tmpOut_unique[,4],col=c('#3182bd99'),add=TRUE,corral = 'wrap',method='swarm',corralWidth=0.25,pch=16,at=2,cex=1.75)
          barplot(sort(c(tmpOut_unique[,2]-tmpOut_unique[,4])),horiz = TRUE,yaxt='n',col=ifelse(sort(c(tmpOut_unique[,2]-tmpOut_unique[,4]))>0,'#de2d26','#3182bd'),border=FALSE,main='Paired Differences in logR')
        }
        
      } else{
        PairedTtest_unique <- list(p.value = NA)
        UnPairedTtest_unique <- list(p.value = NA)
      }



      # why is the p-val unique weird sometimes?
      if(!any(c(length(HLA_A_type1tumorCov_mismatch_unique), length(HLA_A_type2tumorCov_mismatch_unique)) == 0)){
        par(mfrow = c(2,2))
        tmpOut2        <- cbind(missMatchseq1,MultFactor*HLA_A_type1tumorCov[as.character(missMatchseq1)],HLA_A_type1normalCov[as.character(missMatchseq1)],missMatchseq2,MultFactor*HLA_A_type2tumorCov[as.character(missMatchseq2)], HLA_A_type2normalCov[as.character(missMatchseq2)])
        tmpOut_unique2 <- cbind(missMatchseq1,MultFactor*HLA_A_type1tumorCov_mismatch_unique[as.character(missMatchseq1)],HLA_A_type1normalCov_mismatch_unique[as.character(missMatchseq1)],missMatchseq2,MultFactor*HLA_A_type2tumorCov_mismatch_unique[as.character(missMatchseq2)], HLA_A_type2normalCov_mismatch_unique[as.character(missMatchseq2)])
        dup1_unique   <- unique(tmpOut_unique2[duplicated(tmpOut_unique2[,1]),1])
        dup2_unique   <- unique(tmpOut_unique2[duplicated(tmpOut_unique2[,4]),4])
        
        for(duplicationIn1 in dup1_unique)
        {
          tmpOut2[tmpOut2[,1]==duplicationIn1,5] <- median(tmpOut2[tmpOut2[,1]==duplicationIn1,5])
          tmpOut2[tmpOut2[,1]==duplicationIn1,6] <- median(tmpOut2[tmpOut2[,1]==duplicationIn1,6])

          tmpOut_unique2[tmpOut_unique2[,1]==duplicationIn1,5] <- median(tmpOut_unique2[tmpOut_unique2[,1]==duplicationIn1,5])
          tmpOut_unique2[tmpOut_unique2[,1]==duplicationIn1,6] <- median(tmpOut_unique2[tmpOut_unique2[,1]==duplicationIn1,6])
        }
        
        for(duplicationIn2 in dup2_unique)
        {
          tmpOut2[tmpOut2[,4]==duplicationIn1,2] <- median(tmpOut2[tmpOut2[,4]==duplicationIn1,2])
          tmpOut2[tmpOut2[,4]==duplicationIn1,3] <- median(tmpOut2[tmpOut2[,4]==duplicationIn1,3])

          tmpOut_unique2[tmpOut_unique2[,4]==duplicationIn2,2] <- median(tmpOut_unique2[tmpOut_unique2[,4]==duplicationIn2,2])
          tmpOut_unique2[tmpOut_unique2[,4]==duplicationIn2,3] <- median(tmpOut_unique2[tmpOut_unique2[,4]==duplicationIn2,3])
        }
    
        tmpOut2  <- tmpOut2[!duplicated(tmpOut2[,1]),,drop=FALSE]
        tmpOut2  <- tmpOut2[!duplicated(tmpOut2[,4]),,drop=FALSE]    
        tmpOut_unique2  <- tmpOut_unique2[!duplicated(tmpOut_unique2[,1]),,drop=FALSE]
        tmpOut_unique2  <- tmpOut_unique2[!duplicated(tmpOut_unique2[,4]),,drop=FALSE]

        if(nrow(tmpOut2) > 0){
          plot(tmpOut2[,3], tmpOut2[,2], xlab = 'normal coverage', ylab = 'tumor coverage', pch = 16, col = '#de2d2699', main = HLA_A_type1)
          abline(a = 0, b = 1)
          plot(tmpOut2[,6], tmpOut2[,5], xlab = 'normal coverage', ylab = 'tumor coverage', pch = 16, col = '#3182bd99', main = HLA_A_type2)
          abline(a = 0, b = 1)

        }
        if(nrow(tmpOut_unique2) > 0){
          plot(tmpOut_unique2[,3], tmpOut_unique2[,2], xlab = 'normal unique coverage', ylab = 'tumor unique coverage', pch = 16, col = '#de2d2699', main = HLA_A_type1)
          abline(a = 0, b = 1)
          plot(tmpOut_unique2[,6], tmpOut_unique2[,5], xlab = 'normal unique coverage', ylab = 'tumor unique coverage', pch = 16, col = '#3182bd99', main = HLA_A_type2)
          abline(a = 0, b = 1)
        
        }
      }


      if(runWithNormal){
        # rolling mean log tumor/normal
        par(mfrow=c(1,1))
        par(mar=c(5,5,5,2))
        plot(c(log2(rollmean(HLA_A_type2tumorCov,min(500, length(HLA_A_type2tumorCov)))/rollmean(HLA_A_type2normalCov,min(500, length(HLA_A_type2normalCov)))*MultFactor)),ylim=c(-2,2),col='#3182bd',pch=16
             ,xlab='HLA genomic position'
             ,ylab='Log Ratio'
             ,main=c(paste("HLA rolling mean balance",region))
             ,cex=0.75)
        points(c(log2(rollmean(HLA_A_type1tumorCov,min(500, length(HLA_A_type1tumorCov)))/rollmean(HLA_A_type1normalCov,min(500, length(HLA_A_type1normalCov)))*MultFactor)),col='#de2d26',pch=16,cex=0.75)
        legend('bottomright', legend = c(HLA_A_type2,HLA_A_type1) , 
               lty=1, col=c('#3182bd99', '#de2d2699'), bty='n', cex=1,lwd=3)
        
        
        
        # 
        if(extractNONmismatchReads==TRUE)
        {
          plot(c(1:max(HLA_A_type1normal$V2,HLA_A_type2normal$V2,HLA_A_type2tumor$V2,HLA_A_type1tumor$V2)),ylim=c(-2,2),col='#3182bd99',pch=16
               ,xlab='HLA genomic position'
               ,ylab='Log Ratio'
               ,main=c(paste("HLA raw balance",region))
               ,cex=0.75
               ,type='n')
          points(c(HLA_A_type2normal$V2)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2],log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2],col='black',bg='#3182bd',pch=21,cex=1)
          points(c(HLA_A_type1normal$V2)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1],log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1],col='black',bg='#de2d26',pch=21,cex=1)
          abline(h=0,lty='dashed')
          
          points(c(HLA_type2normal_nomissmatch$V2),log2(c(HLA_type2tumor_nomissmatchCov/HLA_type2normal_nomissmatchCov)*MultFactor),col='#3182bd99',pch=16,cex=1)
          points(c(HLA_type1normal_nomissmatch$V2),log2(c(HLA_type1tumor_nomissmatchCov/HLA_type1normal_nomissmatchCov)*MultFactor),col='#de2d2699',pch=16,cex=1)
          
          
          
          # comparing mismatch sites to non-mismatch sites
          par(mfrow=c(1,2))
          par(mar=c(5,5,5,2))
          
          if(length(HLA_type2tumor_nomissmatchCov) > 1 & length(HLA_A_type2normalCov[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2]) > 1){
            Ttest <- t.test(log(c((HLA_type2tumor_nomissmatchCov+0.01)/HLA_type2normal_nomissmatchCov)*MultFactor,2),log(c((HLA_A_type2tumorCov+0.01)/HLA_A_type2normalCov)*MultFactor,2)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2])
            boxplot(log(c((HLA_type2tumor_nomissmatchCov+0.01)/HLA_type2normal_nomissmatchCov)*MultFactor,2),log(c((HLA_A_type2tumorCov+0.01)/HLA_A_type2normalCov)*MultFactor,2)[names(HLA_A_type2normalCov)%in%missMatchPositions$diffSeq2],col=c('#de2d2699','#3182bd99'),boxwex=0.2,ylim=c(-2,2),names=c('No mismatch reads','Mismatch reads'),las=1,main=paste(HLA_A_type2,'\n t.test p.val=',signif(Ttest$p.value,3)),ylab=('logR ratio'),xlab='coverage from')
          }
          
          if(length(HLA_type1tumor_nomissmatchCov) > 1 & length(HLA_A_type1normalCov[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1]) > 1){
            Ttest <- t.test(log(c((HLA_type1tumor_nomissmatchCov+0.01)/HLA_type1normal_nomissmatchCov)*MultFactor,2),log(c((HLA_A_type1tumorCov+0.01)/HLA_A_type1normalCov)*MultFactor,2)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1])
            boxplot(log(c((HLA_type1tumor_nomissmatchCov+0.01)/HLA_type1normal_nomissmatchCov)*MultFactor,2),log(c((HLA_A_type1tumorCov+0.01)/HLA_A_type1normalCov)*MultFactor,2)[names(HLA_A_type1normalCov)%in%missMatchPositions$diffSeq1],col=c('#de2d2699','#3182bd99'),boxwex=0.2,ylim=c(-2,2),names=c('No mismatch reads','Mismatch reads'),las=1,main=paste(HLA_A_type1,'\n t.test p.val=',signif(Ttest$p.value,3)),ylab=('logR ratio'),xlab='coverage from')
          }
          
        }
      }
      
    }


    
    #let's put togehter the output
    HLAtype1Log2MedianCoverage <- median(log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor), na.rm = TRUE)
    HLAtype2Log2MedianCoverage <- median(log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor), na.rm = TRUE)
    HLAtype1Log2MedianCoverageAtSites <- median(log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[names(HLA_A_type1tumorCov)%in%missMatchseq1], na.rm = TRUE)
    HLAtype2Log2MedianCoverageAtSites <- median(log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[names(HLA_A_type2tumorCov)%in%missMatchseq2], na.rm = TRUE)
    PVal                              <- PairedTtest$p.value
    UnPairedPval                      <- UnPairedTtest$p.value
    PVal_unique                       <- PairedTtest_unique$p.value
    UnPairedPval_unique               <- UnPairedTtest_unique$p.value
    LossAllele                 <- c(HLA_A_type1,HLA_A_type2)[ifelse(HLAtype1Log2MedianCoverageAtSites>HLAtype2Log2MedianCoverageAtSites,2,1)]
    KeptAllele                 <- c(HLA_A_type1,HLA_A_type2)[ifelse(HLAtype1Log2MedianCoverageAtSites>HLAtype2Log2MedianCoverageAtSites,1,2)]
    numMisMatchSitesCov        <- length(which(rowSums(misMatchCoveredInBoth)==2))
    tmp                        <- length(which(sort(log(c(HLA_A_type1tumorCov/HLA_A_type1normalCov),2)[as.character(missMatchseq1)]-log(c(HLA_A_type2tumorCov/HLA_A_type2normalCov),2)[as.character(missMatchseq2)])>0))/length(missMatchseq2)
    propSupportiveSites        <- max(tmp,1-tmp)*100
    
    # let's have a look at the additional features we can add in relation to integer copy numbers
    
    HLA_type1copyNum_withoutBAF          <- NA
    HLA_type1copyNum_withoutBAF_lower    <- NA
    HLA_type1copyNum_withoutBAF_upper    <- NA
    
    HLA_type1copyNum_withBAF             <- NA
    HLA_type1copyNum_withBAF_lower       <- NA
    HLA_type1copyNum_withBAF_upper       <- NA
    
    HLA_type2copyNum_withoutBAF          <- NA
    HLA_type2copyNum_withoutBAF_lower    <- NA
    HLA_type2copyNum_withoutBAF_upper    <- NA
    
    HLA_type2copyNum_withBAF             <- NA
    HLA_type2copyNum_withBAF_lower       <- NA
    HLA_type2copyNum_withBAF_upper       <- NA

    HLA_type1copyNum_withoutBAFBin       <- NA
    HLA_type1copyNum_withoutBAFBin_lower       <- NA
    HLA_type1copyNum_withoutBAFBin_upper       <- NA
    
    HLA_type1copyNum_withBAFBin          <- NA
    HLA_type1copyNum_withBAFBin_lower     <- NA
    HLA_type1copyNum_withBAFBin_upper     <- NA
    
    HLA_type2copyNum_withoutBAFBin       <- NA
    HLA_type2copyNum_withoutBAFBin_lower       <- NA
    HLA_type2copyNum_withoutBAFBin_upper       <- NA
    
    HLA_type2copyNum_withBAFBin          <- NA
    HLA_type2copyNum_withBAFBin_lower          <- NA
    HLA_type2copyNum_withBAFBin_upper          <- NA
    
    if(performIntegerCopyNum)
    {
      HLA_type1copyNum_withoutBAF        <- nA_rawVal_withoutBAF
      HLA_type1copyNum_withoutBAF_lower  <- nA_rawVal_withoutBAF_lower
      HLA_type1copyNum_withoutBAF_upper  <- nA_rawVal_withoutBAF_upper
      
      HLA_type1copyNum_withBAF           <- nA_rawVal_withBAF
      HLA_type1copyNum_withBAF_lower           <- nA_rawVal_withBAF_lower
      HLA_type1copyNum_withBAF_upper           <- nA_rawVal_withBAF_upper
      
      HLA_type2copyNum_withoutBAF        <- nB_rawVal_withoutBAF
      HLA_type2copyNum_withoutBAF_lower        <- nB_rawVal_withoutBAF_lower
      HLA_type2copyNum_withoutBAF_upper        <- nB_rawVal_withoutBAF_upper
      
      HLA_type2copyNum_withBAF           <- nB_rawVal_withBAF
      HLA_type2copyNum_withBAF_lower           <- nB_rawVal_withBAF_lower
      HLA_type2copyNum_withBAF_upper           <- nB_rawVal_withBAF_upper

      HLA_type1copyNum_withoutBAFBin     <- nA_rawVal_withoutBAFBin
      HLA_type1copyNum_withoutBAFBin_lower     <- nA_rawVal_withoutBAFBin_lower
      HLA_type1copyNum_withoutBAFBin_upper     <- nA_rawVal_withoutBAFBin_upper
      
      HLA_type1copyNum_withBAFBin        <- nA_rawVal_withBAF_bin
      HLA_type1copyNum_withBAFBin_lower        <- nA_rawVal_withBAF_bin_lower
      HLA_type1copyNum_withBAFBin_upper        <- nA_rawVal_withBAF_bin_upper
      
      HLA_type2copyNum_withoutBAFBin     <- nB_rawVal_withoutBAFBin
      HLA_type2copyNum_withoutBAFBin_lower     <- nB_rawVal_withoutBAFBin_lower
      HLA_type2copyNum_withoutBAFBin_upper     <- nB_rawVal_withoutBAFBin_upper
      
      HLA_type2copyNum_withBAFBin        <- nB_rawVal_withBAF_bin
      HLA_type2copyNum_withBAFBin_lower        <- nB_rawVal_withBAF_bin_lower
      HLA_type2copyNum_withBAFBin_upper        <- nB_rawVal_withBAF_bin_upper
    }

    
    
    
    out <- cbind(HLA_A_type1,HLA_A_type2,HLAtype1Log2MedianCoverage,HLAtype2Log2MedianCoverage,HLAtype1Log2MedianCoverageAtSites,HLAtype2Log2MedianCoverageAtSites
                 ,HLA_type1copyNum_withoutBAF
                 ,HLA_type1copyNum_withoutBAF_lower
                 ,HLA_type1copyNum_withoutBAF_upper
                 ,HLA_type1copyNum_withBAF
                 ,HLA_type1copyNum_withBAF_lower
                 ,HLA_type1copyNum_withBAF_upper
                 ,HLA_type2copyNum_withoutBAF
                 ,HLA_type2copyNum_withoutBAF_lower
                 ,HLA_type2copyNum_withoutBAF_upper
                 ,HLA_type2copyNum_withBAF
                 ,HLA_type2copyNum_withBAF_lower
                 ,HLA_type2copyNum_withBAF_upper
                 ,HLA_type1copyNum_withoutBAFBin
                 ,HLA_type1copyNum_withoutBAFBin_lower
                 ,HLA_type1copyNum_withoutBAFBin_upper
                 ,HLA_type1copyNum_withBAFBin
                 ,HLA_type1copyNum_withBAFBin_lower
                 ,HLA_type1copyNum_withBAFBin_upper
                 ,HLA_type2copyNum_withoutBAFBin
                 ,HLA_type2copyNum_withoutBAFBin_lower
                 ,HLA_type2copyNum_withoutBAFBin_upper
                 ,HLA_type2copyNum_withBAFBin
                 ,HLA_type2copyNum_withBAFBin_lower
                 ,HLA_type2copyNum_withBAFBin_upper
                 ,PVal,UnPairedPval,PVal_unique,UnPairedPval_unique,LossAllele,KeptAllele,numMisMatchSitesCov,propSupportiveSites) 
    HLAoutPut <- rbind(HLAoutPut,out)
    
    
  }
  dev.off()
  
  regionSpecOutPut <- cbind(region,HLAoutPut)
  PatientOutPut    <- rbind(PatientOutPut,regionSpecOutPut)
  
}  

HLAoutLoc <- paste(workDir,full.patient,'.',minCoverageFilter,".DNA.HLAlossPrediction_CI.xls",sep="")
write.table(PatientOutPut,file=HLAoutLoc,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


if(performIntegerCopyNum)
{
  HLABAFsummaryLoc <- paste(workDir,full.patient,'.', minCoverageFilter,".DNA.IntegerCPN_CI.xls",sep="")
  write.table(combinedTable,file=HLABAFsummaryLoc,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
}




if(cleanUp){

  cmd <- paste('rm ', workDir, '*tumor*', sep = '')
  # system(cmd)

  cmd <- paste('rm ', workDir, '*normal*', sep = '')
  # system(cmd)

  cmd <- paste('rm ', workDir, '*/*sam', sep = '')
  system(cmd)

  cmd <- paste('rm ', workDir, '*/*fastq', sep = '')
  system(cmd)

  cmd <- paste('rm ', workDir, '*/*reads', sep = '')
  system(cmd)

  cmd <- paste('rm ', workDir, '*/*temp*bam', sep = '')
  system(cmd)

  cmd <- paste('rm ', workDir, '*/*chr6region.patient.reference.hlas.csorted.bam', sep = '')
  system(cmd)

  cmd <- paste('rm ', workDir, '*/*chr6region.patient.reference.hlas.csorted.noduplicates.bam', sep = '')
  system(cmd)

  cmd <- paste('rm ', workDir, '*/*chr6region.patient.reference.hlas.bam', sep = '')
  system(cmd)

  cmd <- paste('rm ', workDir, '*/*type*[0-9].bam', sep = '')
  # system(cmd)

  cmd <- paste('rm ', workDir, '*/*type*[0-9].bam.bai', sep = '')
  # system(cmd)

}





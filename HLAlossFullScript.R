# before running
system("ml BEDTools/2.26.0-foss-2016b", wait=TRUE)
system("ml SAMtools/1.3.1-foss-2016b", wait=TRUE)
system("ml R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3", wait=TRUE)
system("ml novoalign/3.07.00", wait=TRUE)
system("ml TracerX-Picard-GATK/0.1-Java-1.7.0_80", wait=TRUE)


##########################
# Command line arguments #
##########################

 cmdArgs         <- commandArgs(trailingOnly = TRUE);

# full.patient    <- cmdArgs[1]
# workDir         <- cmdArgs[2]
# normalBAMfile   <- cmdArgs[3]
# BAMDir          <- cmdArgs[4]
# hlaPath         <- cmdArgs[5]
# HLAfastaLoc     <- cmdArgs[6]
# CopyNumLoc      <- cmdArgs[7]
# mapping.step    <- as.logical(cmdArgs[8])
# cleanUp         <- as.logical(cmdArgs[9])
# gc.correction.step <- as.logical(cmdArgs[10])
#TBK

full.patient          <- cmdArgs[1]
workDir               <- paste("/camp/lab/swantonc/working/watkint/scripts/hla_paper/lohhla_runs_2/",full.patient,"/exome/NeoAntigen/LOH/", sep="")
hlaPath               <- paste('/camp/lab/swantonc/working/rosentr/projects/neoantigen/tx100/samples-20160818/',substring(full.patient, 3),'/',substring(full.patient, 3),'.polysolver/winners.hla.txt', sep="")
normalBAMfile         <- paste('/farm/tracerx/lung/release/',full.patient,'/exome/BAM/processed/',full.patient,'_BS_GL.bam', sep="")
BAMDir                <- paste('/farm/tracerx/lung/release/',full.patient,'/exome/BAM/processed/', sep="")
HLAfastaLoc           <- paste("/farm/home/lr-tct-lif/wilson52/installs/polysolver/data/abc_complete.fasta", sep="")
mapping.step          <- FALSE
gc.correction.step <- TRUE
gc.correct.overwrite <- TRUE
cleanUp               <- FALSE
CopyNumLoc            <- paste('/farm/tracerx/lung/release/',full.patient,'/exome/ASCAT/solutions.txt', sep="")
#New gc information:
pt.mpileup.logr.dir <-  paste("/camp/lab/swantonc/working/watkint/scripts/hla_paper/txscratch_2/release/", full.patient, "/exome/mpileup/", sep="")#Perhaps mpileupDir
filt.tx.1000G.positions.path <-  "/camp/lab/swantonc/working/watkint/scripts/hla_paper/txscratch/release/B_LTX038/exome/mpileup/filt_tx_1000G_positions.txt"# bedPosFile
genome.gc.path <- "/camp/lab/swantonc/working/watkint/scripts/gc_correction/SureSelectV5_TRACERx_Edition_padded_1000_genome_no_prob_loci_GC_vals.txt"#GC file path
bed.file.path <- "/camp/lab/swantonc/working/watkint/scripts/gc_correction/SureSelectV5_TRACERx_Edition_padded_1000_genome_no_prob_loci.bed"
run.samtools.depth.flag <- FALSE
hla.names.sizes.path <- "/camp/lab/swantonc/working/watkint/scripts/gc_correction/abc_sizes_genome.txt"
reference <- "/camp/stp/babs/working/data/genomes/homo_sapiens/gatk/hg19/GATK_bundle_2.8/hg19/ucsc.hg19.fasta"
hla.names.sizes <- read.table(hla.names.sizes.path, stringsAsFactors=FALSE, sep="\t", header=FALSE)

hla.reference <- "/camp/lab/swantonc/working/watkint/scripts/gc_correction//abc_complete.fasta"




#############
# libraries #
#############

require(seqinr, quietly = TRUE)
require(Biostrings, quietly = TRUE)
require(beeswarm, quietly = TRUE)
require(zoo, quietly = TRUE)
require(Rsamtools, quietly = TRUE)
library(Rsamtools)#Added for GC correction.
library(GenomicRanges)#Added for GC correction



###########
# inputs #
###########

interactive           <- FALSE
if(interactive)
{
  full.patient          <- "B_LTX038"
  workDir               <- "/camp/lab/swantonc/working/watkint/scripts/hla_paper/lohhla_runs/B_LTX038/"
  hlaPath               <- '/camp/lab/swantonc/working/rosentr/projects/neoantigen/tx100/samples-20160818/LTX038/LTX038.polysolver/winners.hla.txt'
  normalBAMfile         <- '/farm/tracerx/lung/release/B_LTX038/exome/BAM/processed/B_LTX038_BS_GL.bam'
  BAMDir                <- '/farm/tracerx/lung/release/B_LTX038/exome/BAM/processed/'
  HLAfastaLoc           <- "/farm/home/lr-tct-lif/wilson52/installs/polysolver/data/abc_complete.fasta"
  mapping.step          <- FALSE
  gc.correction.step <- TRUE
  gc.correct.overwrite <- TRUE
  cleanUp               <- FALSE
  CopyNumLoc            <- '/farm/tracerx/lung/release/B_LTX038/exome/ASCAT/solutions.txt'


  # full.patient          <- "A_LTX049"
  # workDir               <- "/camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/tx100-noPoly/A_LTX049/rna/NeoAntigen/LOH/"
  # hlaPath               <- '/camp/lab/swantonc/working/rosentr/projects/neoantigen/tx100/samples-20160818/LTX049/LTX049.polysolver/winners.hla.txt'
  # normalBAMfile         <- 'FALSE'
  # BAMDir                <- '/camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/tx100-noPoly/A_LTX049/rna/BAM/'
  # HLAfastaLoc           <- "/camp/lab/swantonc/working/rosentr/data/IMGT/hla_abc_complete.rna.fasta"
  # mapping.step          <- TRUE
  # cleanUp               <- FALSE
  # CopyNumLoc            <- 'FALSE'
}

print(full.patient)
system('echo ${SLURM_JOBID}')

figureDir <- paste(workDir,"/Figures/",sep="")

log.name  <- paste(workDir, "running.hla.loh.exome@", gsub(":", "-", gsub(" +", "_", date())), "_log.txt", sep="")


runWithNormal           <- TRUE
if(normalBAMfile == 'FALSE'){
  runWithNormal         <- FALSE
}

minCoverageFilter       <- 30
extractNONmismatchReads <- TRUE
extractUniqueReads      <- TRUE

performIntegerCopyNum   <- TRUE
useLogRbin              <- TRUE
if(CopyNumLoc == 'FALSE'){
  performIntegerCopyNum <- FALSE
  useLogRbin            <- FALSE

}

gamma                   <- 1
binSize                 <- 150

GATKDir    <- '/camp/apps/eb/software/TracerX-Picard-GATK/0.1-Java-1.7.0_80/bin/'
NOVODir    <- "/camp/apps/eb/software/novoalign/3.07.00/bin/"
HLAexonLoc <- '/camp/lab/swantonc/working/rosentr/data/IMGT/hla.dat'

#############
# functions #
#############


ascat.loadData <- function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22,"X","Y"), gender = NULL, sexchromosomes = c("X","Y"), read.from.disk=FALSE) {

if(read.from.disk){
  # read in SNP array data files
  print.noquote("Reading Tumor LogR data...")
  Tumor_LogR <- read.table(Tumor_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
  print.noquote("Reading Tumor BAF data...")
  Tumor_BAF <- read.table(Tumor_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
 }
 if(!read.from.disk){
  Tumor_LogR <- Tumor_LogR_file
  Tumor_BAF <- Tumor_BAF_file
 }
  #infinite values are a problem - change those
  Tumor_LogR[Tumor_LogR==-Inf]=NA
  Tumor_LogR[Tumor_LogR==Inf]=NA
  Germline_LogR = NULL
  Germline_BAF = NULL

  if(!is.null(Germline_LogR_file)) {
    if(read.from.disk){
    # read in SNP array data files
     print.noquote("Reading Germline LogR data...")
     Germline_LogR <- read.table(Germline_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
     print.noquote("Reading Germline BAF data...")
     Germline_BAF <- read.table(Germline_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
  }
  if(!read.from.disk){
    Germline_LogR <- Germline_LogR_file
    Germline_BAF <- Germline_BAF_file
  }


    #infinite values are a problem - change those
    Germline_LogR[Germline_LogR==-Inf]=NA
    Germline_LogR[Germline_LogR==Inf]=NA
  }

  # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22,X,Y (or whatever is given in the input value of chrs)
  print.noquote("Registering SNP locations...")
  SNPpos <- Tumor_LogR[,1:2]
  SNPpos = SNPpos[SNPpos[,1]%in%chrs,]

  # if some chromosomes have no data, just remove them
  chrs = intersect(chrs,unique(SNPpos[,1]))

  Tumor_LogR = Tumor_LogR[rownames(SNPpos),c(-1,-2),drop=F]
  Tumor_BAF = Tumor_BAF[rownames(SNPpos),c(-1,-2),drop=F]
  # make sure it is all converted to numerical values
  for (cc in 1:dim(Tumor_LogR)[2]) {
    Tumor_LogR[,cc]=as.numeric(as.vector(Tumor_LogR[,cc]))
    Tumor_BAF[,cc]=as.numeric(as.vector(Tumor_BAF[,cc]))
  }

  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[rownames(SNPpos),c(-1,-2),drop=F]
    Germline_BAF = Germline_BAF[rownames(SNPpos),c(-1,-2),drop=F]
    for (cc in 1:dim(Germline_LogR)[2]) {
      Germline_LogR[,cc]=as.numeric(as.vector(Germline_LogR[,cc]))
      Germline_BAF[,cc]=as.numeric(as.vector(Germline_BAF[,cc]))
    }
  }
 
  # sort all data by genomic position
  last = 0;
  ch = list();
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke = SNPpos[SNPpos[,1]==chrs[i],]
    chrpos = chrke[,2]
    names(chrpos) = rownames(chrke)
    chrpos = sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))  
    SNPorder[ch[[i]]] = names(chrpos)
    last = last+length(chrpos)
  }
  SNPpos = SNPpos[SNPorder,]
  Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
  Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]

  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[SNPorder,,drop=F]
    Germline_BAF = Germline_BAF[SNPorder,,drop=F]
  }

  # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
  print.noquote("Splitting genome in distinct chunks...")
  chr = split_genome(SNPpos)

  if (is.null(gender)) {
    gender = rep("XX",dim(Tumor_LogR)[2])
  }
  return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
              Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, 
              Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
              SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs, 
              samples = colnames(Tumor_LogR), gender = gender, 
              sexchromosomes = sexchromosomes,
              failedarrays = NULL))
}


#Calculate logrs:
# note that probes not present in the GCcontentfile will be lost from the results
ascat.GCcorrect <- function(ASCATobj, GCcontentfile = NULL, read.from.disk =FALSE) {
  if(is.null(GCcontentfile)) {
    print.noquote("Error: no GC content file given!")
  }
  else {
    if(read.from.disk){GC_newlist<-read.table(file=GCcontentfile,header=TRUE,as.is=TRUE)}
    if(!read.from.disk){GC_newlist <- GCcontentfile}
    colnames(GC_newlist)[c(1,2)] = c("Chr","Position")
    GC_newlist$Chr<-as.character(GC_newlist$Chr)
    GC_newlist$Position<-as.numeric(as.character(GC_newlist$Position))

    ovl = intersect(row.names(ASCATobj$Tumor_LogR),row.names(GC_newlist))
    
    GC_newlist<-GC_newlist[ovl,]

    SNPpos = ASCATobj$SNPpos[ovl,]
    Tumor_LogR = ASCATobj$Tumor_LogR[ovl,,drop=F]
    Tumor_BAF = ASCATobj$Tumor_BAF[ovl,,drop=F]
    
    Germline_LogR = NULL
    Germline_BAF = NULL
    if(!is.null(ASCATobj$Germline_LogR)) {
      Germline_LogR = ASCATobj$Germline_LogR[ovl,,drop=F]
      Germline_BAF = ASCATobj$Germline_BAF[ovl,,drop=F]
    }

    for (s in 1:length(ASCATobj$samples)) {
      print.noquote(paste("Sample ", ASCATobj$samples[s], " (",s,"/",length(ASCATobj$samples),")",sep=""))
      Tumordata = Tumor_LogR[,s]
      names(Tumordata) = rownames(Tumor_LogR)

      # Calculate weighted correlation
      length_tot<-NULL
      corr_tot<-NULL
      for(chrindex in unique(SNPpos[,1])) {
        GC_newlist_chr<-GC_newlist[GC_newlist$Chr==chrindex,]
        td_chr<-Tumordata[GC_newlist$Chr==chrindex]

        flag_nona<-(complete.cases(td_chr) & complete.cases(GC_newlist_chr)) 
        corr<-cor(GC_newlist_chr[flag_nona,3:ncol(GC_newlist_chr)],td_chr[flag_nona])
        corr_tot<-cbind(corr_tot,corr)
        length_tot<-c(length_tot,length(td_chr))
      }
      corr<-apply(corr_tot,1,function(x) sum(abs(x*length_tot))/sum(length_tot))
      index_1M<-c(which(names(corr)=="X1M"),which(names(corr)=="X1Mb"))
      maxGCcol_short<-which(corr[1:(index_1M-1)]==max(corr[1:(index_1M-1)]))
      maxGCcol_long<-which(corr[index_1M:length(corr)]==max(corr[index_1M:length(corr)]))
      maxGCcol_long<-(maxGCcol_long+(index_1M-1))    
   
      cat("weighted correlation: ",paste(names(corr),format(corr,digits=2), ";"),"\n")   
      cat("Short window size: ",names(GC_newlist)[maxGCcol_short+2],"\n")
      cat("Long window size: ",names(GC_newlist)[maxGCcol_long+2],"\n")

      # Multiple regression 
      flag_NA<-(is.na(Tumordata))|(is.na(GC_newlist[,2+maxGCcol_short]))|(is.na(GC_newlist[,2+maxGCcol_long]))
      td_select<-Tumordata[!flag_NA]
      GC_newlist_select <- GC_newlist[!flag_NA,]
      x1<-GC_newlist_select[,2+maxGCcol_short]
      x2<-GC_newlist_select[,2+maxGCcol_long]
      x3<-(x1)^2
      x4<-(x2)^2
      model<-lm(td_select~x1+x2+x3+x4,y=TRUE)
      GCcorrected<-Tumordata
      GCcorrected[]<-NA     
      GCcorrected[!flag_NA] <- model$residuals
    
      Tumor_LogR[,s] = GCcorrected
    }

  # add some plotting code for each sample while it is generated!!!!
  
    return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
                Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, 
                Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
                SNPpos = SNPpos, ch = ASCATobj$ch, chr = ASCATobj$chr, chrs = ASCATobj$chrs, 
                samples = colnames(Tumor_LogR), gender = ASCATobj$gender, 
                sexchromosomes = ASCATobj$sexchromosomes))
  }  
}


# helper function to split the genome into parts
split_genome <- function(SNPpos) {

  # look for gaps of more than 5Mb (arbitrary treshold to account for big centremeres or other gaps) and chromosome borders
  bigHoles = which(diff(SNPpos[,2])>=5000000)+1
  chrBorders = which(SNPpos[1:(dim(SNPpos)[1]-1),1]!=SNPpos[2:(dim(SNPpos)[1]),1])+1

  holes = unique(sort(c(bigHoles,chrBorders)))

  # find which segments are too small
  #joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)

  # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
  #while (1 %in% joincandidates) {
  #  holes=holes[-1]
  #  joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  #}
  #while ((length(holes)+1) %in% joincandidates) {
  #  holes=holes[-length(holes)]
  #  joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  #}
 
  #while(length(joincandidates)!=0) {
    # the while loop is because after joining, segments may still be too small..

    #startseg = c(1,holes)
    #endseg = c(holes-1,dim(SNPpos)[1])

    # for each segment that is too short, see if it has the same chromosome as the segments before and after
    # the next always works because neither the first or the last segment is in joincandidates now
    #previoussamechr = SNPpos[endseg[joincandidates-1],1]==SNPpos[startseg[joincandidates],1] 
    #nextsamechr = SNPpos[endseg[joincandidates],1]==SNPpos[startseg[joincandidates+1],1]

    #distanceprevious = SNPpos[startseg[joincandidates],2]-SNPpos[endseg[joincandidates-1],2]
    #distancenext = SNPpos[startseg[joincandidates+1],2]-SNPpos[endseg[joincandidates],2]

    # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
    #joins = ifelse(previoussamechr&nextsamechr, 
    #               ifelse(distanceprevious>distancenext, joincandidates, joincandidates-1),
    #               ifelse(nextsamechr, joincandidates, joincandidates-1))

    #holes=holes[-joins]

    #joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  #}
  # if two neighboring segments are selected, this may make bigger segments then absolutely necessary, but I'm sure this is no problem.

  startseg = c(1,holes)
  endseg = c(holes-1,dim(SNPpos)[1])

  chr=list()
  for (i in 1:length(startseg)) {
    chr[[i]]=startseg[i]:endseg[i]
  }
  
  return(chr)
}


##
##Params:
##
##Returns:
##
getHLAPositionsInRef <- function(HLA_gene){

  if (HLA_gene=="hla_a") {
   # hla.hg19.ref.pos <- "chr6:29909037-29913661"
    hla.hg19.ref.pos <- c("chr6", 29909037,29913661)
  } else if (HLA_gene=="hla_b") {
    #hla.hg19.ref.pos <- "chr6:31321649-31324964"
    hla.hg19.ref.pos <- c("chr6", 31321649,31324964)
  } else if (HLA_gene=="hla_c") {
    #hla.hg19.ref.pos <- "chr6:31236526-31239869"
    hla.hg19.ref.pos <- c("chr6", 31236526,31239869)
  }
  return(hla.hg19.ref.pos)
}


##
##Params:
##
##Returns:
##
getHLAWindowGC <- function(hla.pos, HLA_gene,HLA_name, HLA.size, reference, hla.reference, base.change){
  win.start.pos <- hla.pos - base.change
  win.end.pos <- hla.pos + base.change#This can be larger than the size of the hla

  hla.hg19.ref.pos <- getHLAPositionsInRef(HLA_gene)

  if ( win.start.pos < 1 ) {
    pre.hla.pos <- as.numeric(hla.hg19.ref.pos[2]) + win.start.pos
    hla.pre.snp.pos <- 1
  } else {
    pre.hla.pos <- NA
    hla.pre.snp.pos <- win.start.pos
  }
  if ( win.end.pos > HLA.size ) {
    post.hla.pos <- as.numeric(hla.hg19.ref.pos[3]) + (win.end.pos - HLA.size)
    hla.post.snp.pos <- HLA.size
  } else {
    post.hla.pos  <-  NA
    hla.post.snp.pos <- win.end.pos
  }
  # hla.names.sizes <- read.table("/camp/lab/swantonc/working/watkint/scripts/gc_correction/abc_sizes_genome.txt", stringsAsFactors=FALSE, sep="\t", header=FALSE)
  # hla.seqlengths <- setNames(hla.names.sizes$V2, hla.names.sizes$V1)

  #hg19 seqlengths, might need changing for hg38 or be an option that can be passed.
  seqlengths <- setNames(c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566),
                         c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'))

  final.window.seq <- ""
  hla.gr <- GRanges(seqnames = Rle(HLA_name), ranges = IRanges(start=hla.pre.snp.pos, end = hla.post.snp.pos), seqlengths=hla.seqlengths)
  hla.reg.seq <- scanFa(hla.reference, hla.gr)
  final.window.seq <- hla.reg.seq
  #Of any relevant pre-hla region:
  if( ! is.na(pre.hla.pos) ){
    pre.hla.gr <- GRanges(seqnames = Rle("chr6"),
      ranges = IRanges(start=as.numeric(pre.hla.pos) -1 , end = as.numeric(hla.hg19.ref.pos[2]) - 1), 
      seqlengths=seqlengths)
    pre.hla.seq <- scanFa(reference, pre.hla.gr)
    #Update the seq if possible:
    final.window.seq <- paste(as.character(pre.hla.seq), as.character(hla.reg.seq), sep="")
  }
  #Of any relevant post-hla region
  if( ! is.na(post.hla.pos) ){
    post.hla.gr <- GRanges(seqnames = Rle("chr6"),
      ranges = IRanges(start=as.numeric(hla.hg19.ref.pos[3]) + 1, end = post.hla.pos + 1), 
      seqlengths=seqlengths)
    post.hla.seq <- scanFa(reference, post.hla.gr)
    final.window.seq <- paste(final.window.seq, as.character(post.hla.seq), sep="")
  }
  final.window.seq.dss <- DNAStringSet(final.window.seq, start=NA, end=nchar(final.window.seq), width=nchar(final.window.seq))
  alph <- alphabetFrequency(final.window.seq.dss, as.prob=TRUE)
  return(sum(alph[,c("G", "C")]))

}


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
  print("workDir")
  print(workDir)

  if(!override){
    outDir     <- paste(workDir, '/flagstat/', sep = '')
    outDir <- paste("/farm/tracerx/lung/release/", full.patient, "/exome/QC/flagstat/", sep="")
    if( !file.exists(outDir)){
      if( !dir.create(outDir, recursive = TRUE) ){
        stop(paste("Unable to create directory: ",outDir, "!\n", sep = ''))
      }
    }

    BAMs       <- list.files(BAMDir, pattern = 'bam$', full.names = TRUE)

    for(BAM in BAMs){
      region <- unlist(strsplit(BAM, split = '/'))[length(unlist(strsplit(BAM, split = '/')))]
      cmd    <- paste('samtools flagstat ', BAM, ' > ', outDir, region, '.proc.flagstat', sep = '')
      print(cmd)
      #system(cmd)
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

    # samtoolsCMD <- paste("samtools view -f 4 ",BAMDir, '/', BAMfile, ' >> ', regionDir, '/', BAMid, ".chr6region.sam ", sep = "")
    # write.table(paste(samtoolsCMD, '\n', sep = ''), file = log.name, row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    # system(samtoolsCMD)
    
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
      passed.reads <- count.events(paste(regionDir, '/', BAMid, '.type.', allele, '.bam', sep = ''), n = 1)
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
    #system(cmd)
    
  }

  
}
# also extract number of unique reads sequenced in tumor and normal
if(runWithNormal){
  regionUniqMappedRegions <- getUniqMapReads(workDir = workDir, BAMDir = BAMDir, override = FALSE)
  GermLineUniqMappedReads <- regionUniqMappedRegions[[grep("GL",names(regionUniqMappedRegions),value=TRUE)]]
}


####################################
# compare coverage between alleles # 
####################################

###Loading data in for the GC correction step##################
########added by TBKW 20170711################################
if (gc.correction.step) {

  #hla.names.sizes <- read.table("/camp/lab/swantonc/working/watkint/scripts/gc_correction/abc_sizes_genome.txt", stringsAsFactors=FALSE, sep="\t", header=FALSE)#At top of script now (also in function, needs to be rewritten)
  hla.seqlengths <- setNames(hla.names.sizes$V2, hla.names.sizes$V1)
  #reference <- "/camp/stp/babs/working/data/genomes/homo_sapiens/gatk/hg19/GATK_bundle_2.8/hg19/ucsc.hg19.fasta"#at top of script now
  #hla.reference <- "/camp/lab/swantonc/working/watkint/scripts/gc_correction//abc_complete.fasta"#At top of script now
  #Load in regular GC file:
  #genome.gc.path <- "/camp/lab/swantonc/working/watkint/scripts/gc_correction/SureSelectV5_TRACERx_Edition_padded_1000_genome_no_prob_loci_GC_vals.txt"#At top of script now
  orig.genome.gc.df <- read.table(genome.gc.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
  #filt.tx.1000G.positions.path <- "/camp/lab/swantonc/working/watkint/scripts/hla_paper/txscratch/release/B_LTX038/exome/mpileup/filt_tx_1000G_positions.txt"#At top of script now
  pt.logr.file.path <- paste(pt.mpileup.logr.dir, full.patient, "_tumor_LogR.Rdata", sep="")
  #pt.mpileup.logr.dir <- paste("/camp/lab/swantonc/working/watkint/scripts/hla_paper/txscratch/release/", full.patient, "/exome/mpileup/", sep="")#At the top of script now

  if ( (! file.exists(pt.logr.file.path)) |  gc.correct.overwrite){

    print("Patient genome wide LogR file doesn't exist, generating...")
    #This assumes that the samtools depth jobs for all regions being examined have been run.
    #As of 20170711 all regions which made if through to PyClone for the TX100 paper have an "_a_depth.txt" file produced by querying depth based on
    #1000 genomes positions that overlap with the TRACERx capture kit
    mpileup.files <- grep("_a_depth.txt$",list.files(pt.mpileup.logr.dir),value=TRUE)

    if(length(mpileup.files) == 0){
      print('Samtools depth output missing, cannot perform GC correction.')

      if (run.samtools.depth.flag) {

        print('Running samtools depth commands:')

        all.files.vec <- list.files(path=BAMDir, full.names=TRUE)
        bam.files.vec <- grep("bam$", all.files.vec, value=TRUE)
        if(!dir.exists(pt.mpileup.logr.dir)){
          dir.create(pt.mpileup.logr.dir, recursive=TRUE)
        }

        for ( j in 1:length(bam.files.vec) ) {

          depth.command <- paste("samtools depth -b ",bed.file.path," -a ",bam.files.vec[j] , " > ",pt.mpileup.logr.dir,gsub(".bam$", "_a_depth.txt", basename(bam.files.vec[j]))," ", sep="")
          print(depth.command)
          system(depth.command, wait=TRUE)

        }

      }


    } else {

      print("Rerun with samtools depth run enabled")

    }
    #Need to regrep incase these where gust generated
    mpileup.files <- grep("_a_depth.txt$",list.files(pt.mpileup.logr.dir),value=TRUE)

    #germline pileup path:
    germline.pileup <- grep("_GL", mpileup.files, value=TRUE)
    germline.pileup.path <- paste(pt.mpileup.logr.dir, germline.pileup, sep="")
    #tumour pileup paths:
    tumour.pileups <- grep(germline.pileup, mpileup.files, value=TRUE, invert=TRUE)
    tumour.pileup.paths <- paste(pt.mpileup.logr.dir, tumour.pileups, sep="")
    all.mpileup.paths <- c(germline.pileup,tumour.pileups)
    for ( i in 1:length(all.mpileup.paths) ){
      #This call makes a new column that will allow the subsetting by the composite chrom_pos column later on (concatenates chr and position into a new column)
      system(paste("awk -v OFS=\'\t\' \'{print $1\"_\"$2,$1, $2, $3}\' ", pt.mpileup.logr.dir,all.mpileup.paths[i], " > ", pt.mpileup.logr.dir,gsub("_a_depth.txt", "_a_chr_pos_depth.txt", all.mpileup.paths[i]), sep=""))
      #This line makes sure that every position in the depth file matches one from the filtered 1000G positions path, this wasn't done in R as it was too slow.
      system(paste("awk \'NR == FNR{a[$1]=$0; next};($1 in a){print $0; next}\' ",filt.tx.1000G.positions.path," ", pt.mpileup.logr.dir,gsub("_a_depth.txt", "_a_chr_pos_depth.txt",all.mpileup.paths[i]), " > ", pt.mpileup.logr.dir,gsub("_a_depth.txt", "_a_chr_pos_filt_depth.txt",all.mpileup.paths[i]), sep=""))
    }

    germline.pileup.path <- gsub("_a_depth.txt", "_a_chr_pos_filt_depth.txt", germline.pileup.path)
    tumour.pileup.paths <- gsub("_a_depth.txt", "_a_chr_pos_filt_depth.txt", tumour.pileup.paths)

    gl.mp.df <- read.table(germline.pileup.path, stringsAsFactors=FALSE, sep="\t", header=FALSE, fill=TRUE)
    gc.path <- genome.gc.path
    gc.df <- read.table(gc.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
    gc.df$chrom_pos <- paste(gc.df$Chr, "_", gc.df$Position, sep="")
    head(gc.df)
    common.pos <- gl.mp.df[,1]
    tum.pos.list <- list()
    for ( i in 1:length(tumour.pileup.paths) ) {
      temp.df <- read.table(tumour.pileup.paths[[i]], stringsAsFactors=FALSE, sep="\t", header=FALSE)
      common.pos <- common.pos[common.pos%in%temp.df[,1]]
      tum.pos.list[[i]] <- temp.df
      #gc()
    }
    gl.mp.df <- gl.mp.df[gl.mp.df[,1]%in%common.pos,]
    for ( i in 1:length(tum.pos.list) ) {
      temp.df <- tum.pos.list[[i]]
      temp.df <- temp.df[temp.df[,1]%in%common.pos,]
      LogR <- temp.df$V4 / (gl.mp.df$V4+0.0001)
      scaling.factor <- mean(LogR)
      LogR.scaled <- log2(LogR / scaling.factor)
      temp.df$LogR <- LogR.scaled
      tum.pos.list[[i]] <- temp.df
    }
    all.reg.logr <- tum.pos.list
    names(all.reg.logr) <- gsub("_a_depth.txt", "", tumour.pileup.paths)
    print("Saving.....")
    all.logr.file.path <- paste(pt.mpileup.logr.dir, full.patient, "_tumor_LogR.Rdata", sep="")
    save(all.reg.logr, file=all.logr.file.path)

  } else {

    print("Patient genome wide LogR file exists, loading...")
    all.logr.file.path <- paste(pt.mpileup.logr.dir, full.patient, "_tumor_LogR.Rdata", sep="")
    load(all.logr.file.path)
    #Now all.reg.logr should be in the environment
  }
  #load("/camp/lab/swantonc/working/watkint/scripts/gc_correction/B_LTX038_tumor_LogR.Rdata")
}###End of loading data in for the GC correction step##################



normalName <- regions[which(paste(BAMDir, regions, '.bam', sep = '') == normalBAMfile)]

# let's load the winners
# next, we can look at each mpileupFile, and assess whether we see differences in coverage between the two. 
# let's look at a region of interest. 
PatientOutPut <- c()
for (region in regions)
{
  print(region)
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
  
  HLAoutPut <- c()
  if(gc.correction.step){
    pdf(paste(figureDir,region,".minCoverage_",minCoverageFilter,".HLA_GC_corrected.pdf",sep=""),width=10,height=6)
  } else {
    pdf(paste(figureDir,region,".minCoverage_",minCoverageFilter,".HLA.pdf",sep=""),width=10,height=6)
  } 
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
    #TODO remove when sharing with NM and RR
    print(workDir)
    # workDir <- paste("/camp/lab/swantonc/working/watkint/scripts/hla_paper/lohhla_runs_2/",full.patient, "/exome/NeoAntigen/LOH/", sep="")
    # dir.create(workDir, recursive=TRUE)
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
      print(paste(workDir, "/",region,".",HLA_A_type1,".","tumor.mpileup",sep=""))
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
      print(paste(workDir,region,".",HLA_A_type1,".tumor.NoMissMatch.pileup",sep=""))
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
    
    if(extractUniqueReads == TRUE){
      
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
        abline(h=median(HLA_A_type1tumorCov),lty='dashed',col='#de2d26',lwd=1.5)
        abline(h=median(HLA_A_type2tumorCov),lty='dashed',col='#3182bd',lwd=1.5)
        barplot(c(rollmean(HLA_A_type2tumorCov,150)),ylim=c(0,Ymax),xaxt='n',main=HLA_A_type2,las=1,col='#3182bd',border='#3182bd50')
        abline(h=median(HLA_A_type1tumorCov),lty='dashed',col='#de2d26',lwd=1.5)
        abline(h=median(HLA_A_type2tumorCov),lty='dashed',col='#3182bd',lwd=1.5)
        

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
        #LogR calculated:
        for (i in 1:(length(seqToConsider)-1))
        {
          
          PotentialSites   <- as.character(seqToConsider[i]:seqToConsider[i+1])
          combinedBinTumor  <- median(as.numeric(as.numeric(HLA_A_type1tumorCov[names(HLA_A_type1tumorCov)%in%PotentialSites])))+median(as.numeric(as.numeric(HLA_A_type2tumorCov[names(HLA_A_type2tumorCov)%in%PotentialSites])))
          combinedBinNormal <- median(as.numeric(as.numeric(HLA_A_type1normalCov[names(HLA_A_type1normalCov)%in%PotentialSites])))+median(as.numeric(as.numeric(HLA_A_type2normalCov[names(HLA_A_type2normalCov)%in%PotentialSites])))
          combinedBinlogR   <- log2(combinedBinTumor/combinedBinNormal*MultFactor)
          type1BinlogR     <- median(log2(as.numeric(as.numeric(HLA_A_type1tumorCov[names(HLA_A_type1tumorCov)%in%PotentialSites])/as.numeric(HLA_A_type1normalCov[names(HLA_A_type1normalCov)%in%PotentialSites])*MultFactor)))
          type2BinlogR     <- median(log2(as.numeric(as.numeric(HLA_A_type2tumorCov[names(HLA_A_type2tumorCov)%in%PotentialSites])/as.numeric(HLA_A_type2normalCov[names(HLA_A_type2normalCov)%in%PotentialSites])*MultFactor)))
          binLogR <- rbind(binLogR,cbind(seqToConsider[i],seqToConsider[i+1],combinedBinlogR,type1BinlogR,type2BinlogR))
        }

        #Add in GC correction here:
        print("missMatchseq1")
        print(missMatchseq1)
        print("missMatchseq2")
        print(missMatchseq2)
        if (gc.correction.step & (length(missMatchseq1) != 0) ) {
          genome.gc.df <- orig.genome.gc.df
          genome.gc.df <- read.table(genome.gc.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)

          #So now here is where we should remove positions that aren't in the depth files form the utmour regions:
          #At some point this should actually go outside the for loop but for now keep it in.
          #Note all.reg.logr comes from loading in the appropriate patient B_LTX038_tumor_LogR.Rdata equivalent at the beginning of this loop
          region.depth.name <- grep(paste(region,"_", sep=""), names(all.reg.logr), value=TRUE)

          included.chr.pos <- all.reg.logr[[region.depth.name]]$V1
          #Remove any positions that we don't have precomputed GC content for:
          genome.gc.df$chr_pos <- paste(genome.gc.df$Chr, "_", genome.gc.df$Position, sep="")
          genome.gc.df <- genome.gc.df[genome.gc.df$chr_pos%in%included.chr.pos,]

          #Now we need to put in markers for the different hlas:
          #So now we want the closest snp and to add positions to it:
          hla.ref.pos <- getHLAPositionsInRef(HLA_gene)
          chr6.gc.df <- genome.gc.df[genome.gc.df$Chr=="chr6",]
          chr6.gc.df$dist <- abs(chr6.gc.df$Position - as.numeric(hla.ref.pos[2]))
          min.dist.idx <- which(chr6.gc.df$dist==min(chr6.gc.df$dist))
          min.dist.position <- chr6.gc.df[min.dist.idx,]$Position

          #Create the dummy positions that we'll give the SNPs in the HLA regions
          #Reason we do this is because HLAs canbe of varying sizes and won't necessarily fit the area designated in the reference genome
          #In addition the two different HLAs may also be of different sizes and thus can't share the same position,
          #Finally each hla must have the same position added to the dataframe to allow the two HLA specific logrs to undergo correction
          missMatchseq1.gc.df <- chr6.gc.df[rep(min.dist.idx, length(missMatchseq1)),]
          print(seq(1,length(missMatchseq1), 1))
          missMatchseq1.gc.df$Position  <- seq(1,length(missMatchseq1), 1) + as.numeric(chr6.gc.df[min.dist.idx,]$Position)
          missMatchseq2.gc.df <- chr6.gc.df[rep(min.dist.idx, length(missMatchseq2)),]
          missMatchseq2.gc.df$Position  <- seq(1,length(missMatchseq2), 1) + as.numeric(chr6.gc.df[min.dist.idx,]$Position) + length(missMatchseq2)

          #Parameters for GC content calculation:
          base.change <- c(50,100,200,400,800,1600,3200,6400,12800,25600,51200,102400)
          HLA_A_type1.size <- hla.names.sizes[hla.names.sizes$V1==HLA_A_type1,]$V2
          HLA_A_type2.size <- hla.names.sizes[hla.names.sizes$V1==HLA_A_type2,]$V2
          #Here we leave the largest windows unchanged as the variation within a single HLA region has an insignificant effect on the GC content of
          #These extremely large windows.
          #Calculate HLA GC content by window
          for ( j in 1:length(base.change) ) {
            # print(base.change[j])
            missMatchseq1.gc.df[, 2+j ] <- sapply(missMatchseq1, getHLAWindowGC, HLA_gene, HLA_A_type1, HLA_A_type1.size, reference, hla.reference, base.change[j])
            missMatchseq2.gc.df[, 2+j ] <- sapply(missMatchseq2, getHLAWindowGC, HLA_gene, HLA_A_type2,HLA_A_type2.size, reference, hla.reference, base.change[j])
          }
          #Values look like this 
          #  Chr Position    X100bp    X200bp    X400bp    X800bp   X1600bp   X3200bp
          # 1 chr1   762496 0.6930693 0.6815920 0.6508728 0.6329588 0.5883823 0.4995314
          #     X6400bp  X12800bp  X25600bp  X51200bp X102400bp X204800bp       X1M
          # 1 0.4452429 0.4662917 0.4625210 0.4453233 0.4458843 0.4732203 0.4906955
          #         X2M       X5M      X10M
          # 1 0.4724761 0.5113033 0.4910091

          all.missMatchseq1.positions <- c(missMatchseq1.gc.df$Position, missMatchseq2.gc.df$Position)
          new.chr6.gc.df <- chr6.gc.df[ ! chr6.gc.df$Position%in%all.missMatchseq1.positions,]

          #So now insert these into the dataframe:
          #We need the index
          start.idx <- rownames(genome.gc.df[genome.gc.df$Chr=="chr6" & genome.gc.df$Position==min.dist.position,])
          missMatchseq1.gc.df$dist <- NULL
          missMatchseq2.gc.df$dist <- NULL
          #Now we want to remove positions that are not in the tumor logR data but not remove the hla positions that we've just spent ages working out:
          #We also can't remove the chr6 start index, so therefore the 

          hla.gc.df <- rbind(genome.gc.df[1:as.numeric(start.idx),], missMatchseq1.gc.df, missMatchseq2.gc.df, genome.gc.df[(as.numeric(start.idx) + 1):nrow(genome.gc.df),])
          rownames(hla.gc.df) <- 1:nrow(hla.gc.df)

        } else if (!gc.correction.step) {#End of GC correction if statement
          print(HLA_gene)
          print("GC content calculation not performed as was not requested.")
        } else if(length(missMatchseq1) == 0) {
          print(HLA_gene)
          print("No mismatch positions, so no GC content within HLA calculated.")
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


        if (gc.correction.step & (length(missMatchseq1) !=0)) {
          #LogR data from the regions:
          #Modify this to get load in the appropriate region
          tmp <- all.reg.logr[[region.depth.name]]

          ####
          reg.logr <- tmp[,c(2,3,5)]
          reg.logr[,1] <- gsub("chr","", reg.logr[,1])
          colnames(reg.logr) <- c("Chr", "Position", "LogR")

          hla1.logr.df <- data.frame(Chr=rep(6, nrow(tmpOut)), Position=missMatchseq1.gc.df$Position, LogR=tmpOut[,"logR_type1"], stringsAsFactors=FALSE)
          hla2.logr.df <- data.frame(Chr=rep(6, nrow(tmpOut)), Position=missMatchseq2.gc.df$Position, LogR=tmpOut[,"logR_type2"], stringsAsFactors=FALSE)


          #Remove any snps from the precomputed positions that also match the HLA positions that we create as placeholders
          reg.logr$chrom_pos <- paste(reg.logr$Chr, "_", reg.logr$Position, sep="")
          hla1.logr.df$chrom_pos <- paste(hla1.logr.df$Chr, "_", hla1.logr.df$Position, sep="")
          hla2.logr.df$chrom_pos <- paste(hla2.logr.df$Chr, "_", hla2.logr.df$Position, sep="")
          reg.logr <- reg.logr[ ! reg.logr$chrom_pos %in% hla1.logr.df$chrom_pos,]
          reg.logr <- reg.logr[ ! reg.logr$chrom_pos %in% hla2.logr.df$chrom_pos,]
          reg.logr$chrom_pos <- NULL
          hla1.logr.df$chrom_pos <- NULL
          hla2.logr.df$chrom_pos <- NULL

          #Create the new composite logr dataframe that matches the gc content df to pass to gc correction
          reg.and.hla.logr.df <- rbind(reg.logr[1:as.numeric(start.idx),], hla1.logr.df, hla2.logr.df, reg.logr[(as.numeric(start.idx) + 1):nrow(reg.logr),])
          rownames(reg.and.hla.logr.df)  <- 1:nrow(reg.and.hla.logr.df)

          #Initialise the other data that ASCAT needs for GC correction, either empty (BAF) or with 0s for germline values.
          tum.baf <- data.frame(Chr=reg.and.hla.logr.df$Chr, Position=reg.and.hla.logr.df$Position, BAF=rep(0, nrow(reg.and.hla.logr.df)), stringsAsFactors=FALSE)
          gl.logr <- data.frame(Chr=reg.and.hla.logr.df$Chr, Position=reg.and.hla.logr.df$Position, LogR=rep(0, nrow(reg.and.hla.logr.df)), stringsAsFactors=FALSE)
          gl.baf <- data.frame(Chr=reg.and.hla.logr.df$Chr, Position=reg.and.hla.logr.df$Position, BAF=rep(0, nrow(reg.and.hla.logr.df)), stringsAsFactors=FALSE)

          #Set the gender as XY, ASCAT requires a gender but unimportant for calling HLA LOH
          gender <- "XY"
          #Load in data 
          ascat.bc = ascat.loadData(Tumor_LogR_file = reg.and.hla.logr.df,Tumor_BAF_file = tum.baf, Germline_LogR_file=gl.logr,Germline_BAF_file=gl.baf, gender=gender,chrs = c(1:22,"X","Y"),sexchromosomes = c("X","Y"))
          #Reformat for the GC correction step:        
          hla.gc.df$Chr <- gsub("chr", "", hla.gc.df$Chr)
          hla.gc.df$chr_pos <- NULL 
          #Gc correction
          ascat.bc.corr <- ascat.GCcorrect(ascat.bc, hla.gc.df)

          #So now we want to pull pull out the HLA SNPs that we corrected
          #How can we do this, should be able to work out the indices of the HLA snps in the 
          if(length(ascat.bc$Tumor_LogR$LogR) != length(ascat.bc.corr$Tumor_LogR$LogR)) {

            print("SNPs being lost in GC step, check for errors.")

          }

          #Make dataframes to pull out the HLA logr GC corrected values
          gc.uncorrected.df <- data.frame(Chr=ascat.bc$SNPpos$Chr, Position=ascat.bc$SNPpos$Position, GC_correct_LogR=ascat.bc$Tumor_LogR$LogR, stringsAsFactors=FALSE)
          gc.corrected.df <- data.frame(Chr=ascat.bc.corr$SNPpos$Chr, Position=ascat.bc.corr$SNPpos$Position, GC_correct_LogR=ascat.bc.corr$Tumor_LogR$LogR, stringsAsFactors=FALSE)
          #add chrom_pos columns to allow selection of the hla_snps
          gc.uncorrected.df$chrom_pos <- paste(gc.uncorrected.df$Chr, "_", gc.uncorrected.df$Position, sep="") 
          gc.corrected.df$chrom_pos <- paste(gc.corrected.df$Chr, "_", gc.corrected.df$Position, sep="")
          hla1.logr.df$chrom_pos <- paste(hla1.logr.df$Chr, "_", hla1.logr.df$Position, sep="")
          hla2.logr.df$chrom_pos <- paste(hla2.logr.df$Chr, "_", hla2.logr.df$Position, sep="")
          #add the logr to the hla dataframes
          hla1.logr.df$GC_correct_LogR <- gc.corrected.df[gc.corrected.df$chrom_pos%in%hla1.logr.df$chrom_pos,]$GC_correct_LogR
          hla2.logr.df$GC_correct_LogR <- gc.corrected.df[gc.corrected.df$chrom_pos%in%hla2.logr.df$chrom_pos,]$GC_correct_LogR

          #Now take the GC corrected LogR forward for the rest of the analysis.
          tmpOut.old <- tmpOut

          tmpOut[,"logR_type1"] <- hla1.logr.df$GC_correct_LogR
          tmpOut[,"logR_type2"] <- hla2.logr.df$GC_correct_LogR
          
          print("################################")
          print(HLA_gene)

          print("Mean HLA type1 unGCcorrected LogR:")
          print(mean(tmpOut.old[,"logR_type1"]))

          print("Mean HLA type1 GCcorrected LogR:")
          print(mean(tmpOut[,"logR_type1"]))

          print("Mean HLA type2 unGCcorrected LogR:")
          print(mean(tmpOut.old[,"logR_type2"]))

          print("Mean HLA type2 GCcorrected LogR:")
          print(mean(tmpOut[,"logR_type2"]))


          write.table(tmpOut.old,file=paste(workDir,full.patient, "_", HLA_A_type1, "_", region, "_unGCcorrected_LogR.txt", sep="") ,row.names=FALSE, col.names=TRUE,sep="\t", quote=FALSE)
          write.table(tmpOut,file=paste(workDir, full.patient, "_", HLA_A_type2, "_", region, "_GCcorrected_LogR.txt", sep="") ,row.names=FALSE, col.names=TRUE,sep="\t", quote=FALSE)

          print("################################")


         
        } else if (!gc.correction.step) {#End of GC correction if statement
          print(HLA_gene)
          print("GC content calculation not performed as was not requested.")
        } else if(length(missMatchseq1) == 0) {
          print(HLA_gene)
          print("No mismatch positions, so no GC content within HLA calculated.")
        }

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

        nB_rawVal_withBAF <- median(combinedTable$nBcombined)
        nA_rawVal_withBAF <- median(combinedTable$nAcombined)

        rawValsBin      <- funCalcN_withBAF(combinedTable$binlogRCombined,combinedTable$BAFcombined,tumorPloidy,tumorPurity,gamma)
        combinedTable$nAcombinedBin <- rawValsBin[,1]
        combinedTable$nBcombinedBin <- rawValsBin[,2]
        
        nB_rawVal_withBAF     <- median(combinedTable$nBcombined)
        nA_rawVal_withBAF     <- median(combinedTable$nAcombined)
        
        #let's only count non duplicates
        nB_rawVal_withBAF_bin <- median(combinedTable[!duplicated(combinedTable$binNum),]$nBcombinedBin)
        nA_rawVal_withBAF_bin <- median(combinedTable[!duplicated(combinedTable$binNum),]$nAcombinedBin)


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
        
        combinedTable$nAsep <- funCalcN_withoutBAF(combinedTable$logR_type1,tumorPloidy,tumorPurity,gamma)
        combinedTable$nAsepBin <- funCalcN_withoutBAF(combinedTable$binlogRtype1,tumorPloidy,tumorPurity,gamma)
        combinedTable$nBsep <- funCalcN_withoutBAF(combinedTable$logR_type2,tumorPloidy,tumorPurity,gamma)
        combinedTable$nBsepBin <- funCalcN_withoutBAF(combinedTable$binlogRtype2,tumorPloidy,tumorPurity,gamma)

        nB_rawVal_withoutBAF <- median(combinedTable$nBsep)
        nA_rawVal_withoutBAF <- median(combinedTable$nAsep)

        #median(combinedTable[!duplicated(combinedTable$binNum),]$nAcombinedBin)
        
        nB_rawVal_withoutBAFBin <-  median(combinedTable[!duplicated(combinedTable$binNum),]$nBsepBin)
        nA_rawVal_withoutBAFBin <-  median(combinedTable[!duplicated(combinedTable$binNum),]$nAsepBin)
        
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
    HLAtype1Log2MedianCoverage <- median(log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor))
    HLAtype2Log2MedianCoverage <- median(log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor))
    HLAtype1Log2MedianCoverageAtSites <- median(log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[names(HLA_A_type1tumorCov)%in%missMatchseq1])
    HLAtype2Log2MedianCoverageAtSites <- median(log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[names(HLA_A_type2tumorCov)%in%missMatchseq2])
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
    HLA_type1copyNum_withBAF             <- NA
    HLA_type2copyNum_withoutBAF          <- NA
    HLA_type2copyNum_withBAF             <- NA

    HLA_type1copyNum_withoutBAFBin       <- NA
    HLA_type1copyNum_withBAFBin          <- NA
    HLA_type2copyNum_withoutBAFBin       <- NA
    HLA_type2copyNum_withBAFBin          <- NA
    
    if(performIntegerCopyNum)
    {
      HLA_type1copyNum_withoutBAF        <- nA_rawVal_withoutBAF
      HLA_type1copyNum_withBAF           <- nA_rawVal_withBAF
      
      HLA_type2copyNum_withoutBAF        <- nB_rawVal_withoutBAF
      HLA_type2copyNum_withBAF           <- nB_rawVal_withBAF

      HLA_type1copyNum_withoutBAFBin     <- nA_rawVal_withoutBAFBin
      HLA_type1copyNum_withBAFBin        <- nA_rawVal_withBAF_bin
      
      HLA_type2copyNum_withoutBAFBin     <- nB_rawVal_withoutBAFBin
      HLA_type2copyNum_withBAFBin        <- nB_rawVal_withBAF_bin
    }

    
    
    
    out <- cbind(HLA_A_type1,HLA_A_type2,HLAtype1Log2MedianCoverage,HLAtype2Log2MedianCoverage,HLAtype1Log2MedianCoverageAtSites,HLAtype2Log2MedianCoverageAtSites
                 ,HLA_type1copyNum_withoutBAF,HLA_type1copyNum_withBAF
                 ,HLA_type2copyNum_withoutBAF,HLA_type2copyNum_withBAF
                 ,HLA_type1copyNum_withoutBAFBin,HLA_type1copyNum_withBAFBin
                 ,HLA_type2copyNum_withoutBAFBin,HLA_type2copyNum_withBAFBin
                 ,PVal,UnPairedPval,PVal_unique,UnPairedPval_unique,LossAllele,KeptAllele,numMisMatchSitesCov,propSupportiveSites) 
    HLAoutPut <- rbind(HLAoutPut,out)
    
    
  }
  dev.off()
  
  regionSpecOutPut <- cbind(region,HLAoutPut)
  PatientOutPut    <- rbind(PatientOutPut,regionSpecOutPut)
  
}  

if (gc.correction.step) {
  HLAoutLoc <- paste(workDir,full.patient,'.',minCoverageFilter,".DNA.GCcorrection.HLAlossPrediction.xls",sep="")
} else {
  HLAoutLoc <- paste(workDir,full.patient,'.',minCoverageFilter,".DNA.HLAlossPrediction.xls",sep="")
}
write.table(PatientOutPut,file=HLAoutLoc,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


if(performIntegerCopyNum)
{
  if (gc.correction.step) {
    HLABAFsummaryLoc <- paste(workDir,full.patient,'.', minCoverageFilter,".DNA.GCcorrection.IntegerCPN.xls",sep="")
  } else { 
    HLABAFsummaryLoc <- paste(workDir,full.patient,'.', minCoverageFilter,".DNA.IntegerCPN.xls",sep="")
  }
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




                                                                                                    
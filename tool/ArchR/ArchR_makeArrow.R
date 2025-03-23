print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
    library(stringr)
    library(argparse)
})

parse_columns <- function(list_str){
    if(!is.null(list_str)){
        list_toUse <- unlist(strsplit(list_str,","))
    } else{
        list_toUse <- c()
    }
    
    return(list_toUse)
}


parser <- ArgumentParser(description='Make ArchR project and Arrow Files')
parser$add_argument('input_files', metavar='input_files', help='input file of gz fragment filenames; tabbed sample name and filename')
parser$add_argument('--thread_num', metavar='thread_num', default=16, type='integer', help='number of threads to use; default 16')
parser$add_argument('--genome', metavar='genome', default='hg38', choices=c('hg38','hg19'), help='genome; options hg38, hg19; default hg38')
parser$add_argument('--exChr', metavar='exChr', default='chrY,chrM', help='exclude chromosomes; must be in form chr[1-22XYM]; ArchR default: chrY,chrM')
parser$add_argument('project_dir', metavar='project_dir', help='project (aka output) directory - full path')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('--extract_addToFileName', metavar='extract_addToFileName', help='if given, add to the extract filename')
parser$add_argument('--seed', metavar='seed', default=1, type='integer', help='randomization seed; default 1')
args <- parser$parse_args()

input_files <- args$input_files
thread_num <- args$thread_num
genome <- args$genome
exChr <- args$exChr
project_dir <- args$project_dir
wk_dir <- args$wk_dir
extract_dir <- args$extract_dir
extract_addToFileName <- args$extract_addToFileName
seed <- args$seed

print_time("Argument Checking")
if(!all(file.exists(input_files))) stop("Input file doesn't exist.")

excludeChr <- parse_columns(exChr)
if(!all(grepl('^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$|^chr[XYM]$',excludeChr))) stop('Exclusion chromosome formatting wrong.')
cat(paste('Excluding:',paste(collapse=', ',excludeChr),'\n'))

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')
if(substrRight(extract_dir,1)!='/') project_dir<-paste(sep="",extract_dir,'/')
if(!file.exists(extract_dir)) dir.create(extract_dir)

if(!is.null(extract_addToFileName)){
    outPrefix <- paste(sep='',extract_dir,basename(project_dir),'_',extract_addToFileName)
} else {
    outPrefix <- paste(sep='',extract_dir,basename(project_dir))
}
cat(paste("Using:",outPrefix,'\n'))

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_makeArrow_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Set wd")
if(length(wk_dir)==0){
    wk_dir <- paste(sep='',dirname(project_dir),'/')
} else {
    if(substrRight(wk_dir,1)!='/') wk_dir<-paste(sep="",wk_dir,'/')
    if(!file.exists(wk_dir)) wk_dir <- paste(sep='',dirname(project_dir),'/')
}
cat(paste('New working dir:',wk_dir,'\n'))
setwd(wk_dir)


print_time('Set seed and threads')
set.seed(seed)
addArchRThreads(threads = thread_num)

print_time("Data Organization")
inputFiles_tbl <- read.table(input_files,sep="\t",header=F,stringsAsFactors=F)
inputFiles <- inputFiles_tbl$V2
names(inputFiles) <- inputFiles_tbl$V1

addArchRGenome(genome)

print_time("Arrow Files")
#additional QC flags nullified to preserve QC already done
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  minFrags = 0,
  maxFrags = 1e+08,
  excludeChr = excludeChr,
  TileMatParams = list('excludeChr' = excludeChr),
  GeneScoreMatParams = list('excludeChr' = excludeChr)
)

ArrowFiles

print_time("ArchR Project")
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = basename(project_dir),
  copyArrows = TRUE
)

getAvailableMatrices(proj)

print_time("Saving ArchR Project")
proj <- saveArchRProject(ArchRProj = proj)

print_time("Extract Matrices")
saveRDS(getMatrixFromProject(ArchRProj = proj, useMatrix='TileMatrix', binarize = TRUE),
            paste(sep='',outPrefix,'_tile_sumExp.rds'))
saveRDS(getMatrixFromProject(ArchRProj = proj, useMatrix='GeneScoreMatrix', binarize = FALSE),
            paste(sep='',outPrefix,'_archrGAS_sumExp.rds'))

print_time("Saving ArchR Project")
proj <- saveArchRProject(ArchRProj = proj)

print_time("Done.")


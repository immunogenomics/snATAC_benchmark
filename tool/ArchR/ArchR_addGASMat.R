print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

parse_columns <- function(list_str){
    if(!is.null(list_str)){
        list_toUse <- unlist(strsplit(list_str,","))
    } else{
        list_toUse <- c()
    }
    
    return(list_toUse)
}


print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
    library(stringr)
    library(argparse)
})


parser <- ArgumentParser(description='Add Gene-like feature x cells to ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project (aka output) directory - full path')
parser$add_argument('sumExp_file', metavar='sumExp_file', help='SummarizedExperiment file - same cells in ArchR format must be there and the values must be integers!')
parser$add_argument('--exChr', metavar='exChr', default='chrY,chrM', help='exclude chromosomes; must be in form chr[1-22XYM]; ArchR default: chrY,chrM')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('--extract_addToFileName', metavar='extract_addToFileName', help='if given, add to the extract filename')
parser$add_argument('--thread_num', metavar='thread_num', default=8, type='integer', help='number of threads to use; default 8')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('--seed', metavar='seed', default=1, type='integer', help='randomization seed; default 1')
parser$add_argument('--noNameCheck', action='store_true', help='If given, do not exit if matrix name already given.') 
args <- parser$parse_args()

project_dir <- args$project_dir
sumExp_file <- args$sumExp_file
exChr <- args$exChr
extract_dir <- args$extract_dir
extract_addToFileName <- args$extract_addToFileName
thread_num <- args$thread_num
wk_dir <- args$wk_dir
seed <- args$seed
noNameCheck <- args$noNameCheck

print_time("Argument Checking")
if(!all(file.exists(project_dir,sumExp_file))) stop("Input file(s) don't exist.")

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
outFile <- paste(sep="",outPrefix,"_addGASMat_args.txt")
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

print_time("Load GAS SumExp")
gas_se <- readRDS(sumExp_file)

print_time("Load Project")
proj <- loadArchRProject(path=project_dir, showLogo=FALSE)
if(!noNameCheck & ("GeneExpressionMatrix" %in% getAvailableMatrices(proj))) stop('That matrix name is already used!')

print_time("Add Matrix")
getAvailableMatrices(proj)
proj <- addGeneExpressionMatrix(input=proj, seRNA = gas_se, excludeChr=excludeChr)
getAvailableMatrices(proj)

print_time("Extract Matrix")
extract_pxc <- getMatrixFromProject(ArchRProj = proj, useMatrix = "GeneExpressionMatrix", binarize = FALSE)
saveRDS(extract_pxc,paste(sep='',outPrefix,'_GeneExpressionMatrix.rds'))

print_time("Saving ArchR Project")
proj <- saveArchRProject(ArchRProj = proj)

print_time("Done.")


print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(Matrix.utils)
    library(argparse)
})

#COPIED from Signac::BinarizeCounts
makeBinary <- function(object) {
    if (inherits(x = object, what = "CsparseMatrix")) {
        slot(object = object, name = "x") <- rep.int(x = 1, times = length(x = slot(object = object, name = "x")))
    } else {
        object[object > 1] <- 1
    }
    return(object)
}


parser <- ArgumentParser(description='Concatenate matrices')
parser$add_argument('listing_file', metavar='listing_file', help='listing file of matrix files to concatenate')
parser$add_argument('by', metavar='by', choices=c('row','col'), help='row: rownames should be the same and concat columns; col: rownames should be the same and concat columns.')
parser$add_argument('--binarize', action='store_true', help='if given, also save binary matrix.')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

listing_file <- args$listing_file
by <- args$by
binarize <- args$binarize
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(listing_file))) stop("Input file doesn't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_catMat_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}

print_time("Read file list")
file_list <- read.table(listing_file,stringsAsFactors=F)
file_list <- file_list$V1

print_time("Starting concatenation")
for(i in 1:length(file_list)){
    if(!file.exists(file_list[i])) stop(paste("ERROR: File doesn't exist:",file_list[i]))

    if(i==1){
        allCat_df <- readRDS(file_list[i])
        cat(paste(sep="","initial dimensions of index ",i,", row = ",nrow(allCat_df)," col = ",ncol(allCat_df),"\n"))
    } else {
        placeholder <- readRDS(file_list[i])
        
        if(by=='row'){
            if(identical(rownames(allCat_df),rownames(placeholder))){
                allCat_df <- cbind(allCat_df,placeholder)
                cat(paste(sep="","After adding matrix ",i,", new # col= ",ncol(allCat_df),"\n"))
            } else {
                cat(paste("File",i,"rownames are different\n"))
            }
        } else if(by=='col'){
            if(identical(colnames(allCat_df),colnames(placeholder))){
                allCat_df <- rbind(allCat_df,placeholder)
                cat(paste(sep="","After adding matrix ",i,", new # row= ",nrow(allCat_df),"\n"))
            } else {
                cat(paste("File",i,"colnames are different\n"))
            }
        } else {
            stop('by argument should only be row or col')
        }
        
        rm(placeholder)
    }
}
print_time("Stopping concatenation")

cat(paste(sep="","final dimensions: row = ",nrow(allCat_df)," col = ",ncol(allCat_df),"\n"))

saveRDS(allCat_df,paste(sep='',outPrefix,'.rds'))

if(binarize){
    allCat_df <- makeBinary(allCat_df)
    saveRDS(allCat_df,paste(sep='',outPrefix,'_binary.rds'))
}

print_time("Done.")


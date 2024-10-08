suppressMessages({
    library(Matrix)
    library(gtools)
    library(Rmisc)
    library(presto)
    library(ggpubr)
    
    library(plyr)
    library(stringr)
    library(tidyr)
    
    library(ggplot2)
    library(ggrastr)
    library(ggrepel)
    library(viridis)
    library(scales)
    library(RColorBrewer)
    library(grid)
    library(gridExtra)
    library(repr)
})


##UTILITY FUNCTIONS

#SORTING - to allow for the hyphen to be used as a delimiter or a negative number.
#modified from gtools
mixedorder_new <- function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE,
    numeric.type = c("decimal", "roman"), roman.case = c("upper",
        "lower", "both"),keepNegative=FALSE)
{
    numeric.type <- match.arg(numeric.type)
    roman.case <- match.arg(roman.case)
    if (length(x) < 1)
        return(NULL)
    else if (length(x) == 1)
        return(1)
    if (!is.character(x))
        return(order(x, decreasing = decreasing, na.last = na.last))
    delim = "\\$\\@\\$"
    if (numeric.type == "decimal") {
        if(keepNegative)
            regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))"
        else
            regex <- "((?:(?i)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)))"
        numeric <- function(x) as.numeric(x)
    }
    else if (numeric.type == "roman") {
        regex <- switch(roman.case, both = "([IVXCLDMivxcldm]+)",
            upper = "([IVXCLDM]+)", lower = "([ivxcldm]+)")
        numeric <- function(x) roman2int(x)
    }
    else stop("Unknown value for numeric.type: ", numeric.type)
    nonnumeric <- function(x) {
        ifelse(is.na(numeric(x)), toupper(x), NA)
    }
    x <- as.character(x)
    which.nas <- which(is.na(x))
    which.blanks <- which(x == "")
    delimited <- gsub(regex, paste(delim, "\\1", delim, sep = ""),
        x, perl = TRUE)
    step1 <- strsplit(delimited, delim)
    step1 <- lapply(step1, function(x) {x[x > ""]})
    suppressWarnings(step1.numeric <- lapply(step1, numeric))
    suppressWarnings(step1.character <- lapply(step1, nonnumeric))
    maxelem <- max(sapply(step1, length))
    step1.numeric.t <- lapply(1:maxelem, function(i) {sapply(step1.numeric,
        function(x) {x[i]})})
    step1.character.t <- lapply(1:maxelem, function(i) {sapply(step1.character,
        function(x) {x[i]})})
    rank.numeric <- sapply(step1.numeric.t, rank)
    rank.character <- sapply(step1.character.t, function(x) {as.numeric(factor(x))})
    rank.numeric[!is.na(rank.character)] <- 0
    rank.character <- t(t(rank.character) + apply(matrix(rank.numeric),
        2, max, na.rm = TRUE))
    rank.overall <- ifelse(is.na(rank.character), rank.numeric,
        rank.character)
    order.frame <- as.data.frame(rank.overall)
    if (length(which.nas) > 0)
        if (is.na(na.last))
            order.frame[which.nas, ] <- NA
        else if (na.last)
            order.frame[which.nas, ] <- Inf
        else order.frame[which.nas, ] <- -Inf
    if (length(which.blanks) > 0)
        if (is.na(blank.last))
            order.frame[which.blanks, ] <- NA
        else if (blank.last)
            order.frame[which.blanks, ] <- 1e+99
        else order.frame[which.blanks, ] <- -1e+99
    order.frame <- as.list(order.frame)
    order.frame$decreasing <- decreasing
    order.frame$na.last <- NA
    retval <- do.call("order", order.frame)
    return(retval)
}

#SORTING
#modified from gtools
mixedsort_new <- function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE,
    numeric.type = c("decimal", "roman"), roman.case = c("upper",
        "lower", "both"),keepNegative=FALSE)
{
    ord <- mixedorder_new(x, decreasing = decreasing, na.last = na.last,
        blank.last = blank.last, numeric.type = numeric.type,
        roman.case = roman.case, keepNegative = keepNegative)
    x[ord]
}


##MARKER PEAK/GENE FUNCTIONS

#features by hex on UMAP
plot_markerPeaks_norm_hex_v2 <- function(thisMeta, thisGxC, xCol, yCol, plotCol=5, titleSize=30, cutCap=0.01,
                                      colorOpt='viridis',hex_bins=30,exclScale=TRUE,scaleLab=NA,
                                      plot_genes=c('CD3D','GNLY','MS4A1','XBP1','C1QA',
                                                   'VWF','PDPN','PRG4','THY1','NOTCH3'), titleFace='italic'){

    if(!identical(rownames(thisMeta),colnames(thisGxC))) stop("cells have to be ordered the same.")

    toPlot <- thisMeta[1:nrow(thisMeta),1:ncol(thisMeta)]

    myplots <- list()
    for (i in 1:length(plot_genes)) {
        gene <- plot_genes[i]
        
        if(cutCap!=0){
            max.cutoff = quantile(thisGxC[gene, ], 1-cutCap)
            min.cutoff = quantile(thisGxC[gene, ], cutCap)
            tmp <- sapply(X = thisGxC[gene, ], FUN = function(x) {
                return(ifelse(test = x > max.cutoff, yes = max.cutoff, 
                    no = x))
            })
            tmp <- sapply(X = tmp, FUN = function(x) {
                return(ifelse(test = x < min.cutoff, yes = min.cutoff, 
                    no = x))
            })
            toPlot$gene <- as.numeric(tmp)
        } else {
            toPlot$gene <- thisGxC[gene, ]
        }

        ind <- paste("p", i, sep = "")
        ind <- ggplot(data = toPlot,aes(x = !!sym(xCol), y = !!sym(yCol), z = gene)) +
          stat_summary_hex(bins=hex_bins) + scale_fill_viridis(option = colorOpt) +
          labs(x="", y="") + theme_bw(base_size = 15) +
          theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),
            plot.title = element_text(color="black", size=titleSize, face = titleFace)) + labs(title = gene)
        if(exclScale){
            ind <- ind + theme(legend.position = "none")
        } else{
            if(!is.na(scaleLab)) ind <- ind + labs(fill=scaleLab)
        }
        myplots[[i]] <- ind

    }

    p <- arrangeGrob(grobs = myplots, ncol=plotCol)
    
    return(p)

}


scaleGene_forHeatmap <- function(gOrd,ctOrd,gCTnorm,geneCTcutoff=0.15){
    
    if(!all(gOrd %in% rownames(gCTnorm))) stop('genes from order must be in normalized gene expression')
    if(!all(ctOrd %in% colnames(gCTnorm))) stop('cell types from order must be in normalized gene expression')
    
    ll <- apply(gCTnorm[gOrd,],1,max)
    gOrd <- names(ll[ll>=geneCTcutoff])
    #print(cat(paste(length(ll),'genes subset to',length(gOrd),'\n')))

    gCTnorm_subset <- gCTnorm[rev(gOrd),rev(ctOrd)]
        
    gCTnorm_subset_scaled <- t(scale(t(gCTnorm_subset)))
    
    return(gCTnorm_subset_scaled)
    
}


#feature heatmaps
#assumes that the given row/col is the order you want!
pseudobulk_scaled_heatmap <- function(toPlot,xlab,ylab,fillLab,plotTit=NULL,scale_lim=NA,clustColors=NA,includeZero=FALSE,
                                      colorLow='blue',colorMid='white',colorHigh='red'){
    
    toGather <- as.data.frame(toPlot)
    cols_toGather <- colnames(toGather)
    toGather$peak <- rownames(toGather)
    
    gathered <- gather(toGather,'CT','agg',all_of(cols_toGather))
    gathered$peak <- factor(gathered$peak,levels=rev(rownames(toGather)))
    gathered$CT <- factor(gathered$CT,levels=cols_toGather)
    
    if(!is.na(scale_lim)){
        scale_limits <- c(-scale_lim,scale_lim)
    } else {
        if(min(gathered$agg)>0 & includeZero){
            ll <- 0
        } else {
            ll <- min(gathered$agg)
        }
        scale_limits <- c(ll,max(gathered$agg))
    }
    
    g <- ggplot(gathered,aes(x=peak,y=CT,fill=agg)) + geom_tile() +
            theme_classic(base_size=20) + labs(x=xlab,y=ylab,fill=fillLab) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + 
            scale_fill_gradient2(high=colorHigh,mid=colorMid,low=colorLow,midpoint=0,limits=scale_limits)
    if(all(levels(gathered$CT) %in% names(clustColors))){
        suppressWarnings(g <- g + theme(axis.text.y = element_text(color=clustColors[levels(gathered$CT)],face='bold',size=20)))
    }
    if(!is.null(plotTit)) {g <- g + ggtitle(plotTit)}
    
    return(g)
    
}



##ODDS RATIO FUNCTIONS

#OR PLOT ORDER - plot classes/states in order
plot_order <- function(uniq_val){
    if(all(uniq_val %in% grep('[a-zA-Z]+-[0-9]+:',uniq_val,value=TRUE)))
        ret <- mixedsort_new(uniq_val)
    else if(!any(is.na(as.numeric(uniq_val))))
        ret <- as.character(sort(as.numeric(uniq_val)))
    else
        ret <- sort(uniq_val)
    
    return(ret)
}

#OR PLOT ORDER - plot classes/states such that the diagonal has the highest OR
reorder_col_diag_plotOR <- function(fish_df,xCol,yCol,yOrd=NA,mCol='lnOR',op='max'){
    if(!(op %in% c('min','max'))) op='max' #defaulting to max.
    
    if(!any(is.na(yOrd))){fish_df[,yCol] <- factor(fish_df[,yCol],levels=rev(yOrd))}
    
    row_col_maxOR_df <- data.frame('row'=character(),'col'=character(),'val'=numeric(),stringsAsFactors=FALSE)
    for(thisX in unique(fish_df[,xCol])){
        this_subset <- fish_df[which(fish_df[,xCol]==thisX),]
        if(op=='min'){
            this_idx <- which.min(this_subset[,mCol])
        } else {
            this_idx <- which.max(this_subset[,mCol])
        }
        row_col_maxOR_df <- rbind(row_col_maxOR_df,data.frame('row'=this_subset[this_idx,yCol],
                                                              'col'=thisX,'val'=this_subset[this_idx,mCol],
                                                              stringsAsFactors=FALSE))
    }
    #reverse row order since that's what the plot does!
    row_col_maxOR_df <- row_col_maxOR_df[order(row_col_maxOR_df$row,decreasing=TRUE),]
    new_xOrd <- row_col_maxOR_df$col
    return(new_xOrd)
}

#For plotting purposes
label_spacing <- function(s,wiggle=1){
    if(substr(s,nchar(s)-wiggle,nchar(s)-wiggle)=='\n'){substr(s,nchar(s)-wiggle,nchar(s)-wiggle) <- ' '}
    return(s)
}

#For plotting purposes
fix_infinite <- function(this_vec,addVal=1,multVal=1){
    nonInf <- this_vec[which(!is.infinite(this_vec))]
    limit <- ceiling(max(nonInf,abs(min(nonInf)))*multVal)+addVal
    this_vec[which(is.infinite(this_vec) & this_vec>0)] <- limit
    this_vec[which(is.infinite(this_vec) & this_vec<0)] <- -limit
    
    return(this_vec)
}

#OR calculation
calc_OR <- function(toPlot, xCol, yCol, pLim = 0.05,includeNomP = FALSE){
    if(!all(c(xCol,yCol) %in% colnames(toPlot))) stop('xCol and yCol need to be in toPlot columns')
    
    fisher_res_df <- data.frame('xCol'=character(),'yCol'=character(),
                                'inX_inY'=integer(),'inX_noY'=integer(),'noX_inY'=integer(),'noX_noY'=integer(),
                                'pval'=numeric(),'OR'=numeric(),'nullOR'=numeric(),
                                'CI_low'=numeric(),'CI_high'=numeric(),stringsAsFactors=FALSE)
    colnames(fisher_res_df)[1] <- xCol
    colnames(fisher_res_df)[2] <- yCol
    for(xVal in sort(unique(toPlot[,xCol]))){
        for(yVal in sort(unique(toPlot[,yCol]))){
            inX_inY <- nrow(toPlot[which(toPlot[,xCol]==xVal & toPlot[,yCol]==yVal),])
            inX_noY <- nrow(toPlot[which(toPlot[,xCol]==xVal & toPlot[,yCol]!=yVal),])
            noX_inY <- nrow(toPlot[which(toPlot[,xCol]!=xVal & toPlot[,yCol]==yVal),])
            noX_noY <- nrow(toPlot[which(toPlot[,xCol]!=xVal & toPlot[,yCol]!=yVal),])

            fish_df <- data.frame('inX'=c(inX_inY,inX_noY),'noX'=c(noX_inY,noX_noY),
                                  row.names=c("inY", "noY"),
                                  stringsAsFactors = FALSE)

            fish_test <- fisher.test(fish_df)

            placeholder <- data.frame('xCol'=xVal,'yCol'=yVal,
                                      'inX_inY'=inX_inY,'inX_noY'=inX_noY,
                                      'noX_inY'=noX_inY,'noX_noY'=noX_noY,
                                      'pval'=fish_test$p.value,'OR'=fish_test$estimate,
                                      'nullOR'=fish_test$null.value,'CI_low'=fish_test$conf.int[1],
                                      'CI_high'=fish_test$conf.int[2],stringsAsFactors=FALSE)
            colnames(placeholder)[1] <- xCol
            colnames(placeholder)[2] <- yCol

            fisher_res_df <- rbind(fisher_res_df,placeholder)
        }
    }
    rownames(fisher_res_df) <- NULL
    head(fisher_res_df)

    fisher_res_df$sig <- ifelse(fisher_res_df$pval<=pLim,TRUE,FALSE)
    fisher_res_df$padj <- p.adjust(fisher_res_df$pval,method='BH')
    fisher_res_df$sigBH <- ifelse(fisher_res_df$padj<=pLim,TRUE,FALSE)
    fisher_res_df$signif <- 'not'
    if(includeNomP){
        fisher_res_df[which(fisher_res_df$sigBH==FALSE & fisher_res_df$sig==TRUE),'signif'] <- 'nominal'
    }
    fisher_res_df[which(fisher_res_df$sigBH==TRUE),'signif'] <- 'BH'

    fisher_res_df$lnOR <- fix_infinite(log(fisher_res_df$OR))
    
    return(fisher_res_df)

}

#OR - calculate and plot OR values between two columns of a dataframe
plot_OR <- function(toPlot, xCol, yCol, xLab, yLab, xOrder, yOrder, xAxisVert = FALSE, fCol='lnOR', flab='ln(OR)',clustColors=NA,
                    wiggle=1,yLab_charLim=50){
    
    if(!all(c(xCol,yCol,fCol) %in% colnames(toPlot))) stop('xCol and yCol need to be in toPlot columns')
    
    if(!all(unique(toPlot[,xCol]) %in% xOrder)){
        print('xOrder not sufficient, just sorting')
        xOrder <- sort(unique(toPlot[,xCol]))
    } else {
        xOrder <- xOrder[which(xOrder %in% unique(toPlot[,xCol]))]
    }
    
    if(!all(unique(toPlot[,yCol]) %in% yOrder)){
        print('yOrder not sufficient, just sorting')
        yOrder <- sort(unique(toPlot[,yCol]))
    } else {
        yOrder <- yOrder[which(yOrder %in% unique(toPlot[,yCol]))]
    }
    
    toPlot[,yCol] <- factor(toPlot[,yCol],levels=rev(yOrder))
    toPlot[,xCol] <- factor(toPlot[,xCol],levels=xOrder)
    
    xLabels_wSpaces <- str_replace_all(unique(toPlot[,xCol]),' ','\n')
    xLabels_wSpaces <- lapply(xLabels_wSpaces,FUN=label_spacing,wiggle=wiggle)
    names(xLabels_wSpaces) <- unique(toPlot[,xCol])
    
    yLabels_wSpaces <- c()
    for(yy in sort(unique(toPlot[,yCol]))){
        if(nchar(yy)>yLab_charLim){
            yLabels_wSpaces <- c(yLabels_wSpaces,replace_space_newline_afterHalf(yy,wiggle=wiggle))
        } else {
            yLabels_wSpaces <- c(yLabels_wSpaces,yy)
        }
    }
    names(yLabels_wSpaces) <- sort(unique(toPlot[,yCol]))
    
    
    g <- ggplot(toPlot,aes(x=!!sym(xCol),y=!!sym(yCol),fill=!!sym(fCol))) + 
            geom_tile() + theme_classic(base_size=20) + 
            scale_fill_gradient2(low = "grey85", mid = "white", high = "red", midpoint = 0, na.value='grey85') +
            labs(x=xLab,y=yLab,fill=flab) +
            geom_tile(data=toPlot[which(toPlot$signif=='not'),],fill='white') + 
            scale_x_discrete(labels=xLabels_wSpaces[levels(toPlot[,xCol])]) +
            scale_y_discrete(labels=yLabels_wSpaces[levels(toPlot[,yCol])])
    if(all(c(levels(toPlot[,xCol]),levels(toPlot[,yCol])) %in% names(clustColors))){
        suppressWarnings(g <- g + theme(axis.text.x = element_text(color=clustColors[levels(toPlot[,xCol])],
                                                                   face='bold',size=17)) + 
                                  theme(axis.text.y = element_text(color=clustColors[levels(toPlot[,yCol])],
                                                                   face='bold',size=17)))
    }
    
    return(g)
}


## QC

get_2D_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix,iy)
    return(dens$z[ii])
}



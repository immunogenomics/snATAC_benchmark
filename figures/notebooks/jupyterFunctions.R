suppressMessages({
    library(Matrix)
    library(gtools)
    library(Rmisc)
    library(presto)
    library(ggpubr)
    
    library(plyr)
    library(dplyr)
    library(stringr)
    library(tidyr)
    
    library(ggplot2)
    library(ggrastr)
    library(ggrepel)
    library(viridis)
    library(scales)
    library(RColorBrewer)
    library(patchwork)
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



## Linear Models

#to source in jupyter notebook

#NOTE that this does not exclude one column to prevent overdetermination!
create1HotMat <- function(this_meta,this_col){
    if(!is.character(this_meta[,this_col])) {this_meta[,this_col] <- as.character(this_meta[,this_col])}
    donors <- unique(this_meta[,this_col])
    for(idx in 1:length(donors)){
        thisDonor <- donors[idx]
        thisDonorCol <- as.integer(ifelse(this_meta[,this_col]==thisDonor,1,0))

        if(idx==1){
            this_donor_df <- data.frame('donor'=thisDonorCol)
            rownames(this_donor_df) <- rownames(this_meta)
            colnames(this_donor_df) <- thisDonor
        } else {
            this_donor_df[,thisDonor] <- thisDonorCol
        }
    }
    ll <- colSums(this_donor_df)
    rr <- table(this_meta[,this_col])
    ss <- as.numeric(rr)
    names(ss) <- names(rr)
    theseNames <- sort(names(rr))
    if(!identical(ll[theseNames],ss[theseNames])){
        print(ll[theseNames])
        print(ss[theseNames])
        stop("1 hot issues")
    }

    return(this_donor_df)
}

linearModel_both <- function(this_df, ds_col='dataset', ft_col='feature', md_col='method', hm_col='Harmony',
                             ds_ref='dataset1', ft_ref='tile', md_ref='LSI', hm_ref='noHarmony', hm_binary=TRUE, 
                             bio_col='CT_mean', batch_col='sample_mean', bio_perc=0.60, batch_perc=0.40,
                             bio_lowScore_good=TRUE, batch_lowScore_good=TRUE, glm_family='gaussian', verbose=FALSE, 
                             all_levels=c('peak','cCRE','tile','basicGAS','archrGAS',
                                          'PCA','LSI','rLSI','rCCA','itLSI','SA2','CS','pVI','noHarmony','Harmony',
                                          'dataset1','dataset2','dataset3','dataset4','dataset5')){
    
    if(!all(c(ds_col,ft_col,md_col,hm_col,bio_col,batch_col) %in% colnames(this_df))) stop('Not all required columns present.')
    if(!(ds_ref %in% unique(this_df[,ds_col]))) stop('dataset reference not in dataset values')
    if(!(ft_ref %in% unique(this_df[,ft_col]))) stop('feature reference not in feature values')
    if(!(md_ref %in% unique(this_df[,md_col]))) stop('method reference not in method values')
    if(hm_binary) {this_df[,hm_col] <- ifelse(this_df[,hm_col],hm_col,paste(sep='','no',hm_col))}
    if(!(hm_ref %in% unique(this_df[,hm_col]))) stop('correction reference not in correction values')
    
    suppressMessages({
        library(stringr)
        library(dplyr)
        library(stats)
    })
    
    forModel_df <- this_df[,c(ds_col,md_col,ft_col,hm_col)]
    if(verbose) {print(head(forModel_df))}
    forModel_df[,ds_col] <- as.character(forModel_df[,ds_col])
    forModel_df[,ft_col] <- as.character(forModel_df[,ft_col])
    forModel_df[,md_col] <- as.character(forModel_df[,md_col])
    forModel_df[,hm_col] <- as.character(forModel_df[,hm_col])
    forModel_df$bio <- this_df[,bio_col]
    forModel_df$batch <- this_df[,batch_col]
    if(verbose) {print(head(forModel_df))}
    
    forModel_df <- forModel_df %>% 
        group_by(!!sym(ds_col)) %>%
        mutate(mean_bio_rank = rank(if(bio_lowScore_good) {bio} else {-bio},na.last=NA), 
               mean_batch_rank = rank(if(batch_lowScore_good) {batch} else {-batch},na.last=NA)) %>% 
        ungroup() %>% as.data.frame()
    
    forModel_df$overall_score <- forModel_df$mean_bio_rank * bio_perc + forModel_df$mean_batch_rank * batch_perc
    if(verbose) {print(head(forModel_df))}
    
    ds <- create1HotMat(forModel_df,ds_col)
    ft <- create1HotMat(forModel_df,ft_col)
    md <- create1HotMat(forModel_df,md_col)
    hm <- create1HotMat(forModel_df,hm_col)
    
    #to avoid over-determining
    ds <- as.matrix(ds[,which(colnames(ds)!=ds_ref)])
    ft <- as.matrix(ft[,which(colnames(ft)!=ft_ref)])
    md <- as.matrix(md[,which(colnames(md)!=md_ref)])
    if(hm_binary){
        hm <- as.vector(hm[,which(colnames(hm)!=hm_ref)])
    } else {
        hm <- as.matrix(hm[,which(colnames(hm)!=hm_ref)])
    }
    
    overall_score <- as.vector(forModel_df[,'overall_score'])
    
    m1 <- glm(overall_score ~ ds + ft + md + hm, family=glm_family)
    print(m1)
    m1_coeff <- as.data.frame(summary(m1)$coefficient)
    
    if(all(all_levels %in% c(unique(forModel_df[,ds_col]),unique(forModel_df[,md_col]),
                              unique(forModel_df[,ft_col]),unique(forModel_df[,hm_col])))){
        
        varNames <- substr(rownames(m1_coeff),3,nchar(rownames(m1_coeff)))
        varNames[1] <- 'Intercept'
        varNames[length(varNames)] <- 'Harmony'
        if(verbose) {print(varNames)}
        
        subset_levels <- c(all_levels[which(all_levels %in% varNames)],'Intercept')
        
        m1_coeff$varName <- factor(varNames,level=rev(subset_levels))
        colnames(m1_coeff) <- c('Estimate','stdErr','tVal','tSig','varName')
        if(verbose) {print(head(m1_coeff))}
    }
    
    return(list('rank'=forModel_df,'coeff'=m1_coeff))
}


linearModel_sep <- function(this_df, ds_col='dataset', ft_col='feature', md_col='method', hm_col='Harmony',
                            ds_ref='dataset1', ft_ref='tile', md_ref='LSI', hm_ref='noHarmony', hm_binary=TRUE, 
                            bio_col='CT_mean', batch_col='sample_mean',
                            bio_lowScore_good=TRUE, batch_lowScore_good=TRUE, glm_family='gaussian', verbose=FALSE, 
                            all_levels=c('peak','cCRE','tile','basicGAS','archrGAS',
                                         'PCA','LSI','rLSI','rCCA','itLSI','SA2','CS','pVI','noHarmony','Harmony',
                                         'dataset1','dataset2','dataset3','dataset4','dataset5')){
    
    if(!all(c(ds_col,ft_col,md_col,hm_col,bio_col,batch_col) %in% colnames(this_df))) stop('Not all required columns present.')
    if(!(ds_ref %in% unique(this_df[,ds_col]))) stop('dataset reference not in dataset values')
    if(!(ft_ref %in% unique(this_df[,ft_col]))) stop('feature reference not in feature values')
    if(!(md_ref %in% unique(this_df[,md_col]))) stop('method reference not in method values')
    if(hm_binary) {this_df[,hm_col] <- ifelse(this_df[,hm_col],hm_col,paste(sep='','no',hm_col))}
    if(!(hm_ref %in% unique(this_df[,hm_col]))) stop('correction reference not in correction values')
    
    suppressMessages({
        library(stringr)
        library(dplyr)
        library(stats)
    })
    
    forModel_df <- this_df[,c(ds_col,md_col,ft_col,hm_col)]
    if(verbose) {print(head(forModel_df))}
    forModel_df[,ds_col] <- as.character(forModel_df[,ds_col])
    forModel_df[,ft_col] <- as.character(forModel_df[,ft_col])
    forModel_df[,md_col] <- as.character(forModel_df[,md_col])
    forModel_df[,hm_col] <- as.character(forModel_df[,hm_col])
    forModel_df$bio <- this_df[,bio_col]
    forModel_df$batch <- this_df[,batch_col]
    if(verbose) {print(head(forModel_df))}
    
    forModel_df <- forModel_df %>% 
        group_by(!!sym(ds_col)) %>%
        mutate(mean_bio_rank = rank(if(bio_lowScore_good) {bio} else {-bio},na.last=NA), 
               mean_batch_rank = rank(if(batch_lowScore_good) {batch} else {-batch},na.last=NA)) %>% 
        ungroup() %>% as.data.frame()
    if(verbose) {print(head(forModel_df))}
    
    ds <- create1HotMat(forModel_df,ds_col)
    ft <- create1HotMat(forModel_df,ft_col)
    md <- create1HotMat(forModel_df,md_col)
    hm <- create1HotMat(forModel_df,hm_col)
    
    #to avoid over-determining
    ds <- as.matrix(ds[,which(colnames(ds)!=ds_ref)])
    ft <- as.matrix(ft[,which(colnames(ft)!=ft_ref)])
    md <- as.matrix(md[,which(colnames(md)!=md_ref)])
    if(hm_binary){
        hm <- as.vector(hm[,which(colnames(hm)!=hm_ref)])
    } else {
        hm <- as.matrix(hm[,which(colnames(hm)!=hm_ref)])
    }
    
    bio_score <- as.vector(forModel_df$mean_bio_rank)
    batch_score <- as.vector(forModel_df$mean_batch_rank)
    
    bio_mod <- glm(bio_score ~ ds + ft + md + hm, family=glm_family)
    print(bio_mod)
    batch_mod <- glm(batch_score ~ ds + ft + md + hm, family=glm_family)
    print(batch_mod)
    
    bio_coeff <- as.data.frame(summary(bio_mod)$coefficient)
    batch_coeff <- as.data.frame(summary(batch_mod)$coefficient)
    
    if(!identical(rownames(bio_coeff),rownames(batch_coeff))) stop('coefficient dfs do not align')
    
    sepMod_coeff <- batch_coeff[1:nrow(batch_coeff),1:ncol(batch_coeff)]
    colnames(sepMod_coeff) <- paste(sep='_','batch',c('Estimate','stdErr','tVal','tSig'))
    sepMod_coeff <- cbind(sepMod_coeff,bio_coeff)
    colnames(sepMod_coeff) <- c(paste(sep='_','batch',c('Estimate','stdErr','tVal','tSig')),
                                paste(sep='_','bio',c('Estimate','stdErr','tVal','tSig')))
    if(verbose) {print(head(sepMod_coeff))}
    
    if(all(all_levels %in% c(unique(forModel_df[,ds_col]),unique(forModel_df[,md_col]),
                              unique(forModel_df[,ft_col]),unique(forModel_df[,hm_col])))){
        
        varNames <- substr(rownames(sepMod_coeff),3,nchar(rownames(sepMod_coeff)))
        varNames[1] <- 'Intercept'
        varNames[length(varNames)] <- 'Harmony'
        if(verbose) {print(varNames)}
        
        subset_levels <- c(all_levels[which(all_levels %in% varNames)],'Intercept')
        
        sepMod_coeff$varName <- factor(varNames,level=rev(subset_levels))
        sepMod_coeff$varType = factor(c('Model',rep('Dataset',4),rep('Feature',4),rep('Method',7),'Correction'),
                                      levels=c('Dataset','Feature','Method','Correction','Model'))
        if(verbose) {print(head(sepMod_coeff))}
    }
    
    return(list('rank'=forModel_df,'coeff'=sepMod_coeff))
}




## aggrPlots

makeCompDF_wNA <- function(thisDF,ds,catCol,intVal,ctVal,valCol,groupCol,dsCol='dataset',verbose=FALSE,rmNA=TRUE){
    if(!all(c(dsCol,catCol,valCol,groupCol) %in% colnames(thisDF))) stop('Not all required columns in this dataframe')
    if(!all(c(intVal,ctVal) %in% unique(thisDF[,catCol]))) stop('Not all required values in dataframe cat column.')

    toPlot <- thisDF[which(thisDF[,dsCol]==ds & (thisDF[,catCol] %in% c(intVal,ctVal))),]
    toPlot <- tidyr::spread(toPlot, !!sym(catCol), !!sym(valCol))

    summary_sample <- Rmisc::summarySE(toPlot, intVal, groupvars = groupCol, na.rm=rmNA)
    colnames(summary_sample)[(length(groupCol)+1):ncol(summary_sample)] <- c('sample_N','sample_mean',
                                                                             'sample_sd','sample_se','sample_ci')
    if(verbose){print(head(summary_sample))}

    summary_CT <- Rmisc::summarySE(toPlot, ctVal, groupvars = groupCol, na.rm=rmNA)
    colnames(summary_CT)[(length(groupCol)+1):ncol(summary_CT)] <- c('CT_N','CT_mean','CT_sd','CT_se','CT_ci')
    if(verbose){print(head(summary_CT))}

    if(!identical(summary_sample[,1:(length(groupCol))],summary_CT[,1:(length(groupCol))])) stop('not identical columns')
    toPlot <- cbind(summary_sample,summary_CT[,(length(groupCol)+1):ncol(summary_CT)])
    if(verbose){print(head(toPlot))}

    return(toPlot)
}

makeCompDF_wNA_woTR <- function(thisDF,sampCol,phenoCol,groupCol,verbose=FALSE,rmNA=TRUE){
    if(!all(c(sampCol,phenoCol,groupCol) %in% colnames(thisDF))) stop('Not all required columns in this dataframe')

    summary_sample <- Rmisc::summarySE(thisDF, sampCol, groupvars = groupCol, na.rm=rmNA)
    colnames(summary_sample)[(length(groupCol)+1):ncol(summary_sample)] <- c('sample_N','sample_mean',
                                                                             'sample_sd','sample_se','sample_ci')
    if(verbose){print(head(summary_sample))}

    summary_CT <- Rmisc::summarySE(thisDF, phenoCol, groupvars = groupCol, na.rm=rmNA)
    colnames(summary_CT)[(length(groupCol)+1):ncol(summary_CT)] <- c('CT_N','CT_mean','CT_sd','CT_se','CT_ci')
    if(verbose){print(head(summary_CT))}

    if(!identical(summary_sample[,1:(length(groupCol))],summary_CT[,1:(length(groupCol))])) stop('not identical columns')
    toPlot <- cbind(summary_sample,summary_CT[,(length(groupCol)+1):ncol(summary_CT)])
    if(verbose){print(head(toPlot))}

    return(toPlot)
}



## metric math

getDelta <- function(ll,keyCol='Harmony',valCol='NNmet_val'){
    if(!all(c(keyCol,valCol) %in% colnames(ll))) stop('colname issue')
    
    spread_df <- tidyr::spread(ll, keyCol, valCol)
    colnames(spread_df)[(ncol(spread_df)-1):ncol(spread_df)] <- paste(sep='_',keyCol,colnames(spread_df)[(ncol(spread_df)-1):ncol(spread_df)])
    return(spread_df)
}

minMax_scale <- function(x, na.rm = TRUE) {
    return((x- min(x,na.rm=na.rm)) /(max(x,na.rm=na.rm)-min(x,na.rm=na.rm)))
}

missPerc_funct <- function(xx,nn_tot=200){
    return((nn_tot-xx)/nn_tot)
}

stdErr_calc <- function(x, na.rm=TRUE) {
  if(na.rm) {x <- na.omit(x)}
  return(sqrt(var(x)/length(x)))
}



## UMAP from files

substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

getCol_fromFile <- function(df1,fl2,col_get,col_conv,fl2_str=''){
    if(!file.exists(fl2)) stop(paste(fl2_str,'file does not exist'))
    
    for(xx in 1:length(col_conv)){ col_get <- sub(names(col_conv[xx]),unname(col_conv[xx]),col_get) }

    df2 <- readRDS(fl2)
    if(!all(col_get %in% colnames(df2))) stop('meta column name issue')

    if(!identical(rownames(df1),rownames(df2))){
        if(!identical(sort(rownames(df1)),sort(rownames(df2)))){
            stop('meta cell name issue')
        } else {
            df2 <- df2[rownames(df1),]
        }
    }

    df2 <- df2[,col_get,drop=FALSE]
    for(xx in 1:length(col_conv)){ colnames(df2) <- sub(unname(col_conv[xx]),names(col_conv[xx]),colnames(df2)) }
    
    return(df2)
    
}

getUMAP_fromFile <- function(id,od,dataset,feature,method,harmony,cr,verbose=FALSE){
    if(!all(file.exists(c(od,id)))) stop('input or output directory does not exist.')
    if(substrRight(od,1)!='/') od<-paste(sep='',od,'/')
    if(substrRight(id,1)!='/') outDir<-paste(sep='',id,'/')
    
    if(!any(c('Harmony','noHarmony') %in% harmony)) stop('Harmony input only takes in: Harmony, noHarmony')
    
    meta_col <- c('sample','nFrags','ReadsInPromoter','nUMI','nGene','MTperc','ori_phenotype','ori_phenotype_abbr',
                  'phenotype','phenotype_abbr','RNA_cLISI','RNA_iLISI')
    lisi_col <- c('ATAC_cLISI','ATAC_iLISI')
    PH_col <- c('ATAC_PH200_mKNN','ATAC_PH200_lsKLD')
    ST_col <- c('ATAC_ST200_mKNN','ATAC_ST200_lsKLD')
    if(!all(cr %in% c(meta_col,lisi_col,PH_col,ST_col))) stop('color keyword not recognized')
    
    first <- TRUE
    for(ds in dataset){
        for(ft in feature){
            for(md in method){
                for(hm in harmony){
                    if(hm=='Harmony') {hm_name <- 'Harmony_'; hm <- TRUE} else {hm_name <- ''; hm <- FALSE}
                    
                    this_fl <- paste(sep='',od,ds,'/',md,'/',ft,'/UMAP/',ds,'_',ft,'_',md,'_',hm_name,'UMAP.rds')
                    if(!file.exists(this_fl)){ if(verbose){cat(paste('SKIPPING:',ds,md,ft,hm_name,'\n'))}; next}
                    
                    ll <- readRDS(this_fl)
                    ll$cell <- rownames(ll)
                    ll$dataset <- ds
                    ll$feature <- ft
                    ll$method <- md
                    ll$Harmony <- hm
                    
                    if(any(cr %in% meta_col)){
                        toGet <- cr[which(cr %in% meta_col)]
                        toConv <- c('RNA_cLISI'='lisi_phenotype','RNA_iLISI'='lisi_sample')
                        this_fl <- paste(sep='',id,ds,'/',ds,'_meta.rds')
                        
                        rr <- getCol_fromFile(ll,this_fl,toGet,toConv,fl2_str='meta')
                        ll <- cbind(ll,rr)
                    }
                    
                    if(any(cr %in% lisi_col)){
                        toGet <- cr[which(cr %in% lisi_col)]
                        toConv <- c('ATAC_cLISI'='lisi_phenotype','ATAC_iLISI'='lisi_sample')
                        this_fl <- paste(sep='',od,ds,'/',md,'/',ft,'/LISImetrics/',ds,'_',ft,'_',md,'_',hm_name,
                                         'LISImetrics.rds')
                        
                        rr <- getCol_fromFile(ll,this_fl,toGet,toConv,fl2_str='LISI')
                        ll <- cbind(ll,rr)
                    }
                    
                    if(any(cr %in% PH_col)){
                        toGet <- cr[which(cr %in% PH_col)]
                        toConv <- c('ATAC_PH200_mKNN'='mKNN','ATAC_PH200_lsKLD'='lsKLD')
                        this_fl <- paste(sep='',od,ds,'/',md,'/',ft,'/NNmetrics/',ds,'_',ft,'_',md,'_',hm_name,
                                         'RNAPCAHarmony200_NNmetrics.rds')
                        
                        rr <- getCol_fromFile(ll,this_fl,toGet,toConv,fl2_str='PH NN')
                        ll <- cbind(ll,rr)
                    }
                    
                    if(any(cr %in% ST_col)){
                        toGet <- cr[which(cr %in% ST_col)]
                        toConv <- c('ATAC_ST200_mKNN'='mKNN','ATAC_ST200_lsKLD'='lsKLD')
                        this_fl <- paste(sep='',od,ds,'/',md,'/',ft,'/NNmetrics/',ds,'_',ft,'_',md,'_',hm_name,
                                         'RNASeurat200_NNmetrics.rds')
                        
                        rr <- getCol_fromFile(ll,this_fl,toGet,toConv,fl2_str='ST NN')
                        ll <- cbind(ll,rr)
                    }
                    
                    
                    rownames(ll) <- NULL
                    if(first){
                        toPlot <- ll
                        first <- FALSE
                    } else {
                        toPlot <- rbind(toPlot,ll)
                    }
                }
            }
        }
    }
    
    return(toPlot)
}



## job requirements

dicToVec <- function(dic,thisType,thisCol){
    if(!any(thisType %in% unique(dic$type))) stop('That type not in dictionary.')
    if(!(thisCol %in% colnames(dic))) stop('That column not in dictionary.')
    
    ll <- dic[which(dic$type==thisType),]
    vec <- ll[,thisCol]
    names(vec) <- ll$variable
    
    return(vec)
}

conv_hms_to_min <- function(ll,addSys=TRUE){
    if(grepl('[0-9]+:[0-9]+\\.[0-9]+',ll)){
        #nanoseconds - 0:46.93
        split <- str_split_fixed(ll,':',2)
        timeMin <- as.numeric(split[,1])+as.numeric(split[,2])/60
    } else if(grepl('[0-9]+:[0-9]+:[0-9]+',ll)){
        #hms - 2:23:35
        split <- str_split_fixed(ll,':',3)
        timeMin <- as.numeric(split[,1])*60+as.numeric(split[,2])+as.numeric(split[,3])/60
    } else if(grepl('[0-9\\.]+user [0-9\\.]+system',ll)){
        split <- str_split_fixed(ll,' ',2)
        timeMin <- as.numeric(sub('user','',split[,1]))
        if(addSys) {
            sy <- as.numeric(sub('system','',split[,2]))
            timeMin <- timeMin+sy
        }
        timeMin <- timeMin/60
    } else if(grepl('[0-9\\.]+user',ll)){
        timeMin <- as.numeric(sub('user','',ll))/60
    } else {
        print(paste('Unrecognized time format:',ll))
        timeMin <- ll
    }
    return(timeMin)
}

split_vec_function <- function(ll,max_col,typ){
    if(!(typ %in% c('mem','time'))) stop('type must be either time or memory')
    
    rr <- str_split_fixed(ll,'_',max_col)
    ds <- rr[1]
    ft <- rr[2]
    md <- rr[3]
    
    tm <- NA; st <- NA; gp <- NA; hm <- FALSE;
    for(idx in max_col:4){
        if(rr[idx]=="") {next}
        if(is.na(tm)) {tm <- rr[idx]; next}
        if(is.na(st)) {st <- rr[idx]; next}
        if(st=='NNmetrics' & is.na(gp)) {gp <- rr[idx]}
        if(idx==4 & rr[idx]=='Harmony') {hm <- TRUE}
    }
    
    #note that cmds-[0-9]+ steps contain both Harmony and not Harmony, so it should be NA for Harmony
    if(grepl('cmds-[0-9]+',st)) {hm <- NA}
    
    if(typ=='time'){
        tm <- conv_hms_to_min(tm)
    } else if(typ=='mem'){
        tm <- as.numeric(tm)/1000000 #KB -> GB
    } else {
        stop('Should not have gotten here')
    }
    
    return(data.frame('dataset'=ds,'feature'=ft,'method'=md,'Harmony'=hm,
                      'step'=st,'substep'=gp,'val'=tm,stringsAsFactors=FALSE))
}

jobReq_processing <- function(fl,typ,verbose=FALSE){
    if(!(typ %in% c('mem','time'))) stop('type must be either time or memory')
    
    if(verbose) {print("load")}
    vec <- read.table(fl,header=FALSE,sep='\t',stringsAsFactors=FALSE)
    vec <- vec[,1]
    if(typ=='time') {
        vec <- sub('\\.out:','_',str_split_fixed(vec,'/',4)[,4])
        lab <- 'User+system (CPU-min)'
        forTitle <- 'Time'
    } else if(typ=='mem') {
        vec <- sub('\\.out:','_',str_split_fixed(sub('maxresident','',vec),'/',4)[,4])
        lab <- 'Max Resident Memory (GB)'
        forTitle <- 'Memory'
    } else {
        stop('should not have gotten here')
    }

    if(verbose) {print("process")}
    max_col <- max(unlist(lapply(list(vec),function(x){str_count(x,"_")})))+1
    this_df <- data.frame('dataset'=character(),'feature'=character(),'method'=character(),'Harmony'=logical(),
                          'step'=character(),'substep'=character(),'val'=numeric(),stringsAsFactors=FALSE)
    for(xx in vec){
        this_df <- rbind(this_df,split_vec_function(xx,max_col,typ))
    }
    if(verbose) {print(head(this_df))}

    if(verbose) {print("remove previous dnf runs")}
    ll <- this_df[which(grepl('cmds-[0-9]+',this_df$step)),] %>%
            summarise(keep_step = max(step),.by= c(dataset, feature, method, Harmony))
    keep_steps <- ll$keep_step
    rm_steps <- grep('cmds-[0-9]+',this_df$step, value=TRUE)
    rm_steps <- rm_steps[which(!(rm_steps %in% keep_steps))]
    this_df <- this_df[which(!(this_df$step %in% rm_steps)),]
    
    return(this_df)
}







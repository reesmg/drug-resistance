# functions for comparing small-molecule response to baseline gene expression

# load packages
library(ppcor)
library(reshape2)
library(glmnet)
library(magrittr)
library(cdsrmodels)
library(dplyr)

# random forest function
# takes as input:
# X (feature matrix; rows are named features, columns are named cell lines)
# Y (viability matrix; rows are named compounds, columns are named cell lines to match X column names)
# named compound (matching a row name in Y)
rf_feats <- function(X, Y, compound) {
  # select compound and transform X to right shape
  y <- Y[compound, ]
  X <-  t(X)
  
  # fit random forest
  res <- cdsrmodels::random_forest(X, y)
  
  # extract results and return
  rf_table <- res$model_table %>% dplyr::mutate(compound = compound)
  rf_table$feature <- gsub("\\.\\.(.*)\\.", " (\\1)", rf_table$feature)
  rf_table <- dplyr::rename(rf_table, gene=feature, RF.rank=rank)
  
  return(data.frame(rf_table, check.names=F))
}

# elastic net function
# takes as input:
# X (feature matrix; rows are named features, columns are named cell lines)
# Y (viability matrix; rows are named compounds, columns are named cell lines to match X column names)
# named compound (matching a row name in Y)
# an alpha (for lasso/ridge tradeoff)
lasso_feats <- function(X, Y, compound, alpha=0.5) {
  # pre-process data into correct shape (remove missing values)
  y <- Y[compound, ]; y <- y[is.finite(y)]
  X <- t(X)
  overlap <- dplyr::intersect(rownames(X), names(y))
  X <- X[overlap,]; y <- y[overlap]
  
  # fit glmnet
  fit <- glmnet::cv.glmnet(X, y, alpha = alpha)
  pearson = cor(predict(fit, X), y)
  
  # extract coefficients
  a <- coef(fit)
  fit_tab <- tibble(feature = a@Dimnames[[1]][a@i+1],
                    enet.coef = a@x) %>%
    dplyr::filter(feature != "(Intercept)") %>%
    dplyr::mutate(PearsonScore = pearson,
                  rank = rank(desc(abs(enet.coef))))
  
  # fit_tab$feature <- gsub(" .*", "", fit_tab$feature)
  fit_tab$compound <- compound
  fit_tab <- dplyr::rename(fit_tab, gene=feature, enet.rank=rank)
  
  return(data.frame(fit_tab, check.names=F))
}

# standard Pearson correlation function using Fisher's z transformation
# takes as input:
# X (continuous feature matrix; rows are named features, columns are named cell lines)
# Y (viability matrix; rows are named compounds, columns are named cell lines to match X column names)
# named compound (matching a row name in Y)
fishers_z_correlation <- function(X, Y, compound){
  y.match <- Y[compound, intersect(names(Y[compound, !is.na(Y[compound, ])]), colnames(X))]
  x.match <- X[,names(y.match)]
  bonferroni_cutoff <- -qnorm(0.025/nrow(X))
  
  #calculate correlations and scaled correlations using Fisher's z transformation
  cor.xy <- cor(y.match,t(x.match))[1,]
  zcor.xy <- (0.5 * log((1+cor.xy)/(1-cor.xy)) * sqrt(length(y.match) - 3))
  # cors <- data.frame(gene=gsub(" .*", "", names(cor.xy)), r=cor.xy, z_r=zcor.xy, pass_bonferroni=F, check.names=F)
  cors <- data.frame(gene=names(cor.xy), compound=compound, r=cor.xy, z_r=zcor.xy, pass_bonferroni=F, check.names=F)
  cors[which(abs(cors$z_r) > bonferroni_cutoff), "pass_bonferroni"] <- T
  return(cors)
}

# function to filter gene expression to genes whose expression varies greater than 'range_cut' units
# takes as input:
# X (expression matrix; rows are named genes, columns are named cell lines)
expression_range_filter <- function(x, range_cut){
  range.expression <- apply(x, 1, range, na.rm=T)
  y <- x[which((range.expression[2,] - range.expression[1,]) > range_cut), ] 
  return(y)
}

# semipartial correlation function
#often, multiple genes will have correlation coefficients very similar in value, particularly if they are strongly co-expressed. the sd_cut argument allows for keeping genes with a correlation value within a certain number of standard deviations of the top gene 
# takes as input:
# X (feature matrix; rows are named features, columns are named cell lines)
# Y (viability matrix; rows are named compounds, columns are named cell lines to match X column names)
# named compound (matching a row name in Y)
# as the function scales the correlation distribution, scale_cut is the number of standard deviations from the center required to consider a correlation significant
# exp_range: number of units a transcript must vary by (across cell lines) to be included as a hit
# coexp_cut: minimum Pearson correlation of gene coexpression across cell lines to be considered as a single association
iterative_semipartial <- function(X, Y, compound, scale_cut=3, sd_cut=0.25, exp_range=3, coexp_cut=0.8){
  y.match <- Y[compound, intersect(names(Y[compound, !is.na(Y[compound, ])]), colnames(X))]
  x.match <- X[,names(y.match)]
  # threshold for bonferroni adjustment based on number of genes
  bonferroni_cutoff <- -qnorm(0.025/nrow(X))
  
  if(length(y.match) < 9) return("Not enough cell lines")
  
  # apply expression range cutoff
  if(!is.na(exp_range)){
    x.match <- expression_range_filter(x.match, exp_range)
  }
  
  # calculate correlations and scaled correlations using Fisher's z transformation
  cor.xy <- cor(y.match,t(x.match))[1,]
  rank.xy <- rank(-1*abs(cor.xy))
  zcor.xy <- (0.5 * log((1+cor.xy)/(1-cor.xy)) * sqrt(length(y.match) - 3))
  # exclude the most and least sensitive cell lines to ensure correlations not driven by single outlier
  cor.xy.loo <- cor(y.match[-c(which.max(y.match), which.min(y.match))],
                    t(x.match[, -c(which.max(y.match), which.min(y.match))]))[1, ]
  zcor.xy.loo <- (0.5 * log((1+cor.xy.loo)/(1-cor.xy.loo)) * sqrt(length(y.match) - 2 - 3))
  
  # create scaled distribution of correlation coefficients for compound x
  zcor.xy.scale <- scale(zcor.xy)[,1]
  
  # catalog all significantly correlated genes:
  # z-transformed correlation exceeds bonferroni-corrected cutoff for number of genes and is an 'outlier' from the rest of the correlation distribution
  significant.positive.genes <- which(zcor.xy > bonferroni_cutoff & zcor.xy.loo > bonferroni_cutoff & zcor.xy.scale > scale_cut)
  significant.negative.genes <- which(zcor.xy < -1*bonferroni_cutoff & zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy.scale < -scale_cut)
  if(length(significant.positive.genes)==0 & length(significant.negative.genes)==0) return("No significantly correlated genes")
  
  # create output data.frame; flag only the most significant gene (and any genes close to it in terms of significance using the sd_cut variable cutoff)
  hit_genes <- data.frame(compound=compound, gene=c(names(significant.positive.genes), names(significant.negative.genes)), r=cor.xy[c(significant.positive.genes, significant.negative.genes)], unadj_rank=rank.xy[c(significant.positive.genes, significant.negative.genes)], z_r=zcor.xy[c(significant.positive.genes, significant.negative.genes)], z_scale=zcor.xy.scale[c(significant.positive.genes, significant.negative.genes)], pass_sd_cut=F)
  hit_genes[which(hit_genes$gene %in% c(
    names(which(zcor.xy.scale[significant.positive.genes] > (max(abs(zcor.xy.scale), na.rm = T) - sd_cut))),
    names(which(zcor.xy.scale[significant.negative.genes] < (-1*max(abs(zcor.xy.scale), na.rm = T) + sd_cut)))
  )), "pass_sd_cut"] <- T 
  
  # calculate semi-partial correlations adjusting for the top-ranked gene
  # make data frame of viability values and corresponding expression values for top-ranked gene
  semi.partial.xy.max <- data.frame(x=y.match, z1=x.match[names(which.max(abs(zcor.xy.scale))), ], stringsAsFactors = F)
  # add expression values for all remaining genes
  semi.partial.xy.max <- cbind(semi.partial.xy.max, t(x.match[-which(row.names(x.match)==names(which.max(abs(zcor.xy.scale)))), ]))
  # generate semi-partial correlation values between viability and expression of every gene (adjusted for max gene)
  sp.results <- apply(semi.partial.xy.max[, 3:ncol(semi.partial.xy.max)], 2, function(y) spcor.test(semi.partial.xy.max[, "x"], y, semi.partial.xy.max[, "z1"]))
  semi.partial.xy.max <- unlist(lapply(sp.results, `[[`, 1))
  semi.partial.xy.max.pval <- unlist(lapply(sp.results, `[[`, 2))
  semipartial_results <- data.frame(gene=names(semi.partial.xy.max), r1=semi.partial.xy.max, p1=semi.partial.xy.max.pval, rank1=rank(-1*abs(semi.partial.xy.max)))
  rm(sp.results)
  
  # repeat semi-partial correlations, now adjusting for the top-ranked gene and top semi-partial correlated gene
  semi.partial.xy.second <- data.frame(x=y.match, z1=x.match[names(which.max(abs(zcor.xy.scale))), ], z2=x.match[names(which.max(abs(semi.partial.xy.max))), ], stringsAsFactors = F)
  semi.partial.xy.second <- cbind(semi.partial.xy.second, t(x.match[-c(which(row.names(x.match)==names(which.max(abs(zcor.xy.scale)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.max))))), ]))
  sp.results <- apply(semi.partial.xy.second[, (max(grep("z", colnames(semi.partial.xy.second)))+1):ncol(semi.partial.xy.second)], 2, function(y) spcor.test(semi.partial.xy.second[, "x"], y, semi.partial.xy.second[, grep("z", colnames(semi.partial.xy.second))]))
  semi.partial.xy.second <- unlist(lapply(sp.results, `[[`, 1))
  semi.partial.xy.second.pval <- unlist(lapply(sp.results, `[[`, 2))
  semipartial_results <- dplyr::full_join(semipartial_results, 
                                          data.frame(gene=names(semi.partial.xy.second), r2=semi.partial.xy.second, p2=semi.partial.xy.second.pval, rank2=rank(-1*abs(semi.partial.xy.second))))
  rm(sp.results)
  
  # calculate semi-partial correlations adjusting for the top-ranked gene and top 2 semi-partial genes
  semi.partial.xy.third <- data.frame(x=y.match, z1=x.match[names(which.max(abs(zcor.xy.scale))), ], z2=x.match[names(which.max(abs(semi.partial.xy.max))), ], z3=x.match[names(which.max(abs(semi.partial.xy.second))), ], stringsAsFactors = F)
  semi.partial.xy.third <- cbind(semi.partial.xy.third, t(x.match[-c(which(row.names(x.match)==names(which.max(abs(zcor.xy.scale)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.max)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.second))))), ]))
  sp.results <- apply(semi.partial.xy.third[, (max(grep("z", colnames(semi.partial.xy.third)))+1):ncol(semi.partial.xy.third)], 2, function(y) spcor.test(semi.partial.xy.third[, "x"], y, semi.partial.xy.third[, grep("z", colnames(semi.partial.xy.third))]))
  semi.partial.xy.third <- unlist(lapply(sp.results, `[[`, 1))
  semi.partial.xy.third.pval <- unlist(lapply(sp.results, `[[`, 2))
  semipartial_results <- dplyr::full_join(semipartial_results, 
                                          data.frame(gene=names(semi.partial.xy.third), r3=semi.partial.xy.third, p3=semi.partial.xy.third.pval, rank3=rank(-1*abs(semi.partial.xy.third))))
  
  rm(sp.results)
  
  # calculate semi-partial correlations adjusting for the top-ranked gene and top 3 semi-partial genes
  # add error trapping to continue if no significant genes remain
  semi.partial.xy.fourth <- data.frame(x=y.match, z1=x.match[names(which.max(abs(zcor.xy.scale))), ], z2=x.match[names(which.max(abs(semi.partial.xy.max))), ], z3=x.match[names(which.max(abs(semi.partial.xy.second))), ], z4=x.match[names(which.max(abs(semi.partial.xy.third))), ], stringsAsFactors = F)
  semi.partial.xy.fourth <- cbind(semi.partial.xy.fourth, t(x.match[-c(which(row.names(x.match)==names(which.max(abs(zcor.xy.scale)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.max)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.second)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.third))))), ]))
  sp.results <- try(apply(semi.partial.xy.fourth[, (max(grep("z", colnames(semi.partial.xy.fourth)))+1):ncol(semi.partial.xy.fourth)], 2, function(y) spcor.test(semi.partial.xy.fourth[, "x"], y, semi.partial.xy.fourth[, grep("z", colnames(semi.partial.xy.fourth))])), silent = T)
  if(class(sp.results)=="try-error"){
    # create dummy vectors
    semi.partial.xy.fourth <- rep(0, (length(semi.partial.xy.third)-1))
    semi.partial.xy.fourth.pval <- rep(1, (length(semi.partial.xy.third.pval)-1))
    names(semi.partial.xy.fourth) <- names(semi.partial.xy.fourth.pval) <- names(semi.partial.xy.third)[-which.max(abs(semi.partial.xy.third))]
    semipartial_results <- dplyr::full_join(semipartial_results, 
                                            data.frame(gene=names(semi.partial.xy.fourth), r4=NA, p4=NA, rank4=NA))
  }else {
    semi.partial.xy.fourth <- unlist(lapply(sp.results, `[[`, 1))
    semi.partial.xy.fourth.pval <- unlist(lapply(sp.results, `[[`, 2))
    semipartial_results <- dplyr::full_join(semipartial_results, 
                                            data.frame(gene=names(semi.partial.xy.fourth), r4=semi.partial.xy.fourth, p4=semi.partial.xy.fourth.pval, rank4=rank(-1*abs(semi.partial.xy.fourth))))
  }
  rm(sp.results)
  
  # calculate semi-partial correlations adjusting for the top-ranked gene and top 4 semi-partial genes
  # add error trapping to continue if no significant genes remain
  semi.partial.xy.fifth <- data.frame(x=y.match, z1=x.match[names(which.max(abs(zcor.xy.scale))), ], z2=x.match[names(which.max(abs(semi.partial.xy.max))), ], z3=x.match[names(which.max(abs(semi.partial.xy.second))), ], z4=x.match[names(which.max(abs(semi.partial.xy.third))), ], z5=x.match[names(which.max(abs(semi.partial.xy.fourth))), ], stringsAsFactors = F)
  semi.partial.xy.fifth <- cbind(semi.partial.xy.fifth, t(x.match[-c(which(row.names(x.match)==names(which.max(abs(zcor.xy.scale)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.max)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.second)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.third)))), which(row.names(x.match)==names(which.max(abs(semi.partial.xy.fourth))))), ]))
  sp.results <- try(apply(semi.partial.xy.fifth[, (max(grep("z", colnames(semi.partial.xy.fifth)))+1):ncol(semi.partial.xy.fifth)], 2, function(y) spcor.test(semi.partial.xy.fifth[, "x"], y, semi.partial.xy.fifth[, grep("z", colnames(semi.partial.xy.fifth))])), silent = T)
  if(class(sp.results)=="try-error"){
    semi.partial.xy.fifth <- rep(0, (length(semi.partial.xy.fourth)-1))
    semi.partial.xy.fifth.pval <- rep(1, (length(semi.partial.xy.fourth.pval)-1))
    names(semi.partial.xy.fifth) <- names(semi.partial.xy.fifth.pval) <- names(semi.partial.xy.fourth)[-which.max(abs(semi.partial.xy.fourth))]
    semipartial_results <- dplyr::full_join(semipartial_results, 
                                            data.frame(gene=names(semi.partial.xy.fifth), r5=NA, p5=NA, rank5=NA))
  }else {
    semi.partial.xy.fifth <- unlist(lapply(sp.results, `[[`, 1))
    semi.partial.xy.fifth.pval <- unlist(lapply(sp.results, `[[`, 2))
    semipartial_results <- dplyr::full_join(semipartial_results, 
                                            data.frame(gene=names(semi.partial.xy.fifth), r5=semi.partial.xy.fifth, p5=semi.partial.xy.fifth.pval, rank5=rank(-1*abs(semi.partial.xy.fifth))))
  }
  rm(sp.results)
  
  # add semipartial correlations, ranks, pvals to genes already in (unadjusted) hit_genes table
  hit_genes <- dplyr::left_join(hit_genes, semipartial_results)
  
  # store significant genes from any step of the semipartial. Kept genes must be significantly correlated on their own (w/o semipartial) and pass 'outlier' test and significance (adjusted p-value) (and not already be in the top correlated list)
  if(length(names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) > scale_cut)]) > 0){
    sp_hit_genes <- try(
      data.frame(gene=
                   c(
                     # must be significantly correlated on its own:               
                     intersect(names(which(zcor.xy.loo > bonferroni_cutoff & zcor.xy > bonferroni_cutoff)),
                               c(
                                 names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) > scale_cut)],
                                 names(semi.partial.xy.second)[which(scale(semi.partial.xy.second) > scale_cut)],
                                 names(semi.partial.xy.third)[which(scale(semi.partial.xy.third) > scale_cut)],
                                 names(semi.partial.xy.fourth)[which(scale(semi.partial.xy.fourth) > scale_cut)],
                                 names(semi.partial.xy.fifth)[which(scale(semi.partial.xy.fifth) > scale_cut)]
                               )
                     ),
                     # same, but for negative genes:
                     intersect(names(which(zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy < -1*bonferroni_cutoff)),
                               c(
                                 names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) < -scale_cut)],
                                 names(semi.partial.xy.second)[which(scale(semi.partial.xy.second) < -scale_cut)],
                                 names(semi.partial.xy.third)[which(scale(semi.partial.xy.third) < -scale_cut)],
                                 names(semi.partial.xy.fourth)[which(scale(semi.partial.xy.fourth) < -scale_cut)],
                                 names(semi.partial.xy.fifth)[which(scale(semi.partial.xy.fifth) < -scale_cut)]
                               )
                     )
                   )
                 , pass_semipartial_sd_cut=F)
      , silent = T)
    
    #additional step to mark as pass_sd_cut=T those genes that are the most significant (pass sd_cut)
    if(class(sp_hit_genes)!="try-error"){
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           c(
                             intersect(names(which(zcor.xy.loo > bonferroni_cutoff & zcor.xy > bonferroni_cutoff)),
                                       c(
                                         names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) > (max(max(abs(scale(semi.partial.xy.max)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))],
                                         names(semi.partial.xy.second)[which(scale(semi.partial.xy.second) > (max(max(abs(scale(semi.partial.xy.second)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))],
                                         names(semi.partial.xy.third)[which(scale(semi.partial.xy.third) > (max(max(abs(scale(semi.partial.xy.third)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))],
                                         names(semi.partial.xy.fourth)[which(scale(semi.partial.xy.fourth) > (max(max(abs(scale(semi.partial.xy.fourth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))],
                                         names(semi.partial.xy.fifth)[which(scale(semi.partial.xy.fifth) > (max(max(abs(scale(semi.partial.xy.fifth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))]
                                       )
                             ),
                             intersect(names(which(zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy < -1*bonferroni_cutoff)),
                                       c(
                                         names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) < -1*((max(max(abs(scale(semi.partial.xy.max)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))],
                                         names(semi.partial.xy.second)[which(scale(semi.partial.xy.second) < -1*((max(max(abs(scale(semi.partial.xy.second)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))],
                                         names(semi.partial.xy.third)[which(scale(semi.partial.xy.third) < -1*((max(max(abs(scale(semi.partial.xy.third)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))],
                                         names(semi.partial.xy.fourth)[which(scale(semi.partial.xy.fourth) < -1*((max(max(abs(scale(semi.partial.xy.fourth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))],
                                         names(semi.partial.xy.fifth)[which(scale(semi.partial.xy.fifth) < -1*((max(max(abs(scale(semi.partial.xy.fifth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))]
                                       )
                             )
                           )
      ), "pass_semipartial_sd_cut"] <- T
      
      # pass pval in any iteration
      sp_hit_genes$pass_semipartial_pval <- F
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           c(
                             intersect(names(which(zcor.xy.loo > bonferroni_cutoff & zcor.xy > bonferroni_cutoff)),
                                       c(
                                         names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) > (max(max(abs(scale(semi.partial.xy.max)), na.rm = T),(scale_cut+sd_cut)) - sd_cut) & scale(semi.partial.xy.max) > scale_cut & semi.partial.xy.max.pval < 0.025/nrow(X))],
                                         names(semi.partial.xy.second)[which(scale(semi.partial.xy.second) > (max(max(abs(scale(semi.partial.xy.second)), na.rm = T),(scale_cut+sd_cut)) - sd_cut) & scale(semi.partial.xy.second) > scale_cut & semi.partial.xy.second.pval < 0.025/nrow(X))],
                                         names(semi.partial.xy.third)[which(scale(semi.partial.xy.third) > (max(max(abs(scale(semi.partial.xy.third)), na.rm = T),(scale_cut+sd_cut)) - sd_cut) & scale(semi.partial.xy.third) > scale_cut & semi.partial.xy.third.pval < 0.025/nrow(X))],
                                         names(semi.partial.xy.fourth)[which(scale(semi.partial.xy.fourth) > (max(max(abs(scale(semi.partial.xy.fourth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut) & scale(semi.partial.xy.fourth) > scale_cut & semi.partial.xy.fourth.pval < 0.025/nrow(X))],
                                         names(semi.partial.xy.fifth)[which(scale(semi.partial.xy.fifth) > (max(max(abs(scale(semi.partial.xy.fifth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut) & scale(semi.partial.xy.fifth) > scale_cut & semi.partial.xy.fifth.pval < 0.025/nrow(X))]
                                       )
                             ),
                             intersect(names(which(zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy < -1*bonferroni_cutoff)),
                                       c(
                                         names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) < -1*((max(max(abs(scale(semi.partial.xy.max)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)) & scale(semi.partial.xy.max) < -scale_cut & semi.partial.xy.max.pval < 0.025/nrow(X))],
                                         names(semi.partial.xy.second)[which(scale(semi.partial.xy.second) < -1*((max(max(abs(scale(semi.partial.xy.second)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)) & scale(semi.partial.xy.second) < -scale_cut & semi.partial.xy.second.pval < 0.025/nrow(X))],
                                         names(semi.partial.xy.third)[which(scale(semi.partial.xy.third) < -1*((max(max(abs(scale(semi.partial.xy.third)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)) & scale(semi.partial.xy.third) < -scale_cut & semi.partial.xy.third.pval < 0.025/nrow(X))],
                                         names(semi.partial.xy.fourth)[which(scale(semi.partial.xy.fourth) < -1*((max(max(abs(scale(semi.partial.xy.fourth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)) & scale(semi.partial.xy.fourth) < -scale_cut & semi.partial.xy.fourth.pval < 0.025/nrow(X))],
                                         names(semi.partial.xy.fifth)[which(scale(semi.partial.xy.fifth) < -1*((max(max(abs(scale(semi.partial.xy.fifth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)) & scale(semi.partial.xy.fifth) < -scale_cut & semi.partial.xy.fifth.pval < 0.025/nrow(X))]
                                       )
                             )
                           )
      ), "pass_semipartial_pval"] <- T
      
      # which iteration of the semipartial did the result come from?
      sp_hit_genes$iteration <- "none"
      # positive
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo > bonferroni_cutoff & zcor.xy > bonferroni_cutoff)),
                                     names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) > (max(max(abs(scale(semi.partial.xy.max)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))])
      ), "iteration"] <- "1"
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo > bonferroni_cutoff & zcor.xy > bonferroni_cutoff)),
                                     names(semi.partial.xy.second)[which(scale(semi.partial.xy.second) > (max(max(abs(scale(semi.partial.xy.second)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))])
      ), "iteration"] <- "2"
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo > bonferroni_cutoff & zcor.xy > bonferroni_cutoff)),
                                     names(semi.partial.xy.third)[which(scale(semi.partial.xy.third) > (max(max(abs(scale(semi.partial.xy.third)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))])
      ), "iteration"] <- "3"
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo > bonferroni_cutoff & zcor.xy > bonferroni_cutoff)),
                                     names(semi.partial.xy.fourth)[which(scale(semi.partial.xy.fourth) > (max(max(abs(scale(semi.partial.xy.fourth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))])
      ), "iteration"] <- "4"
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo > bonferroni_cutoff & zcor.xy > bonferroni_cutoff)),
                                     names(semi.partial.xy.fifth)[which(scale(semi.partial.xy.fifth) > (max(max(abs(scale(semi.partial.xy.fifth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut))])
      ), "iteration"] <- "5"
      
      # negative
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy < -1*bonferroni_cutoff)),
                                     names(semi.partial.xy.max)[which(scale(semi.partial.xy.max) < -1*((max(max(abs(scale(semi.partial.xy.max)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))])
      ), "iteration"] <- "1neg"
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy < -1*bonferroni_cutoff)),
                                     names(semi.partial.xy.second)[which(scale(semi.partial.xy.second) < -1*((max(max(abs(scale(semi.partial.xy.second)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))])
      ), "iteration"] <- "2neg"
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy < -1*bonferroni_cutoff)),
                                     names(semi.partial.xy.third)[which(scale(semi.partial.xy.third) < -1*((max(max(abs(scale(semi.partial.xy.third)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))])
      ), "iteration"] <- "3neg"
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy < -1*bonferroni_cutoff)),
                                     names(semi.partial.xy.fourth)[which(scale(semi.partial.xy.fourth) < -1*((max(max(abs(scale(semi.partial.xy.fourth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))])
      ), "iteration"] <- "4neg"
      sp_hit_genes[which(sp_hit_genes$gene %in% 
                           intersect(names(which(zcor.xy.loo < -1*bonferroni_cutoff & zcor.xy < -1*bonferroni_cutoff)),
                                     names(semi.partial.xy.fifth)[which(scale(semi.partial.xy.fifth) < -1*((max(max(abs(scale(semi.partial.xy.fifth)), na.rm = T),(scale_cut+sd_cut)) - sd_cut)))])
      ), "iteration"] <- "5neg"
      
      # add genes that pass semipartial sd cut to original df for genes already in the table
      hit_genes <- dplyr::left_join(hit_genes, sp_hit_genes)
      
      # add genes not in the prior table
      if(!all(sp_hit_genes$gene %in% hit_genes$gene)){
        # restrict semipartial hit table to 'new' genes and add underlying semipartial data
        sp_hit_genes <- sp_hit_genes[which(sp_hit_genes$gene %in% setdiff(sp_hit_genes$gene, hit_genes$gene)),] %>%
          dplyr::left_join(semipartial_results)
        # add additional necessary columns
        sp_hit_genes$compound <- compound
        sp_hit_genes$r <-cor.xy[sp_hit_genes$gene]
        sp_hit_genes$z_r <- zcor.xy[sp_hit_genes$gene]
        sp_hit_genes$z_scale <- zcor.xy.scale[sp_hit_genes$gene]
        sp_hit_genes$unadj_rank <- rank.xy[sp_hit_genes$gene]
        
        #add semipartial hit genes to original hits
        hit_genes <- hit_genes %>%
          dplyr::full_join(
            sp_hit_genes
          )
      }
    }
  }
  
  # annotate which genes are highly co-expressed with each other across the matched cell lines using the coexp_cut cutoff
  if(nrow(hit_genes)>1){
    # calculate expression correlations between all genes in the hit_genes df
    total.genes.cor <- cor(t(x.match[hit_genes$gene, ]))
    total.genes.cor <- reshape2::melt(total.genes.cor)
    # remove self-correlation gene pairs and those that fail the co-expression cutoff
    total.genes.cor <- total.genes.cor[which(total.genes.cor[, 1]!=total.genes.cor[, 2] & abs(total.genes.cor[, 3]) > coexp_cut), 1:2]
    #collapse multiple rows of the same gene
    if(nrow(total.genes.cor)>0){
      total.genes.cor <- aggregate(total.genes.cor[,2], by=list(total.genes.cor[,1]), paste, collapse=";") %>%
        dplyr::rename(coexpressed_genes=x)
      hit_genes <- dplyr::left_join(hit_genes, total.genes.cor, by=c("gene"="Group.1"))
    }
  }

  return(hit_genes)
}






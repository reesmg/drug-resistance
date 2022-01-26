# example analysis comparing compound viability to baseline gene expression 
# using GSK-J4 AUC from Lochmann et al. STM DOI: 10.1126/scitranslmed.aao4680
# AUC data from Supplementary Table 1

# import necessary libraries and functions
suppressMessages(source("./src/analysis_functions.R"))

# load viability matrix (rows are named compounds, columns are named cell lines)
# update the path as necessary for the provided file
viability <- readRDS("./src/Lochmann_et_al_GSK-J4_AUC.rds")
viability <- matrix(viability$AUC, nrow=1, dimnames=list("GSK-J4", viability$DepMap_ID))
# load expression matrix from the dependency map 
expression <- data.table::fread("https://ndownloader.figshare.com/files/22897976") 
# reformat: rows are named genes, columns are named cell lines
DepMap_IDs <- expression$V1
expression <- t(data.matrix(expression[, 2:ncol(expression)]))
colnames(expression) <- DepMap_IDs
rm(DepMap_IDs)

# simple correlation analysis
cor_example <- fishers_z_correlation(Y=viability, X=expression, compound="GSK-J4")

# lasso, random forest, and semipartial analyses
lasso_example <- lasso_feats(Y=viability, X=expression, compound="GSK-J4")
rf_example <- rf_feats(Y=viability, X=expression, compound="GSK-J4")
semipartial_example <- iterative_semipartial(Y=viability, compound="GSK-J4", X=expression, scale_cut=3, sd_cut=0.25, exp_range=3, coexp_cut=0.8)

# plot semipartial results
# melt correlation data
p1 <- reshape2::melt(semipartial_example[, c("gene", "r", "r1", "r2", "r3", "r4", "r5")])
# order by rank of standard (unadjusted) correlation and calculate ranks for each 'step' of semipartial
p1$gene <- factor(p1$gene, levels = semipartial_example[order(semipartial_example$unadj_rank), "gene"])
p1 <- p1 %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(by_group_rank=rank(-1*abs(value))) %>%
  dplyr::ungroup()
ggplot(p1,  aes(y = value, x = variable, color=gene)) +
  geom_point() + 
  theme_bw() + #remove gray background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove gridlines
  theme(panel.border = element_blank()) + #remove border lines
  theme(axis.line = element_line(colour = "black")) + #restore axis lines
  xlab("semipartial step") +
  scale_x_discrete(labels=as.character(0:5)) +
  ylab("r") +
  gghighlight(gene %in% dplyr::filter(p1, by_group_rank==1)$gene, label_key = gene, use_direct_label = F)

#merge lasso and rf w/ correlations, create table for export
output_table <- cor_example %>%
  dplyr::full_join(rf_example[, c("gene", "compound", "RF.imp.mean", "RF.imp.sd", "RF.imp.stability", "RF.rank")]) %>%
  dplyr::full_join(lasso_example[, c("gene", "compound", "enet.coef", "enet.rank")]) %>%
  dplyr::full_join(semipartial_example) %>%
  dplyr::filter(gene %in% unique(c(lasso_example$gene, rf_example$gene, semipartial_example$gene)))

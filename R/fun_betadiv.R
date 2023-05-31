#' beta diversity decomposition
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

beta_div_decomp<- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  matrix.pa <- vegan::decostand(df_mean, "pa")
  rownames(matrix.pa) <- meta_mean$arms
  colnames(matrix.pa) <- meta_mean$arms
  
  B.pair.pa <- betapart::beta.pair(matrix.pa, index.family = "jaccard")
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- 1-B.pair.pa$beta.jac
  
  library(reshape2)
  
  #### inter/intra ####
  
  melt(mat.turn)
  df.turn <- melt(mat.turn, varnames = c("row", "col"))
  df.turn <- subset(df.turn, row != col) 
  
  df.turn$row <- substr(df.turn$row, 1, 5)
  df.turn$col <- substr(df.turn$col, 1, 5)
  df.turn$same_value <- ifelse(df.turn$row == df.turn$col, "Yes", "No")
  
  df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  df.nest <- subset(df.nest, row != col)
  
  df.nest$row <- substr(df.nest$row, 1, 5)
  df.nest$col <- substr(df.nest$col, 1, 5)
  df.nest$same_value <- ifelse(df.nest$row == df.nest$col, "Yes", "No")
  
  
  df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  df.jacc <- subset(df.jacc, row != col)
  
  df.jacc$row <- substr(df.jacc$row, 1, 5)
  df.jacc$col <- substr(df.jacc$col, 1, 5)
  df.jacc$same_value <- ifelse(df.jacc$row == df.jacc$col, "Yes", "No")
  
  
  decomp.pa <- as.data.frame(cbind(df.nest$same_value, df.turn$value, df.nest$value  ,df.jacc$value))
  
  colnames(decomp.pa) <- c("intrasite","Turnover", "Nestedness", "Jaccard_diss")
  
  decomp.pa$Jaccard_diss <- abs(as.numeric(decomp.pa$Jaccard_diss))
  decomp.pa$Nestedness <- abs(as.numeric(decomp.pa$Nestedness))
  decomp.pa$Turnover <- abs(as.numeric(decomp.pa$Turnover))
  library(ggplot2)
 
   
  res.aov <- rstatix::anova_test(decomp.pa, Turnover ~ intrasite)
  res.kru <- rstatix::kruskal_test(decomp.pa, Turnover ~ intrasite)
  p.sed <- rstatix::wilcox_test(decomp.pa, Turnover ~ intrasite, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  intrainter = c("between ARMS \n of the same batch", "between ARMS \n of different batch")
  g <- ggplot(decomp.pa, aes(x = fct_relevel(intrasite, "Yes", "No"), y = Turnover)) +
    geom_boxplot(fill =  c("lightblue","lightblue") ) +
    labs(title = "",
         x = "Comparisons",
         y = "Turnover component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intrainter) +
    theme_classic() +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  g

  path_to_turn_intrainter <- paste0("outputs/beta/turnover_intrainter.pdf")
  ggsave(filename = path_to_turn_intrainter, plot = v, width = 6, height = 6)
   
   #### Immersion time ####
   #turnover
   mat.turn
   df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
   df.turn <- subset(df.turn, row != col)
   
   df.turn$row <- substr(df.turn$row, 1, 5)
   df.turn$col <- substr(df.turn$col, 1, 5)
   df.turn$same_value <- ifelse(df.turn$row == df.turn$col, "Yes", "No")
   df.turn <- subset(df.turn, row != col)
  
   df.turn$comp_imm <- ifelse(df.turn$col == "CINA1" & df.turn$row == "CINA2", "six_one",
                                ifelse(df.turn$col == "CINA3" & df.turn$row == "CINA4", "six_one",
                                       ifelse(df.turn$col == "CINA2" & df.turn$row == "RUNA2", "one_two",
                                            ifelse(df.turn$col == "CINA1" & df.turn$row == "RUNA2", "six_two",
                                                   NA))))
   
   df.turn.comp.imm <- df.turn[!is.na(df.turn$comp_imm),]
   
   imm_turn <- df.turn.comp.imm$value
   
   
   #nestedness
   mat.nest
   df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
   df.nest <- subset(df.nest, row != col)
   
   df.nest$row <- substr(df.nest$row, 1, 5)
   df.nest$col <- substr(df.nest$col, 1, 5)
   df.nest$same_value <- ifelse(df.nest$row == df.nest$col, "Yes", "No")
   df.nest <- subset(df.nest, row != col)
   
   df.nest$comp_imm <- ifelse(df.nest$col == "CINA1" & df.nest$row == "CINA2", "six_one",
                              ifelse(df.nest$col == "CINA3" & df.nest$row == "CINA4", "six_one",
                                     ifelse(df.nest$col == "CINA2" & df.nest$row == "RUNA2", "one_two",
                                            ifelse(df.nest$col == "CINA1" & df.nest$row == "RUNA2", "six_two",
                                                   NA))))
   df.nest.comp.imm <- df.nest[!is.na(df.nest$comp_imm),]
   
   
   imm_nest <- df.nest.comp.imm$value
   
   #jac similarity
   mat.jacc
   df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
   df.jacc <- subset(df.jacc, row != col)
   
   df.jacc$row <- substr(df.jacc$row, 1, 5)
   df.jacc$col <- substr(df.jacc$col, 1, 5)
   df.jacc$same_value <- ifelse(df.jacc$row == df.jacc$col, "Yes", "No")
   df.jacc <- subset(df.jacc, row != col)
   
   df.jacc$comp_imm <- ifelse(df.jacc$col == "CINA1" & df.jacc$row == "CINA2", "six_one",
                              ifelse(df.jacc$col == "CINA3" & df.jacc$row == "CINA4", "six_one",
                                     ifelse(df.jacc$col == "CINA2" & df.jacc$row == "RUNA2", "one_two",
                                            ifelse(df.jacc$col == "CINA1" & df.jacc$row == "RUNA2", "six_two",
                                                   NA))))
   df.jacc.comp.imm <- df.jacc[!is.na(df.jacc$comp_imm),]
   
   imm_jacc <- df.jacc.comp.imm$value
   
   df_imm_time <- data.frame(turn = imm_turn,
                             nest = imm_nest,
                             jacc = imm_jacc,
                             comp_imm = df.turn.comp.imm$comp_imm)
   
   
   res.aov <- rstatix::anova_test(df_imm_time, turn ~ comp_imm)
   res.kru <- rstatix::kruskal_test(df_imm_time, turn ~ comp_imm)
   p.sed <- rstatix::wilcox_test(df_imm_time, turn ~ comp_imm, p.adjust.method = "bonferroni")
   p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
   
   library(ggpubr)
   library(forcats)
   comp = c("between six months \n and one year", "between one year \n and two years", "between six months \n and two years")
   
   v <- ggplot(df_imm_time, aes(x = fct_relevel(comp_imm, "six_one", "one_two", "six_two"), y = turn)) +
     geom_boxplot(fill =  c("coral","coral","coral") ) +
     labs(title = "",
          x = "Comparisons",
          y = "Turnover component") +
     theme(legend.position = "none") +
     scale_x_discrete(labels=comp) +
     theme_classic() +
     stat_pvalue_manual(p.sed)  +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   v
   
   res.aov <- rstatix::anova_test(df_imm_time, nest ~ comp_imm)
   res.kru <- rstatix::kruskal_test(df_imm_time, nest ~ comp_imm)
   p.sed <- rstatix::wilcox_test(df_imm_time, nest ~ comp_imm, p.adjust.method = "bonferroni")
   p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
   
   f <- ggplot(df_imm_time, aes(x = fct_relevel(comp_imm, "six_one", "one_two", "six_two"), y = nest)) +
     geom_boxplot(fill =  c("coral","coral","coral") ) +
     labs(title = "",
          x = "Comparisons",
          y = "Nestedness component") +
     theme(legend.position = "none") +
     scale_x_discrete(labels=comp) +
     theme_classic() +
     stat_pvalue_manual(p.sed) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   f
   
   res.aov <- rstatix::anova_test(df_imm_time, jacc ~ comp_imm)
   res.kru <- rstatix::kruskal_test(df_imm_time, jacc ~ comp_imm)
   p.sed <- rstatix::wilcox_test(df_imm_time, jacc ~ comp_imm, p.adjust.method = "bonferroni")
   p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
   
   m <- ggplot(df_imm_time, aes(x = fct_relevel(comp_imm, "six_one", "one_two", "six_two"), y = jacc)) +
     geom_boxplot(fill =  c("coral","coral","coral") ) +
     labs(title = "",
          x = "Comparisons",
          y = "Jaccard similarity") +
     theme(legend.position = "none") +
     scale_x_discrete(labels=comp) +
     theme_classic() +
     stat_pvalue_manual(p.sed) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   m
   
   fin <- cowplot::plot_grid(v,f,m,
                             ncol = 3,
                             nrow = 1)
   
   
   path_to_boxplot <- paste0("outputs/beta/boxplot_imm_tim.pdf")
   ggsave(filename =  path_to_boxplot, plot = fin, width = 25, height = 16)

   #### deployment season ####
   #turnover
   mat.turn
   df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
   df.turn <- subset(df.turn, row != col)
   
   df.turn$row <- substr(df.turn$row, 1, 5)
   df.turn$col <- substr(df.turn$col, 1, 5)
   df.turn$same_value <- ifelse(df.turn$row == df.turn$col, "Yes", "No")
   # df.turn <- subset(df.turn, row != col)
   
   df.turn$comp_deploy <- ifelse(df.turn$col == "CINA1" & df.turn$row == "CINA3", "deployment_hot_cool",
                              ifelse(df.turn$col == "CINA2" & df.turn$row == "CINA4", "deployment_hot_cool",
                                                   ifelse(df.turn$same_value == "Yes", "same_deployment_season",
                                                          NA)))
   
   df.turn.comp.deploy <- df.turn[!is.na(df.turn$comp_deploy),]
   df.turn.comp.deploy <- df.turn.comp.deploy[1:42,]
   
   deploy_turn <- df.turn.comp.deploy$value
   
   #nestedness
   mat.nest
   df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
   df.nest <- subset(df.nest, row != col)
   
   df.nest$row <- substr(df.nest$row, 1, 5)
   df.nest$col <- substr(df.nest$col, 1, 5)
   df.nest$same_value <- ifelse(df.nest$row == df.nest$col, "Yes", "No")
   # df.nest <- subset(df.nest, row != col)
   
   df.nest$comp_deploy <- ifelse(df.nest$col == "CINA1" & df.nest$row == "CINA3", "deployment_hot_cool",
                              ifelse(df.nest$col == "CINA2" & df.nest$row == "CINA4", "deployment_hot_cool",
                                     ifelse(df.nest$same_value == "Yes", "same_deployment_season",
                                            NA)))
   
   df.nest.comp.deploy <- df.nest[!is.na(df.nest$comp_deploy),]
   df.nest.comp.deploy <- df.nest.comp.deploy[1:42,]
   
   deploy_nest <- df.nest.comp.deploy$value
   
   #jaccard similarity
   mat.jacc
   df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
   df.jacc <- subset(df.jacc, row != col)
   
   df.jacc$row <- substr(df.jacc$row, 1, 5)
   df.jacc$col <- substr(df.jacc$col, 1, 5)
   df.jacc$same_value <- ifelse(df.jacc$row == df.jacc$col, "Yes", "No")
   # df.jacc <- subset(df.jacc, row != col)
   
   df.jacc$comp_deploy <- ifelse(df.jacc$col == "CINA1" & df.jacc$row == "CINA3", "deployment_hot_cool",
                              ifelse(df.jacc$col == "CINA2" & df.jacc$row == "CINA4", "deployment_hot_cool",
                                     ifelse(df.jacc$same_value == "Yes", "same_deployment_season",
                                            NA)))
   
   df.jacc.comp.deploy <- df.jacc[!is.na(df.jacc$comp_deploy),]
   df.jacc.comp.deploy <- df.jacc.comp.deploy[1:42,]
   
   deploy_jacc <- df.jacc.comp.deploy$value
   
   
   df_deploy <- data.frame(turn = deploy_turn,
                             nest = deploy_nest,
                             jacc = deploy_jacc,
                             comp_deploy = df.turn.comp.deploy$comp_deploy)
   
 
   p.sed <- rstatix::t_test(df_deploy, turn ~ comp_deploy, p.adjust.method = "bonferroni")
   p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
   depl = c("between same \n deployment season", "between different \n deployment season")
   g <- ggplot(df_deploy, aes(x = fct_relevel(comp_deploy, "same_deployment_season", "deployment_hot_cool"), y = turn)) +
     geom_boxplot(fill =  c("lightgreen","lightgreen") ) +
     labs(title = "",
          x = "Comparisons",
          y = "Turnover component") +
     theme(legend.position = "none") +
     scale_x_discrete(labels=depl) +
     theme_classic() +
     stat_pvalue_manual(p.sed) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   g
}






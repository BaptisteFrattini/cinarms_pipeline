#' permanova & simper
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return 
#' @export
#' 

fun_perm <- function(metadata_data_mean){
  
  
  devtools::load_all()
  library(vegan)
  library(reshape2)
  library(tidyr)
  library(dplyr)
  library(ggplot2)


  # metadata_data_mean = targets::tar_read(mean_metadata_data)

  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  df_mean_pa <- vegan::decostand(df_mean, "pa")
  rownames(df_mean_pa) <- meta_mean$arms
  
  #### Compute permanova ####
  
  meta_mean$imm_rec <- paste(meta_mean$immersion_season, meta_mean$recovery_season)
  
  
  # PERMANOVA ####
  
  perm.bray <- vegan::adonis2(df_mean ~ imm_time*immersion_season, data = meta_mean, method = "bray", permutations = 99999)
  perm.jacc <- vegan::adonis2(df_mean_pa ~ imm_time*immersion_season, data = meta_mean, method = "jaccard", permutations = 99999)
  
  # PERMDISP ####
  
  dis.bray <- vegan::vegdist(df_mean, "bray")
  dis.jacc <- vegan::vegdist(df_mean_pa, "jaccard")
  
  a <- vegan::betadisper(
    dis.bray,
    meta_mean$arms_name,
    type = "median",
    bias.adjust = FALSE,
    sqrt.dist = FALSE,
    add = TRUE
  )
  
  
  b <- vegan::betadisper(
    dis.jacc,
    meta_mean$arms_name,
    type = "median",
    bias.adjust = FALSE,
    sqrt.dist = FALSE,
    add = TRUE
  )
  
  # Pairwise PERMANOVA ####
  
  # pairwise_ado <- pairwise.adonis2(df_mean ~ imm_time, data = meta_mean, method = "bray")

  
  sim_imm_tim <- summary(vegan::simper(df_mean, meta_mean$imm_time))
  
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$p < 0.05) | (sim_imm_tim_1$average > 0.1),]
  
  sim_imm_tim_2 <- sim_imm_tim[[2]]
  contrib_imm_tim_2 <- sim_imm_tim_2[(sim_imm_tim_2$p < 0.05) | (sim_imm_tim_2$average > 0.1),]
  
  sim_imm_tim_3 <- sim_imm_tim[[3]]
  contrib_imm_tim_3 <- sim_imm_tim_3[(sim_imm_tim_3$p < 0.05) | (sim_imm_tim_3$average > 0.1),]
  
  levels(as.factor(c(rownames(contrib_imm_tim_1), rownames(contrib_imm_tim_2), rownames(contrib_imm_tim_3))))
  
  
  sim_imm_season <- summary(vegan::simper(df_mean, meta_mean$immersion_season))
  sim_imm_recovery <- summary(vegan::simper(df_mean, meta_mean$recovery_season))
  
  # Model pour tester l'effet inter ARMS ####
  
  rownames(df_mean) <- meta_mean$arms
  
  mat.bray <- vegan::vegdist(df_mean, "bray")
  
  df.bray <- reshape2::melt(as.matrix(mat.bray), varnames = c("row", "col"))
  
  df.bray$row <- as.character(df.bray$row)
  df.bray$col <- as.character(df.bray$col)
  df.bray <- subset(df.bray, row != col)
  df.bray <- df.bray %>%
              mutate(
              row = as.character(row),  # Convert 'row' to character
              col = as.character(col),  # Convert 'col' to character
              combined = paste(pmin(row, col), pmax(row, col), sep = "_")
              )
  
  df.bray <- df.bray[!duplicated(df.bray$combined), ]
  df.bray$ARMS1 <- substr(df.bray$row, 1, 6)
  df.bray$ARMS2 <- substr(df.bray$col, 1, 6)
  
  df.bray <- data.frame(ARMS1 = df.bray$ARMS1, 
                        ARMS2 = df.bray$ARMS2, 
                        value = as.numeric(df.bray$value))
  
  meta_mean$comb <- paste0(substr(meta_mean$imm_time,1,1), "_", substr(meta_mean$immersion_season,5,8))
  
  df.bray <- df.bray %>%
    left_join(meta_mean %>% select(arms, comb), by = c("ARMS1" = "arms")) %>%
    rename(Mod1 = comb) %>%
    left_join(meta_mean %>% select(arms, comb), by = c("ARMS2" = "arms")) %>%
    rename(Mod2 = comb)
  
  step = 0.025*(max(df.bray$value)-min(df.bray$value))
  brek = seq(min(df.bray$value),max(df.bray$value),  step)
  
  d = ggplot(df.bray, aes(x = value)) +
      geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
      labs(title = paste0("Frequency of Bray-Curtis values"), x = "BC index", y = "Frequency") +
      theme_minimal()
  
  # path_to_boxplot <- paste0("outputs/bray_distrib.pdf")
  # ggsave(filename =  path_to_boxplot, plot = d, width = 18.5, height = 7.5)
  
  # # Modèle linéaire mixte avec effets aléatoires
  # Model <- lme4::lmer(value ~ Mod2 + (1 | ARMS1), data = df.bray)
  # summary(Model)
  # anova(Model)
  # 
  # # Visualisation de l'effet du temps sur la similarité moyenne
  # ggplot(df.bray, aes(x = Mod2, y = value, group = ARMS1)) +
  #   geom_point() +
  #   geom_smooth(method = "lm", se = FALSE) +
  #   theme_minimal()
  # 
  # m2 <- lme4::lmer(value ~ Mod2 + (1 | ARMS1)+(1|ARMS2), df.bray, REML=FALSE)
  # m1 <- stats::update(m2,.~Mod2 + (1 | ARMS1))
  # m0 <- lm(sim ~ Mod2,mod.db)
  # anova(m2,m1,m0)
  # ## compare m0 and m1
  # exactLRT(m1,m0)
  
  return(NULL)  
    
}

  
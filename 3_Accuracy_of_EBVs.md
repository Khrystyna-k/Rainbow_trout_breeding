# Load required libraries
library(tidyverse)
library(ggpubr)
library(ggExtra)
library(gridExtra)
library(here)


# This section evaluates the accuracy of mid-parent Estimated Breeding Values (EBVs) for rainbow trout by comparing them with adjusted phenotypes (y_star) from year class 2019 (first and second measurements). Accuracy is visualized via scatter plots with marginal histograms.

# Load adjusted phenotype predictions (y_star) from the first measurement
yhat_1st <- read.table("yhat_residual_19_l_1st", header = TRUE)
colnames(yhat_1st) <- c("level", "y_star", "y_hat", "residual")

# Load corresponding pedigree and ID information
ped_pred <- read.table("renadd04_19_pr.ped", header = FALSE)
colnames(ped_pred) <- c("level", "2", "3", "remove", "5", "6", "7", "8", "9", "id")

# Merge to map individual IDs to y_star
adj_pheno_2019_1st <- yhat_1st %>% inner_join(ped_pred, by = "level") %>% select("id", "y_star")
adj_pheno_2019_f_1st <- adj_pheno_2019_1st %>% inner_join(id_fam_19, by = c("id" = "Tag_ID"))
adj_pheno_2019_f_1st <- adj_pheno_2019_f_1st[adj_pheno_2019_f_1st$Family != 0, ]

# Extract family and parental pedigree
PED_all <- read.table("RT_pedigree_16_19.txt", header = FALSE)
colnames(PED_all) <- c("id", "Sire", "Dam")
PED_19 <- PED_all[3990:5634, ]
PED_19_f <- PED_19 %>% inner_join(id_fam_19, by = c("id" = "Tag_ID"))
PED_19_f$Family <- as.numeric(as.factor(str_c(PED_19_f$Sire, PED_19_f$Dam, sep = ".")))

# Calculate mid-parent EBVs using EBVs from 2016
parents19_ebv16l_1st <- ebv_ped_16_TRl_1st[, c("id", "solution")]
parents19_ebv16l_1st2 <- PED_19_f %>%
  inner_join(parents19_ebv16l_1st, by = c("Sire" = "id")) %>%
  inner_join(parents19_ebv16l_1st, by = c("Dam" = "id"))
colnames(parents19_ebv16l_1st2)[5:6] <- c("Sire_EBV", "Dam_EBV")
parents19_ebv16l_1st2$Mid_EBV <- rowMeans(parents19_ebv16l_1st2[, c("Sire_EBV", "Dam_EBV")])

ystar_per_ebv_1st <- parents19_ebv16l_1st2 %>% inner_join(adj_pheno_2019_1st, by = "id")

# Repeat for second measurement
yhat_2 <- read.table("yhat_residual_19_l_2nd", header = FALSE)
colnames(yhat_2) <- c("level", "y_star", "y_hat", "residual")

adj_pheno_2019_2 <- yhat_2 %>% inner_join(ped_19_pred, by = "level") %>% select("id", "y_star")
adj_pheno_2019_f_2 <- adj_pheno_2019_2 %>% inner_join(id_fam_19, by = c("id" = "Tag_ID"))
adj_pheno_2019_f_2 <- adj_pheno_2019_f_2[adj_pheno_2019_f_2$Family != 0, ]

parents19_ebv16l <- ebv_ped_16_TRl[, c("id", "solution")]
parents19_ebv16l_2 <- PED_19_f %>%
  inner_join(parents19_ebv16l, by = c("Sire" = "id")) %>%
  inner_join(parents19_ebv16l, by = c("Dam" = "id"))
colnames(parents19_ebv16l_2)[5:6] <- c("Sire_EBV", "Dam_EBV")
parents19_ebv16l_2$Mid_EBV <- rowMeans(parents19_ebv16l_2[, c("Sire_EBV", "Dam_EBV")])

ystar_per_ebv_2 <- parents19_ebv16l_2 %>% inner_join(adj_pheno_2019_f_2, by = "id")

# Compute mean adjusted phenotype and EBVs per family (first measurement)
means_1st <- ystar_per_ebv_1st %>%
  group_by(Family) %>%
  summarise(
    y_star_mean = mean(y_star, na.rm = TRUE),
    Sire_EBV = mean(Sire_EBV),
    Dam_EBV = mean(Dam_EBV),
    Mid_EBV = mean(Mid_EBV)
  )
cor(means_1st$y_star_mean, means_1st$Mid_EBV, method = "spearman")

# Compute for second measurement
means <- ystar_per_ebv_2 %>%
  group_by(Family.y) %>%
  summarise(
    y_star_mean = mean(y_star, na.rm = TRUE),
    Sire_EBV = mean(Sire_EBV),
    Dam_EBV = mean(Dam_EBV),
    Mid_EBV = mean(Mid_EBV)
  )
cor(means$y_star_mean, means$Mid_EBV, method = "spearman")

# Plot correlation for first measurement
g_1st <- ggplot(means_1st, aes(Mid_EBV, y_star_mean)) +
  geom_point(shape = 21, fill = "darkgreen", color = "black", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text", label = "Correlation: 0.41, p < 0.01", x = -1.6, y = 35, size = 5, colour = "black", fontface = "italic") +
  annotate("text", label = "Accuracy: 0.66", x = -3.5, y = 31, size = 5, colour = "black", fontface = "italic") +
  theme_pubr() +
  theme(legend.position = "none", text = element_text(size = 14)) +
  labs(x = "\n Mid-parental EBV", y = "Mean progeny phenotype (TL, mm)", title = "A) Year class 2016")
g_1st <- ggMarginal(g_1st, type = "histogram", fill = "grey", xparams = list(bins = 10))

# Plot correlation for second measurement
g_2 <- ggplot(means, aes(Mid_EBV, y_star_mean)) +
  geom_point(shape = 21, fill = "darkgreen", color = "black", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text", label = "Correlation: 0.31, p < 0.01", x = 3.9, y = 35, size = 5, colour = "black", fontface = "italic") +
  annotate("text", label = "Accuracy: 0.65", x = 2, y = 30, size = 5, colour = "black", fontface = "italic") +
  theme_pubr() +
  theme(legend.position = "none", text = element_text(size = 14)) +
  labs(x = "\n Mid-parental EBV", y = "Mean progeny phenotype (TL, mm)", title = "B) Year class 2019")
g_2 <- ggMarginal(g_2, type = "histogram", fill = "grey", xparams = list(bins = 10))

# Combine and export the accuracy plots
acc <- gridExtra::grid.arrange(g_1st, g_2, ncol = 2)
ggsave(path = " ", "Accuracy.png", acc, width = 11, height = 6, dpi = 300)

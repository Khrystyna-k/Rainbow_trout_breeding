# Plot means for fecundity and egg survival traits in generation 0 (G0) and generation 1 (G1)

```{r}
# Load libraries
library(tidyverse)
library(ggpubr)
library(viridis)
library(here)
library(psych)


# Load data
surv_data <- read_tsv(here("data/pheno_survival.txt"))

# Recode year to generation
surv_data <- surv_data %>%
  mutate(Generation = recode(as.character(Year), `2016` = "G0", `2019` = "G1"))

# Quick check
glimpse(surv_data)

# Plotting fecundity and survival
plot_fecundity <- ggboxplot(
  data = surv_data,
  x = "Generation", y = "Fecundity", fill = "Generation"
) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.7) +
  theme_pubr() +
  theme(legend.position = 'none') +
  labs(x = "Generation", y = "Fecundity (No. of eggs/female)")

plot_survival <- ggboxplot(
  data = surv_data,
  x = "Generation", y = "Egg_survival_perc", fill = "Generation"
) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.7) +
  theme_pubr() +
  theme(legend.position = 'none') +
  labs(x = "Generation", y = "Egg survival (%)")

# Combine and export
combined_plot <- ggarrange(plot_fecundity, plot_survival, nrow = 1)
ggsave(here("results/Fecundity_Survival.pdf"), combined_plot, width = 8, height = 4, dpi = 300)

```

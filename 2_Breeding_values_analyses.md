
# Analyses of breeding values estimated with blupf90
```{r}
# Load required packages
library(tidyverse)
library(here)



# Helper Function: Load and filter EBVs
load_and_filter_ebv <- function(path, effect_value) {
  read_table(here("data", path), col_names = TRUE) %>%
    filter(effect == effect_value)
}

# Helper Function: Add pedigree and filter individuals
merge_pedigree <- function(ebv_df, pedigree_df, keep_col = "remove", keep_val = 1) {
  ebv_df %>%
    inner_join(pedigree_df, by = "level") %>%
    filter(!!sym(keep_col) == keep_val) %>%
    select(id, solution, `s.e.`)
}


# Load pedigree files
ped_16 <- read_table(here("data", "renadd05_16.ped"), col_names = FALSE)
ped_19 <- read_table(here("data", "renadd04_19.ped"), col_names = FALSE)

colnames(ped_16) <- colnames(ped_19) <- c("level", "V2", "V3", "remove", "V5", "V6", "keep", "V8", "V9", "id")

# Read EBVs for 2016 and 2019 (length, weight, 1st measurement)
ebv_16_l_1 <- load_and_filter_ebv("solutions_16_l_1st", effect_value = 5)
ebv_19_l_1 <- load_and_filter_ebv("solutions_19_l_1st", effect_value = 4)

ebv_16_w_1 <- load_and_filter_ebv("solutions_16_w_1st", effect_value = 5)
ebv_19_w_1 <- load_and_filter_ebv("solutions_19_w_1st", effect_value = 4)

# Join with pedigree
ebv_16_l_1_clean <- merge_pedigree(ebv_16_l_1, ped_16)
ebv_19_l_1_clean <- merge_pedigree(ebv_19_l_1, ped_19, keep_col = "keep")

ebv_16_w_1_clean <- merge_pedigree(ebv_16_w_1, ped_16)
ebv_19_w_1_clean <- merge_pedigree(ebv_19_w_1, ped_19, keep_col = "keep")

# Save cleaned files
write_table_safe <- function(df, filename) {
  write.table(df, here("results", filename), col.names = TRUE, row.names = FALSE, quote = FALSE)
}

write_table_safe(ebv_16_l_1_clean, "ebv_16_length_1st.txt")
write_table_safe(ebv_16_w_1_clean, "ebv_16_weight_1st.txt")
write_table_safe(ebv_19_l_1_clean, "ebv_19_length_1st.txt")
write_table_safe(ebv_19_w_1_clean, "ebv_19_weight_1st.txt")

# Load libraries
library(tidyverse)
library(ggpubr)
library(psych)
library(here)

# Define input filenames
length_files <- list(
  ped_ebv_16_19_l_1_f = "ped_ebv_16_19_l_1_f.txt",
  ped_ebv_16_19_l_2_f = "ped_ebv_16_19_l_2_f.txt",
  ped_ebv_16_l_1_f = "ped_ebv_16_l_1_f.txt",
  ped_ebv_16_l_2_f = "ped_ebv_16_l_2_f.txt",
  ped_ebv_19_l_1_f = "ped_ebv_19_l_1_f.txt",
  ped_ebv_19_l_2_f = "ped_ebv_19_l_2_f.txt"
)

weight_files <- list(
  ped_ebv_16_19_w_1_f = "ped_ebv_16_19_w_1_f.txt",
  ped_ebv_16_19_w_2_f = "ped_ebv_16_19_w_2_f.txt",
  ped_ebv_16_w_1_f = "ped_ebv_16_w_1_f.txt",
  ped_ebv_16_w_2_f = "ped_ebv_16_w_2_f.txt",
  ped_ebv_19_w_1_f = "ped_ebv_19_w_1_f.txt",
  ped_ebv_19_w_2_f = "ped_ebv_19_w_2_f.txt"
)

# Read all tables and bind into a named list
read_ebv_files <- function(file_list) {
  map(file_list, ~ read.table(here("data", .x), header = TRUE))
}

ebv_length <- read_ebv_files(length_files)
ebv_weight <- read_ebv_files(weight_files)

# Add generation and age info
add_metadata <- function(df, generation, age) {
  df %>% mutate(Generation = generation, Age = age)
}

# Annotate each EBV table
ebv_length <- list(
  l16_1 = add_metadata(ebv_length$ped_ebv_16_l_1_f, "G0", "9 months"),
  l16_2 = add_metadata(ebv_length$ped_ebv_16_l_2_f, "G0", "24 months"),
  l19_1 = add_metadata(ebv_length$ped_ebv_19_l_1_f, "G1", "9 months"),
  l19_2 = add_metadata(ebv_length$ped_ebv_19_l_2_f, "G1", "27 months"),
  l16_19_1 = add_metadata(ebv_length$ped_ebv_16_19_l_1_f, "Combined G0 and G1", "9 months"),
  l16_19_2 = ebv_length$ped_ebv_16_19_l_2_f %>%
    mutate(Generation = "Combined G0 and G1",
           Age = case_when(
             str_starts(id, "G2") ~ "24 months",
             str_starts(id, "G3") ~ "27 months",
             TRUE ~ NA_character_
           ))
)

length_all <- bind_rows(ebv_length)

# Group Means by Family
length_summary <- length_all %>%
  group_by(Family, Generation, Age) %>%
  summarise(mean_sol_length = round(mean(solution), 2), .groups = "drop")

# Plotting function
plot_family_histogram <- function(df, gen_label, palette, y_lim = c(-30, 20)) {
  ggplot(df, aes(x = reorder(Family, mean_sol_length), y = mean_sol_length, fill = Age)) +
    geom_col(position = "dodge") +
    theme_classic(base_size = 14) +
    scale_fill_manual(values = palette) +
    labs(x = "Family", y = "Estimated breeding value", fill = "Age") +
    ggtitle(paste("EBVs for", gen_label)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(y_lim)
}

# Generate plots

cbPalette <- c("orangered1", "#0072B2")

g0_plot <- length_summary %>% filter(Generation == "G0") %>% plot_family_histogram("G0", cbPalette)
g1_plot <- length_summary %>% filter(Generation == "G1") %>% plot_family_histogram("G1", cbPalette)
combined_plot <- length_summary %>% filter(Generation == "Combined G0 and G1") %>% plot_family_histogram("Combined", cbPalette)


# Save plots
ggsave(here("results/ebv_2016.png"), g0_plot, width = 13, height = 12, dpi = 300)
ggsave(here("results/ebv_2019.png"), g1_plot, width = 13, height = 12, dpi = 300)
ggsave(here("results/ebv_combined.png"), combined_plot, width = 13, height = 12, dpi = 300)


# Descriptive statistics by Age and Generation
describeBy(length_summary %>% filter(Generation == "G0"), group = "Age")
describeBy(length_summary %>% filter(Generation == "G1"), group = "Age")
describeBy(length_summary %>% filter(Generation == "Combined G0 and G1"), group = "Age")

```



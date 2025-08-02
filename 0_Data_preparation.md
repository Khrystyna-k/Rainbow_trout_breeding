# scripts/plot_parent_growth.R

# Load required libraries
library(tidyverse)
library(ggpubr)
library(ggridges)
library(ggbeeswarm)


#Standartized data
#Data growth of both YC
growth16_19 <- read.table("C:/Users/rabu0002/Desktop/blupf/New calculations/rt.growth16_19.txt", header = FALSE, sep = "")

colnames(growth16_19) <- c("id", "dat_1", "lenght_1", "weight_1", "tank", "dat_2", "lenght_2", "weight_2", "pond", "yearclass")

growth16_19[growth16_19 == 0] <- NA

growth16_norm <- growth16_19 %>% filter(yearclass == 16) %>% select("id", "dat_1", "lenght_1", "weight_1", "tank", "dat_2", "lenght_2", "weight_2", "pond", "yearclass")

growth19_norm <- growth16_19 %>% filter(yearclass == 19) %>% select("id", "dat_1", "lenght_1", "weight_1", "tank", "dat_2", "lenght_2", "weight_2", "pond", "yearclass")
  
#Normalize the data
growth16_norm$length_scaled<-round(scale(growth16_norm$lenght_2),2)
growth16_norm$weight_scaled<-round(scale(growth16_norm$weight_2),2)

growth19_norm$length_scaled<-round(scale(growth19_norm$lenght_2),2)
growth19_norm$weight_scaled<-round(scale(growth19_norm$weight_2),2)

#Join for blupf90
growth16_19_norm <- growth16_norm %>% 
  bind_rows(growth19_norm)


write.table(growth16_19_norm, file =  "growth16_19_norm.txt", append = FALSE, sep =" ", dec = ".",row.names = F, col.names = F)




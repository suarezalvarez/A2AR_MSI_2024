library(ggplot2)
library(tidyverse)

dist_adn = read.table("data/DISTANCE_SMD/TRP246CZ_vs_ADNO5.txt",
                      col.names = c("frame", "distance"))

dist_mod_2 = read.table("data/DISTANCE_SMD/MOD_2_TRP246CZ3_vs_UNKO5.txt",
                      col.names = c("frame", "distance"))

dist_mod_2 = dist_mod_2 %>% filter(frame < 101)

dist_mod_1 = read.table("data/DISTANCE_SMD/MOD_1_TRP246CZ3_vs_UNKO5.txt",
                          col.names = c("frame", "distance"))

dist_mod_1 = dist_mod_1 %>% filter(frame < 101)




global_dist = cbind(dist_adn, dist_mod_1$distance, dist_mod_2$distance)

colnames(global_dist) = c("frame", "ADN", "Modification 1" , "Modification 2")

global_dist = global_dist %>% pivot_longer(cols = -frame, names_to = "Ligand", values_to = "Distance")

global_dist_plot = global_dist %>% 
  ggplot() +
  aes(x = frame, y = Distance, color = Ligand) +
  geom_line() +
  labs(
       x = "Frame",
       y = "Distance between TRP246 CZ and O5 from ligands (Å)",
       color = "") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c('magenta','orange','cyan3'))



# save plot

jpeg('../../project/distance_SMD.jpeg',
     res = 500,
     units = "in",
     width = 6,
     height = 4)

global_dist_plot 
dev.off()  

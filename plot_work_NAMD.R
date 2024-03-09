library(ggplot2)
library(tidyverse)


# read files
w = readLines('data/PMF/w_ADN.txt')
c = readLines('data/PMF/c_ADN.txt')

# comma separated string
w = gsub(' ' , ',', w)
c = gsub(' ' , ',' , c)


# vectors
w = as.numeric(strsplit(w, ',')[[1]])

c = as.numeric(strsplit(c, ',')[[1]])


# data frame
w_vs_c = data.frame(work = w, extension = c)

w_vs_c$ligand = "ADN"

w_vs_c = w_vs_c %>% filter(extension < 53.7)


### mod_1

# read files

w_mod_1 = readLines('data/PMF/w_mod_1.out')
c_mod_1 = readLines('data/PMF/c_mod_1.out')

# comma separated string

w_mod_1 = gsub(' ' , ',', w_mod_1)
c_mod_1 = gsub(' ' , ',' , c_mod_1)


# vectors

w_mod_1 = as.numeric(strsplit(w_mod_1, ',')[[1]])

c_mod_1 = as.numeric(strsplit(c_mod_1, ',')[[1]])

w_mod_1 = w_mod_1[1:600]

# data frame

w_vs_c_mod_1 = data.frame(work = w_mod_1, extension = c_mod_1)

w_vs_c_mod_1$ligand = "Modification 1"
w_vs_c_mod_1 = w_vs_c_mod_1 %>% filter(extension < 53.7)


### mod_2

# read files
w_mod_2 = readLines('data/PMF/w_mod_2.out')
c_mod_2 = readLines('data/PMF/c_mod_2.out')

# comma separated string
w_mod_2 = gsub(' ' , ',', w_mod_2)
c_mod_2 = gsub(' ' , ',' , c_mod_2)


# vectors
w_mod_2 = as.numeric(strsplit(w_mod_2, ',')[[1]])

c_mod_2 = as.numeric(strsplit(c_mod_2, ',')[[1]])
w_mod_2 = w_mod_2[1:600]

# data frame
w_vs_c_mod_2 = data.frame(work = w_mod_2, extension = c_mod_2)

w_vs_c_mod_2 = w_vs_c_mod_2 %>% filter(extension < 53.7)
w_vs_c_mod_2$ligand = "Modification 2"






# plot


w_vs_c_merged = rbind(w_vs_c, w_vs_c_mod_1,w_vs_c_mod_2)

w_vs_c_merged %>% 
  ggplot(.,  aes(x = extension , y = work,
                 col = ligand)) +
  geom_line(lwd = .5,
            linetype = 1) +
  xlab("Extension (Å)") +
  ylab(expression('Work (kcal)')) +
  scale_x_continuous(n.breaks = 10) +
  theme_minimal() +
  scale_color_manual(values = c('magenta','orange','cyan3')) +
  
  theme(text = element_text(family = "Sans Sharif"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  
  labs(color = "Ligand")
  
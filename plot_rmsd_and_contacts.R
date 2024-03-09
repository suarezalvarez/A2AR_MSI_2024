library(ggplot2)
library(tidyverse)
library(gganimate)
library(patchwork)
library(ggrepel)
#######################
### ORIGINAL LIGAND ###
#######################

### RMSD over time

ligand_rmsd = read.delim('data/RMSD/ProtAlign_ADN.rmsd', 
                         header=T,
                         sep = '',
                         col.names = c('frame' , 'rmsd'))


rmsd_lig_plot = ligand_rmsd %>% 
  ggplot() +
  aes(x = frame, y = rmsd) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.83 , linetype = 'dashed', color = 'red') +
  xlim(c(0,17)) +
  xlab('Frame') + 
  ylab('RMSD (A)') +
  
  geom_text(aes(x = 13.5, y = 0.84, label = 'Average RMSD (0.83 A)'), 
            hjust = 0, vjust = 0, size = 2.5, color = 'red') +
  
  theme_minimal() 

  


jpeg('project/rmsd_ligand.jpeg' , units = 'in',
     res = 500, width = 8, height = 4)
rmsd_lig_plot
dev.off()


### Contacts over time

lines = readLines('data/CONTACTS/contacts_3A_ADN.txt')

# Read the file line by line


# Initialize an empty data frame
contacts_ADN <- data.frame(frame = 0:17)

# Loop over the lines
for (i in seq_along(lines)) {
  line <- lines[i]
  
  # Check if the line starts with "freeSelLabel"
  if (startsWith(line, "freeSelLabel")) {
    # Extract the label from the line
    label <- str_extract(line, "\\w+ \\d+")
    
    # Initialize a vector for the data
    data <- integer(18)
    
    # Loop over the next 17 lines and add the data to the vector
    for (j in 1:18) {
      data[j] <- as.integer(str_split(lines[i + j + 1], " ")[[1]][2])
    }
    
    # Add the data to the data frame
    contacts_ADN[[label]] <- data
  }
}

contacts_ADN = contacts_ADN %>% select(!contains('ADN'))



# Pivot  longer

contacts_ADN = contacts_ADN %>% 
  pivot_longer(cols = -frame,names_to = 'residue', values_to = 'contacts')



contacts_plot = contacts_ADN %>% 
  ggplot(.) + 
  aes(x = frame , y = contacts) +
  geom_line() +
  geom_point(size = 2) +
  xlab('Frame') +
  ylab('Contacts') +
  theme_minimal() +
  facet_wrap(~residue) +
  theme(axis.title = element_text(size = 15),
        text = element_text(family = "Sans sharif")) 

  # animation


contacts_plot + 
  transition_reveal(frame) +
  ease_aes('linear') 

  

animate(contacts_plot , fps = 20)
  
animate(contacts_plot, duration = 5, fps = 20, width = 600, height = 400, renderer = gifski_renderer())
anim_save('project/contacts.gif') 



# binary and percent

contacts_ADN_bin = contacts_ADN %>% 
  mutate(contacts = ifelse(contacts > 0, 1, 0)) %>% 
  group_by(residue) %>% 
  summarise(contacts = sum(contacts)) %>% 
  mutate(contacts = contacts/18*100) 



jpeg('project/contacts_binary_ADN.jpeg', units = 'in',
     res = 500, width = 4, height = 8)

ggplot(contacts_ADN_bin) + 
  aes(x = reorder(residue, -contacts), y = contacts) +
  geom_point(size = 8,
             alpha = .7) +
  xlab('Residue') +
  ylab('Contacts (%)') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 17)) +
  scale_y_continuous(breaks = seq(0,100,10)) +
  coord_flip()


dev.off()










#############################
### MODIFIED LIGAND mod_1 ###
#############################


### RMSD over time

mod_mod_1_rmsd = read.delim('data/RMSD/ProtAlign_mod_1.rmsd', 
                         header=T,
                         sep = '',
                         col.names = c('frame' , 'rmsd'))
                         
(mod_1_rmsd_plot = mod_mod_1_rmsd %>% 
  ggplot() +
  aes(x = frame, y = rmsd) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 1.16 , linetype = 'dashed', color = 'red') +
  xlim(c(0,17)) +
  xlab('Frame') +
  ylab('RMSD (A)') +
  theme_minimal() +
  geom_text(aes(x = 10 ,y = 1.17, label = 'Average RMSD (1.15 A)'), 
            hjust = 0, vjust = 0, size = 2.5, color = 'red'))


jpeg('project/mod_1_rmsd.jpeg' , units = 'in',
     res = 500, width = 8, height = 4)
mod_1_rmsd_plot
dev.off()



### Contacts over time

lines_mod_1 = readLines('data/CONTACTS/contacts_3A_mod_1.txt')



# Initialize an empty data frame
contacts_UNK_mod_1 <- data.frame(frame = 0:17)

# Loop over the lines
for (i in seq_along(lines_mod_1)) {
  line <- lines_mod_1[i]
  
  # Check if the line starts with "freeSelLabel"
  if (startsWith(line, "freeSelLabel")) {
    # Extract the label from the line
    label <- str_extract(line, "\\w+ \\d+")
    
    # Initialize a vector for the data
    data <- integer(18)
    
    # Loop over the next 17 lines and add the data to the vector
    for (j in 1:18) {
      data[j] <- as.integer(str_split(lines_mod_1[i + j + 1], " ")[[1]][2])
    }
    
    # Add the data to the data frame
    contacts_UNK_mod_1[[label]] <- data
  }
}


contacts_UNK_mod_1 = contacts_UNK_mod_1 %>% select(!contains('UNK'))

contacts_UNK_mod_1 = contacts_UNK_mod_1 %>% 
  pivot_longer(cols = -frame,names_to = 'residue', values_to = 'contacts')



contacts_mod_1_plot = contacts_UNK_mod_1 %>% 
  ggplot(.) + 
  aes(x = frame , y = contacts) +
  geom_line() +
  geom_point(size = 2) +
  xlab('Frame') +
  ylab('Contacts') +
  theme_minimal() +
  facet_wrap(~residue) +
  theme(axis.title = element_text(size = 15),
        text = element_text(family = "Sans sharif"))
  
  # animation (add plus sign at the end of previous line to continue the code in the next line)

  contacts_mod_1_plot +
  transition_reveal(frame) +
  ease_aes('linear')

animate(contacts_mod_1_plot, duration = 5, fps = 20, width = 600, height = 400, renderer = gifski_renderer())
anim_save('project/contacts_mod_1.gif') 


# binary and percent

contacts_mod_1_bin = contacts_UNK_mod_1 %>% 
  mutate(contacts = ifelse(contacts > 0, 1, 0)) %>% 
  group_by(residue) %>% 
  summarise(contacts = sum(contacts)) %>% 
  mutate(contacts = contacts/18*100)


jpeg('project/contacts_binary_mod_1.jpeg', units = 'in',
     res = 500, width = 4, height = 8)

ggplot(contacts_mod_1_bin) + 
  aes(x = reorder(residue, -contacts), y = contacts) +
  geom_point(size = 8,
             alpha = .7) +
  xlab('Residue') +
  ylab('Contacts (%)') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 17)) +
  scale_y_continuous(breaks = seq(0,100,10)) +
  coord_flip()


dev.off()








###############################
### MODIFIED LIGAND mod_2 ###
###############################



### RMSD over time

mod_2_rmsd = read.delim('data/RMSD/ProtAlign_mod_2.rmsd', 
                         header=T,
                         sep = '',
                         col.names = c('frame' , 'rmsd'))



(mod_2_rmsd_plot = mod_2_rmsd %>% 
  ggplot() +
  aes(x = frame, y = rmsd) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 1.24 , linetype = 'dashed', color = 'red') +
  xlim(c(0,17)) +
  xlab('Frame') +
  ylab('RMSD (A)') +
  theme_minimal() +
  geom_text(aes(x = 10 ,y = 1.25, label = 'Average RMSD (1.24 A)'), 
              hjust = 0, vjust = 0, size = 2.5, color = 'red'))


jpeg('project/mod_2_rmsd.jpeg' , units = 'in',
     res = 500, width = 8, height = 4)
mod_2_rmsd_plot
dev.off()



### Contacts over time

lines_mod_2 = readLines('data/CONTACTS/contacts_3A_mod_2.txt')

# Initialize an empty data frame

contacts_UNK_mod_2 <- data.frame(frame = 0:17)

# Loop over the lines

for (i in seq_along(lines_mod_2)) {
  line <- lines_mod_2[i]
  
  # Check if the line starts with "freeSelLabel"
  if (startsWith(line, "freeSelLabel")) {
    # Extract the label from the line
    label <- str_extract(line, "\\w+ \\d+")
    
    # Initialize a vector for the data
    data <- integer(18)
    
    # Loop over the next 17 lines and add the data to the vector
    for (j in 1:18) {
      data[j] <- as.integer(str_split(lines_mod_2[i + j + 1], " ")[[1]][2])
    }
    
    # Add the data to the data frame
    contacts_UNK_mod_2[[label]] <- data
  }
}


contacts_UNK_mod_2 = contacts_UNK_mod_2 %>% select(!contains('UNK'))


contacts_UNK_mod_2 = contacts_UNK_mod_2 %>% pivot_longer(-frame,
                                                             names_to = 'residue',
                                                             values_to = 'contacts')
contacts_UNK_mod_2_plot = contacts_UNK_mod_2 %>%
  ggplot(.) + 
  aes(x = frame , y = contacts) +
  geom_line() +
  geom_point(size = 2) +
  xlab('Frame') +
  ylab('Contacts') +
  theme_minimal() +
  facet_wrap(~residue) +
  scale_y_continuous(n.breaks = 4) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15),
        text = element_text(family = "Sans sharif"))



  # animation (add plus sign at the end of previous line to continue the code in the next line)
  
contacts_UNK_mod_2_plot +
  transition_reveal(frame) +
  ease_aes('linear') 



animate(contacts_UNK_mod_2_plot, duration = 4, fps = 20, width = 600, height = 400, renderer = gifski_renderer())
anim_save('project/contacts_mod_2.gif') 




# binary and percent

contacts_mod_2_bin = contacts_UNK_mod_2 %>% 
  mutate(contacts = ifelse(contacts > 0, 1, 0)) %>% 
  group_by(residue) %>% 
  summarise(contacts = sum(contacts)) %>% 
  mutate(contacts = contacts/18*100)


jpeg('project/contacts_binary_mod_2.jpeg', units = 'in',
     res = 500, width = 4, height = 8)
ggplot(contacts_mod_2_bin) + 
  aes(x = reorder(residue, -contacts), y = contacts) +
  geom_point(size = 8,
             alpha = .7) +
  xlab('Residue') +
  ylab('Contacts (%)') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 17)) +
  scale_y_continuous(breaks = seq(0,100,10)) +
  ylim(c(0,100)) 

dev.off()







####################
### GLOBAL PLOTS ###
####################


# gathered contacts plot

global_contacts = rbind(contacts_ADN %>% mutate(ligand = "ADN"),
                        contacts_UNK_mod_1 %>% mutate(ligand = "Modification 1"),
                        contacts_UNK_mod_2 %>% mutate(ligand = "Modification 2"))



global_contacts_plot = global_contacts %>% 
  ggplot() + 
  aes(x = frame , y = contacts, col = ligand) +
  
  geom_line(lwd = .7) +
 
  xlab('Frame') +
  ylab('Number of atoms\nwithin 4 Å of the ligand') +
  
  theme_minimal() +
  
  facet_wrap(~residue,
             ncol = 6,
             nrow = 3) +
  
  scale_y_continuous(n.breaks = 4) +
  
  theme(legend.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 9),
        text = element_text(family = "Sans sherif"),
        legend.position = "bottom") +
  
  scale_color_manual(values = c('magenta','orange','cyan3')) +
  labs(col = '') 








# gathered rmsd plot

global_rmsd = rbind(ligand_rmsd %>% mutate(ligand = "ADN"),
                    mod_mod_1_rmsd %>% mutate(ligand = "Modification 1"),
                    mod_mod_2_rmsd %>% mutate(ligand = "Modification 2"))





global_rmsd_plot = global_rmsd %>%
  ggplot() +
  aes(x = frame, y = rmsd, col = ligand) +
  
  geom_line(lwd = 1) +
  geom_point(size = 2) +
  
  geom_hline(yintercept = 0.83 , linetype = 'dashed', color = 'magenta') +
  geom_hline(yintercept = 1.15 , linetype = 'dashed', color = 'orange') +
  geom_hline(yintercept = 1.24 , linetype = 'dashed', color = 'cyan3') +
  
  geom_text(aes(x = 1.5 ,y = .84, label = 'Average RMSD (0.83 Å)'), 
            hjust = 0, vjust = 0, size = 2, color = 'magenta') +
  geom_text(aes(x = 9.8 ,y = 1.16, label = 'Average RMSD (1.15 Å)'), 
            hjust = 0, vjust = 0, size = 2, color = 'orange') +
  geom_text(aes(x = 9.8 ,y = 1.25, label = 'Average RMSD (1.24 Å)'), 
            hjust = 0, vjust = 0, size = 2, color = 'cyan3') +
  
  xlab('Frame') +
  ylab('RMSD (Å)') +
  
  theme_minimal() +
  
  theme(legend.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 9),
        text = element_text(family = "Sans sherif"),
        legend.position = "bottom") +

  scale_color_manual(values = c('magenta','orange','cyan3')) +
  
  xlim(c(0,16)) +
  labs(col = "") 






# gather binary results


global_contacts_bin = rbind(contacts_ADN_bin %>% mutate(ligand = "ADN"),
                            contacts_mod_1_bin %>% mutate(ligand = "Modification 1"),
                            contacts_mod_2_bin %>% mutate(ligand = "Modification 2"))



global_contacts_bin_plot = global_contacts_bin %>% 
  ggplot() +
  aes(x = contacts, 
      y = residue,
      shape = ligand,
      size = ligand) +
  
  geom_point() +
  
  theme(legend.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 9),
        text = element_text(family = "Sans sherif"),
        legend.position = "bottom") +
  
  scale_shape_manual(values = c(0,1,2)) +
  scale_size_manual(values = c(3,5,7)) +

  
  labs(col = '',
       shape = '',
       size = '') +
  
  xlab('Contacts (%)') +
  ylab('Residue') +
  
  theme_minimal()
  


# save plots

jpeg('project/global_contacts.jpeg', units = 'in',
     res = 500, width = 6, height = 5)

global_contacts_plot

dev.off()




jpeg('project/global_rmsd.jpeg', units = 'in',
     res = 500, width = 6, height = 4)

global_rmsd_plot

dev.off()



jpeg('project/global_contacts_bin_plot', units = 'in',
     res = 500, width = 8, height = 7)

global_contacts_bin_plot

dev.off()



combined <- global_contacts_plot + global_contacts_bin_plot


jpeg('project/combined_pannel.jpeg', units = 'in',
     res = 500, width = 8, height = 12)
combined + plot_layout(ncol = 1, guides = "collect") + plot_annotation(tag_levels = "A")
dev.off()



library(vegan)
library(ggplot2)

pc = read.csv("Compositional_pattern/compositional_pattern_table.txt", header= TRUE, sep="\t")


# Step 1 Plot all seasons
## Make community matrix - extract columns with abundance information
com = pc[,3:ncol(pc)]
## Turn abundance data frame into a matrix
m_com = as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds
plot(nmds)
## Conduct ANOSIM Test
ano = anosim(m_com, pc$season, distance = "bray", permutations = 9999)
ano

## Extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds)$sites)

## Add columns to data frame 
data.scores$Sample = pc$head
data.scores$Season = pc$season

head(data.scores)

## Define the desired order of the Season levels
season_order <- c("Spring", "Clearwater", "Early Summer", "Late Summer", "Fall", "Ice-on")

xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = Season), alpha = 0.75) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Season", y = "NMDS2")  + 
  scale_colour_manual(values = c("#C8C5E2", "#A7C4D2", "#D2AE6D", "#F5874F", "#84716B", "#00FFFF"), 
                      breaks = season_order, 
                      labels = season_order, 
                      limits = season_order, 
                      drop = FALSE)

xx


# Step 2 Plot seasons of "Fall", "Ice-on", "Spring"
pc1 <- subset(pc, pc[, 2] %in% c("Ice-on", "Spring", "Fall"))
#pc1 <- subset(pc1, !(pc1[, 1] %in% c("3300042313", "3300042358"))) # Exclude two outlier samples
## Make community matrix - extract columns with abundance information
com1 = pc1[,3:ncol(pc)]
## Turn abundance data frame into a matrix
m_com1 = as.matrix(com1)
set.seed(123)
nmds1 = metaMDS(m_com1, distance = "bray")
nmds1
plot(nmds1)
## Conduct ANOSIM Test
ano1 = anosim(m_com1, pc1$season, distance = "bray", permutations = 9999)
ano1

## Extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores1 = as.data.frame(scores(nmds1)$sites)

## Add columns to data frame 
data.scores1$Sample = pc1$head
data.scores1$Season = pc1$season

head(data.scores1)

## Define the desired order of the Season levels
season_order <- c("Fall", "Ice-on", "Spring")

xx1 = ggplot(data.scores1, aes(x = NMDS1, y = NMDS2, label = Sample)) +
  geom_point(size = 4, aes(colour = Season), alpha = 0.75) +
  #geom_text(size = 3, show.legend = FALSE) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", colour = "Season", y = "NMDS2") +
  scale_colour_manual(values = c("#84716B", "#00FFFF", "#C8C5E2"),
                      breaks = season_order,
                      labels = season_order,
                      limits = season_order,
                      drop = FALSE)

xx1

# Set the dimensions for the PDF file
xx1 <- xx1 + coord_fixed(ratio = 1)

# Save xx1 plot to PDF
ggsave("Compositional_pattern/Rplot-F-I-S.pdf", xx1)


# Step 3 Plot seasons of "Spring", "Clearwater", and "Early Summer"
pc2 <- subset(pc, pc[, 2] %in% c("Spring", "Clearwater", "Early Summer"))
#pc2 <- subset(pc2, !(pc2[, 1] %in% c("3300042358"))) # Exclude one outlier sample
## Make community matrix - extract columns with abundance information
com2 = pc2[,3:ncol(pc)]
## Turn abundance data frame into a matrix
m_com2 = as.matrix(com2)
set.seed(123)
nmds2 = metaMDS(m_com2, distance = "bray")
nmds2
plot(nmds2)
## Conduct ANOSIM Test
ano2 = anosim(m_com2, pc2$season, distance = "bray", permutations = 9999)
ano2

## Extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores2 = as.data.frame(scores(nmds2)$sites)

## Add columns to data frame 
data.scores2$Sample = pc2$head
data.scores2$Season = pc2$season

head(data.scores2)

## Define the desired order of the Season levels
season_order <- c("Spring", "Clearwater", "Early Summer")

xx2 = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2, label = Sample)) +
  geom_point(size = 4, aes(colour = Season), alpha = 0.75) +
  #geom_text(size = 3, show.legend = FALSE) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", colour = "Season", y = "NMDS2") +
  scale_colour_manual(values = c("#C8C5E2", "#A7C4D2", "#D2AE6D"),
                      breaks = season_order,
                      labels = season_order,
                      limits = season_order,
                      drop = FALSE)

xx2

# Set the dimensions for the PDF file
xx2 <- xx2 + coord_fixed(ratio = 1)

# Save xx2 plot to PDF
ggsave("Compositional_pattern/Rplot-S-C-E.pdf", xx2)


# Step 4 Plot seasons of "Early Summer", "Late Summer", and "Fall"
pc3 <- subset(pc, pc[, 2] %in% c("Early Summer", "Late Summer", "Fall"))
#pc3 <- subset(pc3, !(pc3[, 1] %in% c("3300034094"))) # Exclude two outlier samples
## Make community matrix - extract columns with abundance information
com3 = pc3[,3:ncol(pc)]
## Turn abundance data frame into a matrix
m_com3 = as.matrix(com3)
set.seed(123)
nmds3 = metaMDS(m_com3, distance = "bray")
nmds3
plot(nmds3)
## Conduct ANOSIM Test
ano3 = anosim(m_com3, pc3$season, distance = "bray", permutations = 9999)
ano3

## Extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores3 = as.data.frame(scores(nmds3)$sites)

## Add columns to data frame 
data.scores3$Sample = pc3$head
data.scores3$Season = pc3$season

head(data.scores3)

## Define the desired order of the Season levels
season_order <- c("Early Summer", "Late Summer", "Fall")

xx3 = ggplot(data.scores3, aes(x = NMDS1, y = NMDS2, label = Sample)) +
  geom_point(size = 4, aes(colour = Season), alpha = 0.75) +
  #geom_text(size = 3, show.legend = FALSE) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", colour = "Season", y = "NMDS2") +
  scale_colour_manual(values = c("#D2AE6D", "#F5874F", "#84716B"),
                      breaks = season_order,
                      labels = season_order,
                      limits = season_order,
                      drop = FALSE)

xx3

# Set the dimensions for the PDF file
pdf_width <- 7.5
pdf_height <- 7.5
xx3 <- xx3 + coord_fixed(ratio = 1)

# Save xx3 plot to PDF
ggsave("Compositional_pattern/Rplot-E-L-F.pdf", xx3, width = pdf_width, height = pdf_height)

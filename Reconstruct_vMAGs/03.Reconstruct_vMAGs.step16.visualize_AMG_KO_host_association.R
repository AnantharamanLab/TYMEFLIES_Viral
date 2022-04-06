library(tidyverse)
library(ggplot2)
library(ggpmisc)

# Aim: Visulize the AMG KO to host association

# Step 1 Visualize the trend for table 1 (host is o__Cyanobacteriales;f__Nostocaceae)
## Read table
table <-read.table("./AMG_analysis/Viral_host_association/Viral_host_association_table_1.txt",sep = "\t",head=F, row.names = 1)
table <- t(table)
table <- as.data.frame(table)
# Change column names
colnames(table) <- c("Month","Host.family","K00507.desC","K01627.kdsA")

# Give the range of two y-axes
scaleFactor <- max(table$Host.family) / max(table$K00507.desC)
scaleFactor <- scaleFactor * 0.5

p1 <- ggplot(table,aes(x=Month)) + 
     geom_line(aes(y = Host.family,colour="Host.family"), size = 0.5)+
     geom_line(aes(y = K00507.desC*scaleFactor,colour="K00507.desC"), size = 0.5)+
     geom_line(aes(y = K01627.kdsA*scaleFactor,colour="K01627.kdsA"), size = 0.5)+
     scale_y_continuous(name="Host abundance", sec.axis=sec_axis(~./scaleFactor, name="Viral abundance")) + 
     scale_x_continuous(breaks=seq(1,12,1))+
     theme(
       legend.position="bottom",
     )+
     labs(colour="")
p1
ggsave(p1, file="./AMG_analysis/Viral_host_association/Viral_host_association_table_1.pdf", width = 10, height = 5, units = "cm")

# Step 2 Visualize the trend for table 2 (host is o__Flavobacteriales;f__Flavobacteriaceae)
## Read table
table2 <-read.table("./AMG_analysis/Viral_host_association/Viral_host_association_table_2.txt",sep = "\t",head=F, row.names = 1)
table2 <- t(table2)
table2 <- as.data.frame(table2)
# Change column names
colnames(table2) <- c("Month","Host.family","K00798.pduO","K01734.mgsA","K15895.pseC")

# Give the range of two y-axes
scaleFactor <- max(table2$Host.family) / max(table2$K01734.mgsA)
scaleFactor <- scaleFactor * 0.7

p2 <- ggplot(table2,aes(x=Month)) + 
  geom_line(aes(y = Host.family,colour="Host.family"), size = 0.5)+
  geom_line(aes(y = K00798.pduO*scaleFactor,colour="K00798.pduO"), size = 0.5)+
  geom_line(aes(y = K01734.mgsA*scaleFactor,colour="K01734.mgsA"), size = 0.5)+
  geom_line(aes(y = K15895.pseC*scaleFactor,colour="K15895.pseC"), size = 0.5)+
  scale_y_continuous(name="Host abundance", sec.axis=sec_axis(~./scaleFactor, name="Viral abundance")) + 
  scale_x_continuous(breaks=seq(1,12,1))+
  theme(
    legend.position="bottom",
  )+
  labs(colour="")
p2
ggsave(p2, file="./AMG_analysis/Viral_host_association/Viral_host_association_table_2.pdf", width = 10, height = 5, units = "cm")

# Step 3 Visualize the trend for table 3 (host is o__Kapabacteriales;f__UBA961)
## Read table
table3 <-read.table("./AMG_analysis/Viral_host_association/Viral_host_association_table_3.txt",sep = "\t",head=F, row.names = 1)
table3 <- t(table3)
table3 <- as.data.frame(table3)
# Change column names
colnames(table3) <- c("Month","Host.family","K01666.mhpE")

# Give the range of two y-axes
scaleFactor <- max(table3$Host.family) / max(table3$K01666.mhpE)
scaleFactor <- scaleFactor * 0.8

p3 <- ggplot(table3,aes(x=Month)) + 
  geom_line(aes(y = Host.family,colour="Host.family"), size = 0.5)+
  geom_line(aes(y = K01666.mhpE*scaleFactor,colour="K01666.mhpE"), size = 0.5)+
  scale_y_continuous(name="Host abundance", sec.axis=sec_axis(~./scaleFactor, name="Viral abundance")) + 
  scale_x_continuous(breaks=seq(1,12,1))+
  theme(
    legend.position="bottom",
  )+
  labs(colour="")
p3
ggsave(p3, file="./AMG_analysis/Viral_host_association/Viral_host_association_table_3.pdf", width = 10, height = 5, units = "cm")

# Step 4 Visualize the trend for table 4 (host is o__PCC-6307;f__Cyanobiaceae)
## Read table
table4 <-read.table("./AMG_analysis/Viral_host_association/Viral_host_association_table_4.txt",sep = "\t",head=F, row.names = 1)
table4 <- t(table4)
table4 <- as.data.frame(table4)
# Change column names
colnames(table4) <- c("Month","Host.family","K02703.psbA","K02706.psbD")

# Give the range of two y-axes
scaleFactor <- max(table4$Host.family) / max(table4$K02703.psbA)
scaleFactor <- scaleFactor * 0.8

p4 <- ggplot(table4,aes(x=Month)) + 
  geom_line(aes(y = Host.family,colour="Host.family"), size = 0.5)+
  geom_line(aes(y = K02703.psbA*scaleFactor,colour="K02703.psbA"), size = 0.5)+
  geom_line(aes(y = K02706.psbD*scaleFactor,colour="K02706.psbD"), size = 0.5)+
  scale_y_continuous(name="Host abundance", sec.axis=sec_axis(~./scaleFactor, name="Viral abundance")) + 
  scale_x_continuous(breaks=seq(1,12,1))+
  theme(
    legend.position="bottom",
  )+
  labs(colour="")
p4
ggsave(p4, file="./AMG_analysis/Viral_host_association/Viral_host_association_table_4.pdf", width = 10, height = 5, units = "cm")
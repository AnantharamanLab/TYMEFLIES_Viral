library(tidyverse)
library(reshape2)
library(ggpubr)

# Read VIBRANT result table
VIBRANT_result <- read.table ("../VIBRANT_result_summary.txt", sep = "\t", head = T, na.strings = c("NA"))

# Read read number and base number
Read_num_n_base_num <- read.table ("../Read_num_n_base_num.txt", sep = "\t", head = F, na.strings = c("NA"))

# Combined two tables
VIBRANT_result.2 <- left_join(VIBRANT_result, Read_num_n_base_num, by =c("IMG.ID" = "V1"))
colnames(VIBRANT_result.2)[9] <- "Read.num"
colnames(VIBRANT_result.2)[10] <- "Base.num"

head(VIBRANT_result.2)

# Normalize all parameters by 100M reads / metagenome
VIBRANT_result.2  <- VIBRANT_result.2 %>% mutate(Total.scaffold.num.normlized = Total.scaffold.num / (Read.num /100000000))
VIBRANT_result.2  <- VIBRANT_result.2 %>% mutate(Scaffold.num.over.min2000.normlized = Scaffold.num.over.min2000 / (Read.num /100000000))
VIBRANT_result.2  <- VIBRANT_result.2 %>% mutate(Phage.num.in.total.normlized = Phage.num.in.total / (Read.num /100000000))
VIBRANT_result.2  <- VIBRANT_result.2 %>% mutate(Lytic.phage.num.normlized = Lytic.phage.num / (Read.num /100000000))
VIBRANT_result.2  <- VIBRANT_result.2 %>% mutate(Lysogenic.phage..excluding.prophage..num.normlized = Lysogenic.phage..excluding.prophage..num / (Read.num /100000000))
VIBRANT_result.2  <- VIBRANT_result.2 %>% mutate(Prophage.num.normlized = Prophage.num / (Read.num /100000000))
VIBRANT_result.2  <- VIBRANT_result.2 %>% mutate(Complete.phage.num.normlized = Complete.phage.num / (Read.num /100000000))

# Plot boxplot
## Reshape the table and draw the bar plot
VIBRANT_result.3 <- VIBRANT_result.2[,c(1,11)] # Only use the 1st column
VIBRANT_result.4 <- melt(VIBRANT_result.3,
                         id.vars = c('IMG.ID'),
                         variable.name = "Scaffold.type",
                         value.name = 'Number.of.scaffolds')


p.part1 <- ggplot(data=VIBRANT_result.4,aes(x=Scaffold.type,y=Number.of.scaffolds))+
  geom_boxplot(aes(fill=Scaffold.type),outlier.shape = NA)+
  labs(title="VIBRANT result summary part 1",x="Scaffold assignment", y = "Scaffold number (Normalized by 100M reads/metagenome)")+
  #geom_jitter(width = 0.25, alpha = 0.1, color = 'black')+
  stat_summary(fun=mean, colour="darkred", geom="point", hape=18, size=3,show.legend = FALSE)+
  stat_summary(fun=mean, colour="red", geom="text", show.legend = FALSE, 
               vjust=-1, aes( label=round(..y.., digits=1)))+
  scale_y_continuous(labels = scales::scientific)+
  scale_x_discrete(labels=c("Total scaffold number"))+
  scale_fill_discrete(name = "Scaffold assignment", labels = c("Total scaffold number"))

p.part1.stat <- ggplot_build(p.part1)$data  
p.part1.stat <- as.data.frame(p.part1.stat)
rownames(p.part1.stat)[1] <- "Total scaffold number"


VIBRANT_result.3 <- VIBRANT_result.2[,c(1,12)] # Only use the later 5 columns
VIBRANT_result.4 <- melt(VIBRANT_result.3,
                         id.vars = c('IMG.ID'),
                         variable.name = "Scaffold.type",
                         value.name = 'Number.of.scaffolds')


p.part2 <- ggplot(data=VIBRANT_result.4,aes(x=Scaffold.type,y=Number.of.scaffolds))+
  geom_boxplot(aes(fill=Scaffold.type),outlier.shape = NA)+
  labs(title="VIBRANT result summary part 2",x="Scaffold assignment", y = "Scaffold number (Normalized by 100M reads/metagenome)")+
  #geom_jitter(width = 0.25, alpha = 0.1, color = 'black')+
  stat_summary(fun=mean, colour="darkred", geom="point", hape=18, size=3,show.legend = FALSE)+
  stat_summary(fun=mean, colour="red", geom="text", show.legend = FALSE, 
               vjust=-1, aes( label=round(..y.., digits=1)))+
  scale_y_continuous(labels = scales::scientific)+
  scale_x_discrete(labels=c("Scaffold num over 2000 bp"))+
  scale_fill_discrete(name = "Scaffold assignment", labels = c("Scaffold num over 2000 bp"))

p.part2.stat <- ggplot_build(p.part2)$data  
p.part2.stat <- as.data.frame(p.part2.stat)
rownames(p.part2.stat)[1] <- "Scaffold num over 2000 bp"

VIBRANT_result.3 <- VIBRANT_result.2[,c(1,13:17)] # Only use the later 5 columns
VIBRANT_result.4 <- melt(VIBRANT_result.3,
                         id.vars = c('IMG.ID'),
                         variable.name = "Scaffold.type",
                         value.name = 'Number.of.scaffolds')


p.part3 <- ggplot(data=VIBRANT_result.4,aes(x=Scaffold.type,y=Number.of.scaffolds))+
     geom_boxplot(aes(fill=Scaffold.type),outlier.shape = NA)+
     labs(title="VIBRANT result summary part 3",x="Scaffold assignment", y = "Scaffold number (Normalized by 100M reads/metagenome)")+
     #geom_jitter(width = 0.25, alpha = 0.1, color = 'black')+
     stat_summary(fun=mean, colour="darkred", geom="point", hape=18, size=3,show.legend = FALSE)+
     stat_summary(fun=mean, colour="red", geom="text", show.legend = FALSE, 
               vjust=-1, aes( label=round(..y.., digits=1)))+
     scale_y_continuous(labels = scales::scientific)+   
     scale_x_discrete(labels=c("Phage num in total", "Lytic phage num",	"Lysogenic phage (excluding prophage) num",	"Prophage num",	"Complete phage num"))+
     scale_fill_discrete(name = "Scaffold assignment", labels = c("Phage num in total", "Lytic phage num",	"Lysogenic phage (excluding prophage) num",	"Prophage num",	"Complete phage num"))

p.part3.stat <- ggplot_build(p.part3)$data  
p.part3.stat <- as.data.frame(p.part3.stat)
rownames(p.part3.stat)[1] <- "Phage num in total"
rownames(p.part3.stat)[2] <- "Lytic phage num"
rownames(p.part3.stat)[3] <- "Lysogenic phage (excluding prophage) num"
rownames(p.part3.stat)[4] <- "Prophage num"
rownames(p.part3.stat)[5] <- "Complete phage num"

ggarrange(p.part1,p.part2,p.part3,labels=c("A","B","C"), ncol=3, nrow=1, widths = c(1,1.1,3))

ggsave(file=paste0("VIBRANT_result_summary.pdf"), width = 60, height = 24, units = "cm")
ggsave(file=paste0("VIBRANT_result_summary.png"), width = 60, height = 24, units = "cm")


p.stat <- rbind(p.part1.stat, p.part2.stat, p.part3.stat)
p.stat <- subset(p.stat, select = -outliers)
p.stat.2 <- p.stat[,c(2:6,42)]
colnames(p.stat.2)[6] <- "mean"
write.table(p.stat.2,"VIBRANT_result_summary.tsv",na = "NA",sep = "\t",quote=F)

data<-read.csv("spearman-corr-heatmap-data.for_Cyanobacteria.csv")

library(tidyverse)
library(RColorBrewer)

head(data)

str(data$Environment)

data$Environment <- as.factor(data$Environment)
my_order <- c("Temperature","Dissolved oxygen","Secchi depth","Chlorophyll a","Chloride","Sulfate","Calcium","Magnesium","Sodium","Potassium","Iron","Manganese","pH","Dissolved inorganic carbon","Total inorganic carbon","Dissolved organic carbon","Total organic carbon","Nitrate plus nitrite","Ammonium","Total phosphorus","Soluble reactive phosphorus")

length(my_order)

str(data$Environment)

data$Environment <- factor(data$Environment, levels=my_order)
str(data$Environment)

colnames(data)[1] <- "Biology.Entity"

p <- ggplot(data, aes(x=Biology.Entity, y=Environment))+
  geom_tile(aes(fill=Correlation.Coefficient, col=is_significant), size=1)+
  scale_fill_gradient2(low = "#349ed3", 
                       midpoint = 0, 
                       mid = "white", 
                       high = "#ff5b4d")+
  scale_color_manual(values=c("black","white"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=0, vjust=0),
        panel.grid = element_blank(),
        panel.border = element_blank())+
  
  scale_x_discrete(position = "top")+
  xlab("Biological entity")+
  ylab("Environmental variable")+
  scale_y_discrete(limits=rev)+
  geom_text(data=data,aes(x=Biology.Entity, y=Environment,label=Correlation.Coefficient))

p

ggsave("spearman-corr-heatmap-data.for_Cyanobacteria.pdf", plot = p, device = "pdf", width = 6, height = 8)

library(tidyverse)
graph1.data  <- read.csv("Cyanobacteria_group2host_abun.csv")
colnames(graph1.data)[1] <- "Days_since_early_summer_start"
graph1.data <- graph1.data %>% mutate(sum = Cyanobiaceae+Microcystis+Planktothrix+Other_Cyanobacteria)

for (i in 2:5){
  for (j in 1:nrow(graph1.data)){
    print(j)
    print(i)
    graph1.data[j,i] <- graph1.data[j,i]/graph1.data$sum[j]
  }
}


graph1.data <- graph1.data %>% gather("taxa","rel_abund",2:5)

host.colors <- c(Cyanobiaceae="#0099CC",
                      Microcystis="#336666",
                      Planktothrix="#FF6666",
                      Other_Cyanobacteria="#d3d3d3")
graph1.data$taxa <- factor(graph1.data$taxa, 
                           levels = c("Other_Cyanobacteria",
                                      "Planktothrix", 
                                      "Microcystis",
                                      "Cyanobiaceae"))

ggplot(graph1.data,aes(x=Days_since_early_summer_start, 
                       y=rel_abund, 
                       fill=taxa))+
         geom_area()+
  scale_fill_manual(values=host.colors)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom")+
  ylab("Relative abundance")+
  xlab("Days since Early Summer start")+
  geom_vline(xintercept=0)+
  labs(fill='Host taxonomy')+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks = seq(from = -50, to = 200, by = 25))

ggsave("Cyanobacteria_group2host_abun.pdf", width = 12, height = 4, units = "cm")


#### GRAPH 2 ####
graph2.data  <- read.csv("Cyanobacteria_group2virus_abun.csv")
colnames(graph2.data)[1] <- "Days_since_early_summer_start"
graph2.data <- graph2.data %>% mutate(sum = Cyanobiaceae+Microcystis+Planktothrix+Other_Cyanobacteria)

for (i in 2:5){
  for (j in 1:nrow(graph2.data)){
    print(j)
    print(i)
    graph2.data[j,i] <- graph2.data[j,i]/graph2.data$sum[j]
  }
}


graph2.data <- graph2.data %>% gather("taxa","rel_abund",2:5)

host.colors <- c(Cyanobiaceae="#0099CC",
                 Microcystis="#336666",
                 Planktothrix="#FF6666",
                 Other_Cyanobacteria="#d3d3d3")


graph2.data$taxa <- factor(graph2.data$taxa, 
                           levels = c("Other_Cyanobacteria",
                                      "Planktothrix", 
                                      "Microcystis",
                                      "Cyanobiaceae"))

ggplot(graph2.data,aes(x=Days_since_early_summer_start, 
                       y=rel_abund, 
                       fill=taxa))+
  geom_area()+
  scale_fill_manual(values=host.colors)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom")+
  ylab("Relative abundance")+
  xlab("Days since Early Summer start")+
  geom_vline(xintercept=0)+
  labs(fill='Host taxonomy')+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks = seq(from = -50, to = 200, by = 25))

ggsave("Cyanobacteria_group2host_virus.pdf", width = 12, height = 4, units = "cm")


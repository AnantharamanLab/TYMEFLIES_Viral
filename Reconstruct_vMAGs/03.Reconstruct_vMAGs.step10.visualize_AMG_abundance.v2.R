library(tidyverse)
library(ggpubr)

# Note:
# Before running the datasets, in the input files,
# all the seasons have been converted to numbers as:
# "Spring" => 1
# "Clearwater" => 2
# "Early Summer" => 3
# "Late Summer" => 4
# "Fall" => 5
# "Ice-on" => 6

# Step 1. Make AMG trend R plots for each season (written by Chao)
# Read table
table <-read.table("AMG_analysis/KO2season2abun.txt",sep = "\t",head=F, row.names = 1)
table.tbl <- as_tibble(table)
ko_details <-read.table("AMG_analysis/KO2month_ko_details.txt",sep = "\t",head=F)

# Print the AMG trend plots within a loop
for (i in 2:579){ # The total row number of table
  table.tbl.subset <- table.tbl %>% slice(1,i)
  table.tbl.subset <- t(table.tbl.subset)
  table.tbl.subset <- as.data.frame(table.tbl.subset)

  j <- i - 1
  ko_id <- ko_details[j,1]
  ko_detail <- ko_details[j,2]

  p <- ggplot() + geom_line(aes(y = V2, x = V1, color = "brown"), size = 1.5,
              data = table.tbl.subset, stat="identity")+
              theme(legend.position="bottom")+
              ggtitle(paste0(ko_id, " ", ko_detail)) + labs(x="Season", y="KO abundance")
  p1 <- p + scale_x_continuous(breaks=seq(1,6,1))
  p1
  ggsave(p1, file=paste0("AMG_trend_R_plots/PDFs_3/",ko_id,".pdf"), width = 20, height = 10, units = "cm")
  ggsave(p1, file=paste0("AMG_trend_R_plots/PNGs_3/",ko_id,".png"), width = 20, height = 10, units = "cm")
}

# Step 2. Make AMG trend R plots for each year_season (for example: 2000-Spring)
# 20 years will be plotted in 20 facets in a single figure
## Read KO abundance for each year each season 
table <-read.table("AMG_analysis/KO2year_season2abun.txt",sep = "\t",head=F, row.names = 1)
table.tbl <- as_tibble(table)

## Read KO details
ko_details <-read.table("AMG_analysis/KO2month_ko_details.txt",sep = "\t",head=F)

## Read and store metagenome number for each year each season
mg_per_season_per_year <- read.table("AMG_analysis/Year_season2num_of_metagenomes.txt", sep="\t")
mg_per_season_per_year <- mg_per_season_per_year %>% separate(V1, sep="-", into=c("Year","season"), remove=FALSE)
colnames(mg_per_season_per_year)[4] <- "nb_mg"

# Print the AMG trend plots within a loop
for (i in 2:579){ # The total row number of table (total number AMG KOs)
  #i <- 2
  table.tbl.subset <- table.tbl %>% slice(1,i) # The subset of the row i 
  
    # Select table.tbl.subset by column
    table.tbl.subset <- t(table.tbl.subset)
    table.tbl.subset <- as.data.frame(table.tbl.subset)
    str(table.tbl.subset)
    
    table.tbl.subset <- table.tbl.subset %>% separate(V1, sep="-", into=c("Year","season"),remove=FALSE)
    table.tbl.subset$V2 <- as.numeric(table.tbl.subset$V2)

    # Store KO details
    j <- i - 1
    ko_id <- ko_details[j,1]
    print(ko_id)
    
    ko_detail <- ko_details[j,2]
    print(ko_detail)
    
    title_of_plot <- paste(ko_id, ko_detail)
    print(title_of_plot)
    
    # Change the season column into integer format
    table.tbl.subset$season <- as.integer(table.tbl.subset$season)
    
    str(table.tbl.subset)
    str(mg_per_season_per_year)
    
    mg_per_season_per_year$season <- as.integer(mg_per_season_per_year$season)
    
    # Combine two tables 
    table.tbl.subset <- left_join(table.tbl.subset, mg_per_season_per_year)
    head(table.tbl.subset)
    colnames(table.tbl.subset)[4] <- "KO"
    
    head(table.tbl.subset)
    
    # Plot the Panel A heatmap:
    plot1 <- ggplot(table.tbl.subset, aes(x=season, y=Year, fill=nb_mg))+
             geom_tile(color="black")+
             theme_bw()+
             scale_fill_continuous(low="white",high="red",name="Metagenome number")+
             geom_text(aes(label = nb_mg), size=2.5)+
             theme(text=element_text(size=11))+
             scale_x_continuous(breaks=seq(1,6,1))+
             labs(x="Season", y="Year")

    # Plot Panel B
    table.tbl.subset  <- table.tbl.subset %>% group_by(season) %>% mutate(KO_valid = case_when(nb_mg != 0 ~ KO))
  
    plot2 <- ggplot() + geom_line(aes(y = KO_valid, x = season), color="brown", size = 1.5,
                                  data = table.tbl.subset)+
             geom_point(aes(y = KO_valid, x = season), color="brown", size = 1.5,
                        data = table.tbl.subset)+
             theme(legend.position="bottom", 
                   text=element_text(size=12))+
             ggtitle(paste0(ko_id, " ", ko_detail)) + labs(x="Season", y="KO abundance")+
             facet_wrap(Year~., scales = "free", nrow = 10, ncol=2)+ 
             scale_x_continuous(breaks=seq(1,6,1))
    
    # Plot Panel C
    MG2season <- mg_per_season_per_year %>% group_by(season) %>% tally(nb_mg)
    
    table.tbl.subset.C <-left_join(table.tbl.subset, MG2season, by="season")
    
    table.tbl.subset.C <- table.tbl.subset.C %>% group_by(season) %>% mutate(KO_by_season=sum(KO * nb_mg)/n)

    plot3 <- ggplot(table.tbl.subset.C, aes(x=season, y=KO_by_season))+
             geom_line(size=1.5)+ geom_point(size=1.5)+
             ggtitle(paste0(ko_id, " ", ko_detail), 
             subtitle="Average by season across all years")+
             labs(x="Season", y="KO abundance") + 
             scale_x_continuous(breaks=seq(1,6,1))+
             theme(text=element_text(size=12))

    ggarrange(plot1,plot2,plot3,labels=c("A","B","C"), ncol=1, nrow=3, heights = c(0.2,1,0.2))
    
    ggsave(file=paste0("AMG_trend_R_plots/PDFs_4/",ko_id,".pdf"), width = 8.5, height = 20, units = "in")
    ggsave(file=paste0("AMG_trend_R_plots/PNGs_4/",ko_id,".png"), width = 8.5, height = 20, units = "in")

    print(i)
}


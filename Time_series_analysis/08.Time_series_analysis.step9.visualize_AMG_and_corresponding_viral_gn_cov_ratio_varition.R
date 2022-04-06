library(tidyverse)
library(ggpubr)

# Step 1 Draw AMG gene cov ratio plots
# 20 years will be plotted in 20 facets in a single figure
## Read AMG gene cov ratio for each year each month 
table <-read.table("MetaPop/AMG_gene2year_month2cov_ratio.txt",sep = "\t",head=F, row.names = 1)
table.tbl <- as_tibble(table)

## Read KO details
amg_gene2ko_n_detail <- read.table("MetaPop/AMG_gene2ko_n_detail.txt",sep = "\t",head=F)

## Read and store metagenome number for each year each month
mg_per_month_per_year <- read.table("MetaPop/Year_month2num_of_metagenomes.txt", sep="\t")
mg_per_month_per_year <- mg_per_month_per_year %>% separate(V1, sep="-", into=c("Year","month"), remove=FALSE)
colnames(mg_per_month_per_year)[4] <- "nb_mg"

# Print the AMG gene cov ratio plots within a loop
for (i in 2:122){ # The total row number of table (total number AMG gene)
  
  table.tbl.subset <- table.tbl %>% slice(1,i) # The subset of the row i 
  
    # Select table.tbl.subset by column
    table.tbl.subset <- t(table.tbl.subset)
    table.tbl.subset <- as.data.frame(table.tbl.subset)
    str(table.tbl.subset)
    
    table.tbl.subset <- table.tbl.subset %>% separate(V1, sep="-", into=c("Year","month"),remove=FALSE)
    table.tbl.subset$V2 <- as.numeric(table.tbl.subset$V2)

    # Store KO details
    j <- i - 1
    amg_gene <- amg_gene2ko_n_detail[j,1]
    ko_n_detail <- amg_gene2ko_n_detail[j,2]
    
    # Change the month column into integer format
    table.tbl.subset$month <- as.integer(table.tbl.subset$month)
    
    str(table.tbl.subset)
    str(mg_per_month_per_year)
    
    mg_per_month_per_year$month <- as.integer(mg_per_month_per_year$month)
    
    # Combine two tables 
    table.tbl.subset <- left_join(table.tbl.subset, mg_per_month_per_year)
    head(table.tbl.subset)
    colnames(table.tbl.subset)[4] <- "KO"
    
    head(table.tbl.subset)
    
    # Plot the Panel A heatmap:
    plot1 <- ggplot(table.tbl.subset, aes(x=month, y=Year, fill=nb_mg))+
             geom_tile(color="black")+
             theme_bw()+
             scale_fill_continuous(low="white",high="red",name="Metagenome number")+
             geom_text(aes(label = nb_mg), size=2.5)+
             theme(text=element_text(size=11))+
             scale_x_continuous(breaks=seq(1,12,1))+
             labs(x="Month", y="Year")

    # Plot Panel B
    table.tbl.subset  <- table.tbl.subset %>% group_by(month) %>% mutate(KO_valid = case_when(nb_mg != 0 ~ KO))
  
    plot2 <- ggplot() + geom_line(aes(y = KO_valid, x = month), color="brown", size = 1.5,
                                  data = table.tbl.subset)+
             geom_point(aes(y = KO_valid, x = month), color="brown", size = 1.5,
                        data = table.tbl.subset)+
             theme(legend.position="bottom", 
                   text=element_text(size=12))+
             ggtitle(amg_gene, subtitle=ko_n_detail) + 
             labs(x="Month", y="AMG coverage ratio")+
             facet_wrap(Year~., scales = "free", nrow = 10, ncol=2)+ 
             scale_x_continuous(breaks=seq(1,12,1))
    
    # Plot Panel C
    MG2month <- mg_per_month_per_year %>% group_by(month) %>% tally(nb_mg)
    
    table.tbl.subset.C <-left_join(table.tbl.subset, MG2month, by="month")
    
    table.tbl.subset.C <- table.tbl.subset.C %>% group_by(month) %>% mutate(KO_by_month=sum(KO * nb_mg)/n)

    plot3 <- ggplot(table.tbl.subset.C, aes(x=month, y=KO_by_month))+
             geom_line(size=1.5)+ geom_point(size=1.5)+
             ggtitle(amg_gene, subtitle=ko_n_detail)+
             labs(x="Month", y="AMG coverage ratio") + 
             scale_x_continuous(breaks=seq(1,12,1))+
             theme(text=element_text(size=12))

    ggarrange(plot1,plot2,plot3,labels=c("A","B","C"), ncol=1, nrow=3, heights = c(0.2,1,0.2))
    
    ggsave(file=paste0("AMG_cov_ratio_R_plots/PDFs/",j,".pdf"), width = 8.5, height = 20, units = "in")
    ggsave(file=paste0("AMG_cov_ratio_R_plots/PNGs/",j,".png"), width = 8.5, height = 20, units = "in")

    print(i)
}


# Step 2 Draw AMG_gene_containing_viral_gn coverage plots
# Note: Before running, clean all the files in Environment and Console
# 20 years will be plotted in 20 facets in a single figure
## Read AMG gene cov ratio for each year each month 
table <-read.table("MetaPop/AMG_gene_containing_viral_gn2year_month2cov_ratio.txt",sep = "\t",head=F, row.names = 1)
table.tbl <- as_tibble(table)

## Read KO details
AMG_gene_containing_viral_gn2kos <- read.table("MetaPop/AMG_gene_containing_viral_gn2KOs.txt",sep = "\t",head=F)

## Read and store metagenome number for each year each month
mg_per_month_per_year <- read.table("MetaPop/Year_month2num_of_metagenomes.txt", sep="\t")
mg_per_month_per_year <- mg_per_month_per_year %>% separate(V1, sep="-", into=c("Year","month"), remove=FALSE)
colnames(mg_per_month_per_year)[4] <- "nb_mg"

# Print the AMG_gene_containing_viral_gn cov plots within a loop
for (i in 2:108){ # The total row number of table (total number AMG_gene_containing_viral_gn)
  
  table.tbl.subset <- table.tbl %>% slice(1,i) # The subset of the row i 
  
  # Select table.tbl.subset by column
  table.tbl.subset <- t(table.tbl.subset)
  table.tbl.subset <- as.data.frame(table.tbl.subset)
  str(table.tbl.subset)
  
  table.tbl.subset <- table.tbl.subset %>% separate(V1, sep="-", into=c("Year","month"),remove=FALSE)
  table.tbl.subset$V2 <- as.numeric(table.tbl.subset$V2)
  
  # Store KO details
  j <- i - 1
  AMG_gene_containing_viral_gn <- AMG_gene_containing_viral_gn2kos[j,1]
  kos <- AMG_gene_containing_viral_gn2kos[j,2]
  
  # Change the month column into integer format
  table.tbl.subset$month <- as.integer(table.tbl.subset$month)
  
  str(table.tbl.subset)
  str(mg_per_month_per_year)
  
  mg_per_month_per_year$month <- as.integer(mg_per_month_per_year$month)
  
  # Combine two tables 
  table.tbl.subset <- left_join(table.tbl.subset, mg_per_month_per_year)
  head(table.tbl.subset)
  colnames(table.tbl.subset)[4] <- "KO"
  
  head(table.tbl.subset)
  
  # Plot the Panel A heatmap:
  plot1 <- ggplot(table.tbl.subset, aes(x=month, y=Year, fill=nb_mg))+
    geom_tile(color="black")+
    theme_bw()+
    scale_fill_continuous(low="white",high="red",name="Metagenome number")+
    geom_text(aes(label = nb_mg), size=2.5)+
    theme(text=element_text(size=11))+
    scale_x_continuous(breaks=seq(1,12,1))+
    labs(x="Month", y="Year")
  
  # Plot Panel B
  table.tbl.subset  <- table.tbl.subset %>% group_by(month) %>% mutate(KO_valid = case_when(nb_mg != 0 ~ KO))
  
  plot2 <- ggplot() + geom_line(aes(y = KO_valid, x = month), color="brown", size = 1.5,
                                data = table.tbl.subset)+
    geom_point(aes(y = KO_valid, x = month), color="brown", size = 1.5,
               data = table.tbl.subset)+
    theme(legend.position="bottom", 
          text=element_text(size=12))+
    ggtitle(AMG_gene_containing_viral_gn, subtitle=kos) + 
    labs(x="Month", y="AMG gene containing viral genome coverage")+
    facet_wrap(Year~., scales = "free", nrow = 10, ncol=2)+ 
    scale_x_continuous(breaks=seq(1,12,1))
  
  # Plot Panel C
  MG2month <- mg_per_month_per_year %>% group_by(month) %>% tally(nb_mg)
  
  table.tbl.subset.C <-left_join(table.tbl.subset, MG2month, by="month")
  
  table.tbl.subset.C <- table.tbl.subset.C %>% group_by(month) %>% mutate(KO_by_month=sum(KO * nb_mg)/n)
  
  plot3 <- ggplot(table.tbl.subset.C, aes(x=month, y=KO_by_month))+
    geom_line(size=1.5)+ geom_point(size=1.5)+
    ggtitle(AMG_gene_containing_viral_gn, subtitle=kos)+
    labs(x="Month", y="AMG gene containing viral genome coverage") + 
    scale_x_continuous(breaks=seq(1,12,1))+
    theme(text=element_text(size=12))
  
  ggarrange(plot1,plot2,plot3,labels=c("A","B","C"), ncol=1, nrow=3, heights = c(0.2,1,0.2))
  
  ggsave(file=paste0("AMG_gene_containing_viral_gn_cov_R_plots/PDFs/",j,".pdf"), width = 8.5, height = 20, units = "in")
  ggsave(file=paste0("AMG_gene_containing_viral_gn_cov_R_plots/PNGs/",j,".png"), width = 8.5, height = 20, units = "in")
  
  print(i)
}
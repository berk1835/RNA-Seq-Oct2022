
## Creating a biotype pie chart for DEGs
## 
## 
##
## 
##
## Author: Rebekah White Dec 2022, Weadick lab, University of 
## Exeter



## Set wd

setwd("C:/Users/rw617/OneDrive - University of Exeter/rw617/08_Transcriptomics")


## Load in packages

#options("install.lock"=FALSE) #prevents Failed To Lock Directory /00LOCK error
#install.packages("qdapTools")

#library(qdapTools)



## EDIT: State study

e <- "dr DEGs biotype.png"



## EDIT: Load in data

dat <-  read.csv("./01_PpaDr/degs/tidy_dr_degs.csv")





## Seperate upreg from downreg

up <-  subset(dat, log2FoldChange > 0,
              select=c(locus,	baseMean,	log2FoldChange,	pvalue,	padj,	pfam_id.y,	
                       Description,	WBID))


down <-  subset(dat, log2FoldChange < 0,
                select=c(locus,	baseMean,	log2FoldChange,	pvalue,	padj,	pfam_id.y,	
                         Description,	WBID))






## Set biotype labels

uplabels <- up$Description

downlabels <- down$Description







# Convert labels list into data frame of counts

upcounts <- table(uplabels)
downcounts <- table(downlabels)

upcounts <- data.frame(upcounts)
downcounts <- data.frame(downcounts)






# Create pie charts of everything 

par(mfrow = c(1, 2), oma=c(-2,-2,-2,-2))
pie(upcounts$Freq, upcounts$uplabels, cex=0.5, main="Upregulated in DR")
pie(downcounts$Freq, downcounts$downlabels, cex=0.5, main="Downregulated in DR")







# Seperate biotypes where frequency is less than 2

upcless <-  subset(upcounts, Freq < 2,
                  select=c(uplabels, Freq))

upcmore <-  subset(upcounts, Freq > 2,
               select=c(uplabels, Freq))



downcless <-  subset(downcounts, Freq < 2,
                   select=c(uplabels, Freq))

downcmore <-  subset(downcounts, Freq > 2,
                   select=c(uplabels, Freq))







# NOT UPDATED # Count number of freq = 1 and create column 'other'

other_freq <- sum(cless$Freq)

new_row <- list("other", other_freq)
new_row <- data.frame(new_row)

colnames(new_row)[1] ="labels"
colnames(new_row)[2] ="Freq"

tidyc <- rbind(cmore, new_row)





# Create pie chart (for checking)

pie(tidyc$Freq, tidyc$labels, cex=0.5)

#par(mfrow = c(1, 2))
#pie(upcounts$Freq, upcounts$uplabels, cex=0.5, main="Upregulated in ya")
#pie(downcounts$Freq, downcounts$downlabels, cex=0.5, main="Downregulated in ya")


  
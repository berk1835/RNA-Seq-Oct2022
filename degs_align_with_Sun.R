## Align DEGs from DESeq2 with supplementary material table 8 
## from Sun et al. 2021
##
## This will give locus functional descriptions that are not 
## available elsewhere (i.e., SimpleMine for WB_ID)
##
## For new data, edit start and end sections. 
##
## Author: Rebekah White Dec 2022, Weadick lab, University of 
## Exeter



## Set wd

setwd("C:/Users/rw617/OneDrive - University of Exeter/rw617/08_Transcriptomics")


## Load in Sun et al. database

sun <- read.csv("sun_et_al_tableS8_useful.csv")



# rename ppa_locus to locus

colnames(sun)[2] ="locus"



# EDIT: Load in the DEGs from DESeq2 (or alternative)
# name the first column locus

deg <- read.csv("./02_PpaAges/degs/ages_mid-vs-old_padj01_log2fold.csv") 

colnames(deg)[1] ="locus"


## Subset loci in common

loc <- sun[sun$locus %in% deg$locus, ]






#### Add information from Sun back into DEGs file

# merge loc and degs, and remove unnecessary column

new_degs <- merge(deg, loc, by = 'locus')
new_degs <- subset(new_degs, select = -c(evalue_full_seq))




## Create list of unmatched loci

non <- deg[!(deg$locus %in% sun$locus), ]




# Add same columns as new_degs with value n/a, and rename differing ones

non['Description'] <- NA
non['pfam_id'] <- NA





# Add unmatched loci to end of matched file 

new_degs2 <- rbind(new_degs, non)





# Check merge and tidy up if necessary 

#View(new_degs2)
new_degs_tidy = subset(new_degs2, select = c(locus, baseMean, log2FoldChange, pvalue, padj,
                                        pfam_id, Description))





# EDIT: Save to new file

write.csv(new_degs_tidy, "./02_PpaAges/degs/tidy_ages-mid_old_degs.csv", row.names=FALSE)








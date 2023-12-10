library(readxl)
library(tidyverse)
library(ape)
library(phytools)
library(geiger)

stpt_chr <- read_xlsx("/Users/hu_zhehao/Downloads/Charactersheet ggf. für ancestral state analyse.xlsx", range = "A2:I31")
stpt_tree <- read.nexus("/Users/hu_zhehao/Library/Mobile Documents/com~apple~CloudDocs/UHH/M. Sc. Biologie/EvoSys/Exercise/Project/ASR/IQtree.nex")
# HOW TO IMPORT TREES -> https://cran.r-project.org/web/packages/TreeTools/vignettes/load-trees.html
plot(stpt_tree)

row.names(stpt_chr) <- stpt_chr$Characters

# subset
aqua <- stpt_chr$Aquatisch
names(aqua) <- row.names(stpt_chr)
str(aqua)

# Discret data 
# isolate the binary character
#find the best fitting model
fitDiscrete(stpt_tree, aqua, model = "ER")
# AICc = 13.722218, lnL = -5.761109
fitDiscrete(stpt_tree, aqua, model = "SYM")
# AICc = 13.722218, lnL = -5.761109
fitDiscrete(stpt_tree, aqua, model = "ARD")
# AICc = 15.911984, lnL = -5.640202
fitDiscrete(stpt_tree, aqua, model = "meristic")
# AICc = 13.522218, lnL = -5.761109
# If tie, choose the model with least parameters
# meristic has lowest AICc

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(aqua, stpt_tree, type = "discret", model = "ER")
# Error while using meristic

# Plot the reconstructed tree
cols <- c("blue", "red")
plotTree(stpt_tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:stpt_tree$Nnode+Ntip(stpt_tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(aqua[stpt_tree$tip.label],levels(aqua)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)

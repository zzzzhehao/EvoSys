library(readxl)
library(tidyverse)
library(ape)
library(phytools)
library(geiger)

#### BEAST TREE ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
stpt_chr <- read_xlsx("Charactersheet ggf. für ancestral state analyse.xlsx", range = "A2:I31") %>% mutate(Lebenszyklus = case_when(
  Lebenszyklus == "1" ~ "haplont",
  Lebenszyklus == "2" ~ "diplont",
  Lebenszyklus == "3" ~ "haplodiplont"
))
stpt_tree <- read.nexus("beast4 NEXtree")
# HOW TO IMPORT TREES -> https://cran.r-project.org/web/packages/TreeTools/vignettes/load-trees.html
plot(stpt_tree)

row.names(stpt_chr) <- stpt_chr$Characters

# subset
lifecycle <- as.factor(stpt_chr$Lebenszyklus)
names(lifecycle) <- row.names(stpt_chr)
str(lifecycle)

# Discret data 
# isolate the binary character
#find the best fitting model
fitDiscrete(stpt_tree, lifecycle, model = "ER")
# AICc = 55.882959, fp = 1
fitDiscrete(stpt_tree, lifecycle, model = "SYM")
# AICc = 52.972050, fp = 3
fitDiscrete(stpt_tree, lifecycle, model = "ARD")
# AICc = 60.954145, fp = 6
fitDiscrete(stpt_tree, lifecycle, model = "meristic")
# AICc = 64.524249, fp = 2
# If tie, choose the model with least parameters
# SYM has lowest AICc

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(lifecycle, stpt_tree, type = "discret", model = "ARD")


# Plot the reconstructed tree
cols <- c("blue", "red", "green")
plotTree(stpt_tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:stpt_tree$Nnode+Ntip(stpt_tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(lifecycle[stpt_tree$tip.label],levels(lifecycle)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)


#### ML TREE ####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
stpt_chr <- read_xlsx("Charactersheet ggf. für ancestral state analyse.xlsx", range = "A2:I31")
stpt_tree <- read.nexus("IQTREE NEXtree")
# HOW TO IMPORT TREES -> https://cran.r-project.org/web/packages/TreeTools/vignettes/load-trees.html
plot(stpt_tree)

row.names(stpt_chr) <- stpt_chr$Characters

# subset
lifecycle <- as.factor(stpt_chr$Lebenszyklus)
names(lifecycle) <- row.names(stpt_chr)
str(lifecycle)

# Discret data 
# isolate the binary character
#find the best fitting model
fitDiscrete(stpt_tree, lifecycle, model = "ER")
# AICc = 47.207749, fp = 1
fitDiscrete(stpt_tree, lifecycle, model = "SYM")
# AICc = 44.524491, fp = 3
fitDiscrete(stpt_tree, lifecycle, model = "ARD")
# AICc = 53.645877, fp = 6
fitDiscrete(stpt_tree, lifecycle, model = "meristic")
# AICc = 49.864008, fp = 2
# If tie, choose the model with least parameters
# SYM has lowest AICc

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(lifecycle, stpt_tree, type = "discret", model = "SYM")
# Error while using meristic

# Plot the reconstructed tree
cols <- c("blue", "red", "green")
plotTree(stpt_tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:stpt_tree$Nnode+Ntip(stpt_tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(lifecycle[stpt_tree$tip.label],levels(lifecycle)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)



library(readxl)
library(tidyverse)
library(ape)
library(phytools)
library(geiger)

#### Lifecycle ####

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




#### aqua ####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
stpt_chr <- read_xlsx("Charactersheet ggf. für ancestral state analyse.xlsx", range = "A2:I31")
stpt_tree <- read.nexus("IQTREE NEXtree")
# HOW TO IMPORT TREES -> https://cran.r-project.org/web/packages/TreeTools/vignettes/load-trees.html
plot(stpt_tree)

row.names(stpt_chr) <- stpt_chr$Characters

# subset
aqua <- as.factor(stpt_chr$Aquatisch)
names(aqua) <- row.names(stpt_chr)
str(aqua)

# Discret data 
# isolate the binary character
#find the best fitting model
fitDiscrete(stpt_tree, aqua, model = "ER")
# AICc = 14.056289, fp = 1
fitDiscrete(stpt_tree, aqua, model = "SYM")
# AICc = 14.056289, fp = 1
fitDiscrete(stpt_tree, aqua, model = "ARD")
# AICc = 15.984681, fp = 2
fitDiscrete(stpt_tree, aqua, model = "meristic")
# AICc = 14.056289, fp = 1
# If tie, choose the model with least parameters
# ER, SYM, meristic tie

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(aqua, stpt_tree, type = "discret", model = "ER")
# Error while using meristic

# Plot the reconstructed tree
cols <- c("blue", "red", "green")
plotTree(stpt_tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:stpt_tree$Nnode+Ntip(stpt_tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(aqua[stpt_tree$tip.label],levels(aqua)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)




#### vas ####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
stpt_chr <- read_xlsx("Charactersheet ggf. für ancestral state analyse.xlsx", range = "A2:I31")
stpt_tree <- read.nexus("IQTREE NEXtree")
# HOW TO IMPORT TREES -> https://cran.r-project.org/web/packages/TreeTools/vignettes/load-trees.html
plot(stpt_tree)

row.names(stpt_chr) <- stpt_chr$Characters

# subset
vas <- as.factor(stpt_chr$Vaskulärsystem)
names(vas) <- row.names(stpt_chr)
str(vas)

# Discret data 
# isolate the binary character
#find the best fitting model
fitDiscrete(stpt_tree, vas, model = "ER")
# AICc = 21.166001, fp = 1
fitDiscrete(stpt_tree, vas, model = "SYM")
# AICc = 21.166001, fp = 1
fitDiscrete(stpt_tree, vas, model = "ARD")
# AICc = 22.684097, fp = 2
fitDiscrete(stpt_tree, vas, model = "meristic")
# AICc = 21.166001, fp = 1
# If tie, choose the model with least parameters
# ER, SYM, meristic tie

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(vas, stpt_tree, type = "discret", model = "ER")
# Error while using meristic

# Plot the reconstructed tree
cols <- c("blue", "red", "green")
plotTree(stpt_tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:stpt_tree$Nnode+Ntip(stpt_tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(vas[stpt_tree$tip.label],levels(vas)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)



#### rhizoid ####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
stpt_chr <- read_xlsx("Charactersheet ggf. für ancestral state analyse.xlsx", range = "A2:I31")
stpt_tree <- read.nexus("IQTREE NEXtree")
# HOW TO IMPORT TREES -> https://cran.r-project.org/web/packages/TreeTools/vignettes/load-trees.html
plot(stpt_tree)

row.names(stpt_chr) <- stpt_chr$Characters

# subset
rz <- as.factor(stpt_chr$Rhizoide)
names(rz) <- row.names(stpt_chr)
str(rz)

# Discret data 
# isolate the binary character
#find the best fitting model
fitDiscrete(stpt_tree, rz, model = "ER")
# AICc = 15.414243, fp = 1
fitDiscrete(stpt_tree, rz, model = "SYM")
# AICc = 15.414243, fp = 1
fitDiscrete(stpt_tree, rz, model = "ARD")
# AICc = 17.604169, fp = 2
fitDiscrete(stpt_tree, rz, model = "meristic")
# AICc = 15.414243, fp = 1
# If tie, choose the model with least parameters
# ER, SYM, meristic tie

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(rz, stpt_tree, type = "discret", model = "ER")
# Error while using meristic

# Plot the reconstructed tree
cols <- c("blue", "red", "green")
plotTree(stpt_tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:stpt_tree$Nnode+Ntip(stpt_tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(rz[stpt_tree$tip.label],levels(rz)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)



#### gso ####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
stpt_chr <- read_xlsx("Charactersheet ggf. für ancestral state analyse.xlsx", range = "A2:I31")
stpt_tree <- read.nexus("IQTREE NEXtree")
# HOW TO IMPORT TREES -> https://cran.r-project.org/web/packages/TreeTools/vignettes/load-trees.html
plot(stpt_tree)

row.names(stpt_chr) <- stpt_chr$Characters

# subset
gso <- as.factor(stpt_chr$Geschlechtsorgane)
names(gso) <- row.names(stpt_chr)
str(gso)

# Discret data 
# isolate the binary character
#find the best fitting model
fitDiscrete(stpt_tree, gso, model = "ER")
# AICc = 40.544342, fp = 1
fitDiscrete(stpt_tree, gso, model = "SYM")
# AICc = 40.544342, fp = 1
fitDiscrete(stpt_tree, gso, model = "ARD")
# AICc = 42.822219, fp = 2
fitDiscrete(stpt_tree, gso, model = "meristic")
# AICc = 40.544342, fp = 1
# If tie, choose the model with least parameters
# ER, SYM, meristic tie

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(gso, stpt_tree, type = "discret", model = "ER")
# Error while using meristic

# Plot the reconstructed tree
cols <- c("blue", "red", "green")
plotTree(stpt_tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:stpt_tree$Nnode+Ntip(stpt_tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(gso[stpt_tree$tip.label],levels(gso)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)



#### fer ####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
stpt_chr <- read_xlsx("Charactersheet ggf. für ancestral state analyse.xlsx", range = "A2:I31")
stpt_tree <- read.nexus("IQTREE NEXtree")
# HOW TO IMPORT TREES -> https://cran.r-project.org/web/packages/TreeTools/vignettes/load-trees.html
plot(stpt_tree)

row.names(stpt_chr) <- stpt_chr$Characters

# subset
fer <- as.factor(stpt_chr$Befruchtungstyp)
names(fer) <- row.names(stpt_chr)
str(fer)

# Discret data 
# isolate the binary character
#find the best fitting model
fitDiscrete(stpt_tree, fer, model = "ER")
# AICc = 56.620408, fp = 1
fitDiscrete(stpt_tree, fer, model = "SYM")
# AICc = 80.695153, fp = 10
fitDiscrete(stpt_tree, fer, model = "ARD")
# AICc = 188.969667, fp = 20
fitDiscrete(stpt_tree, fer, model = "meristic")
# AICc = 61.295717, fp = 4
# If tie, choose the model with least parameters
# ER

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(fer, stpt_tree, type = "discret", model = "ER")
# Error while using meristic

# Plot the reconstructed tree
cols <- c("blue", "red", "green","yellow","cyan")
plotTree(stpt_tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:stpt_tree$Nnode+Ntip(stpt_tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(fer[stpt_tree$tip.label],levels(fer)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)



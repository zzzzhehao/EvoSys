library(ape)
library(phytools)
library(geiger)

#### Continuous data ####

data("mammal.data")
data("mammal.tree")
plot(mammal.tree)

# subset data
bodymass <- log(mammal.data$bodyMass)
names(bodymass) <- row.names(mammal.data)

homerange <- log(mammal.data$homeRange)
names(homerange) <- row.names(mammal.data)

str(mammal.tree)
# - attr(*, "class")= chr "phylo"

# calculate PIC continuous data
pic.bodymass <- pic(bodymass, mammal.tree, scaled = TRUE, var.contrasts = TRUE)
pic.homerange <- pic(homerange, mammal.tree, scaled = TRUE, var.contrasts = TRUE)

# plot PICed value and call the regression
fit.pic <- lm(pic.homerange~pic.bodymass + 0)
plot(pic.bodymass, pic.homerange)
abline(fit.pic, lwd = 2, lty = "dashed", col = "red")

# find the best-fitting evolutionary model
fitContinuous(mammal.tree, bodymass, model = "BM")
# AICc = 154.417886, smallest value ~ best-fitting model
fitContinuous(mammal.tree, bodymass, model = "OU")
# AICc = 155.815161
fitContinuous(mammal.tree, bodymass, model = "EB")
# AICc = 156.690459
fitContinuous(mammal.tree, bodymass, model = "rate_trend")
# AICc = 155.920991

# reconstruct the ancestral states of bodymass, specify the model
ace.bodymass <- ace(bodymass, mammal.tree, model = "BM")

# plot the reconstructed tree
contMap(mammal.tree, bodymass, method="user", anc.states = ace.bodymass$ace)

plot(mammal.tree)
nodelabels()

#### Discret data ####

data("eel.data")
data("eel.tree")
plot(eel.tree)

View(eel.data)

# isolate the binary character
feed_mode <- setNames(eel.data[,1], rownames(eel.data))
names(feed_mode) <- row.names(eel.data)
maxTL <- eel.data$Max_TL_cm
names(maxTL) <- row.names(eel.data)

#find the best fitting model
fitDiscrete(eel.tree, feed_mode, model = "ER")
# AICc = 76.115158
fitDiscrete(eel.tree, feed_mode, model = "SYM")
# AICc = 76.115158
fitDiscrete(eel.tree, feed_mode, model = "ARD")
# AICc = 78.212670
fitDiscrete(eel.tree, feed_mode, model = "meristic")
# AICc = 76.115158
# If tie, choose the model with least parameters
# ER only has 1 free parameter

# reconstruct the ancestral states of feed_mode, specify the data type as discret
fit.ace <- ace(feed_mode, eel.tree, type = "discret", model = "ER")

# Plot the reconstructed tree
cols <- c("blue", "red")
plotTree(eel.tree, fsize = 0.7, ftype = "i", lwd = 1, offset = 0.5)

nodelabels(node=1:eel.tree$Nnode+Ntip(eel.tree),pie=fit.ace$lik.anc,piecol=cols,cex=0.4)

tiplabels(pie=to.matrix(feed.mode[eel.tree$tip.label],levels(feed.mode)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)

#### Simulate BM model ####
# how does sig2 effects the trait value
simBMphylo(61, 500, 1)
simBMphylo(61, 500, 10)
simBMphylo(61, 500, 100)



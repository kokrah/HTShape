## ----setup, include = FALSE, echo=FALSE------------------------
# R options
options(width = 65)
options(continue=" ")
options(warn=-1)
set.seed(42)

# R chunk options.
require(knitr)
opts_chunk$set(concordance=TRUE, prompt=TRUE, comment=NA, size="small")

# smoothScatter(rnorm(100), rnorm(100))

## ----datasets--------------------------------------------------
library(shape)
data(examplesData) # Load datasets.
names(examplesData)

## ----filterPickrell1-------------------------------------------
filterCounts <- function (pcounts, thresh, minSamples) {
  cpm <- t(t(pcounts) /  colSums(pcounts)) * 1e+06
  keep <- rowSums(cpm > thresh) >= minSamples
  filteredPcounts <- pcounts[keep, ]
  filteredPcounts
}

## ----filterPickrell2-------------------------------------------
counts <- examplesData$pickrell$exprs
gender <- examplesData$pickrell$cond 
(tab <- table(gender))
pcounts <- counts + 1 # pseudo-counts
minSamples <- min(tab)
dim(pcounts) # Before filteration.
pcounts <- filterCounts(pcounts, 1, minSamples)
dim(pcounts) # After filteration. 

## ----normalizePickrell-----------------------------------------
ref <- exp(rowMeans(log(pcounts)))
deseqScal <- apply(pcounts / ref, 2, median)
pcounts <- t(t(pcounts) / deseqScal)
y <- log2(pcounts)

## ----shapePickrell---------------------------------------------
# Compute the L-skew (t3) and L-kurt (t4) of each gene.
res <- fitShape(y) 
class(res)
lapply(res, head, n=3)

## ----dvalsPickrell---------------------------------------------
# Compute d-values
t3 <- res$lrats[, "t3"] # Grab L-skew estimates. 
t4 <- res$lrats[, "t4"] # Grab L-kurt estimates. 
dvals <- computeDvals(t3, t4)

## ----soplotPickrell, fig.width=5, fig.height=3.5, fig.align='center'----
# Symmetry-Outlier plot.
plotSO(t3, t4, dataName="Pickrell (No Gender Adjustment)", verbose = TRUE)

# Pick volatile / outlier genes.
sel <- which(dvals < 0.0001) # select 0.01% cutoff

# Seperate outlier genes into 2 groups for illustration purposes
blueGroup <- sel[abs(t3[sel]) > 0.3]
redGroup <- sel[abs(t3[sel]) <= 0.3]
points(t3[blueGroup], t4[blueGroup], cex=0.5, col="blue")
points(t3[redGroup], t4[redGroup], cex=0.5, col="red")

## ----redGroupPickrell, fig.width=5, fig.height=2.3, fig.align='center', echo=FALSE----
par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5))

exprsRed <- y[redGroup, ]

# Shift and scale epression levels.
exprsRed <- t(scale(t(exprsRed), TRUE, TRUE))

# Add canvas.
plot(0, 0, pch="",
     xlim=c(0.75, nrow(exprsRed)), ylim=range(exprsRed), 
     ylab="Standard Expr. Level", xlab="Genes", 
     cex.axis=0.7, cex.lab=0.7, xaxt="n")
axis(1, 1:nrow(exprsRed), 1:nrow(exprsRed), cex.axis=0.7)
abline(h=0, col="gray", lty=3)

# Add color.
gcol <- ifelse(gender == "female", "black", "red")

for(k in 1:nrow(exprsRed)){ # Add data points.
  points(jitter(rep(k, ncol(exprsRed)), amount=0.2), 
         exprsRed[k,], cex=0.35, col=gcol)
} 
text(1:nrow(exprsRed)-0.45, rep(0, nrow(exprsRed)), 
     rownames(exprsRed), srt=90, cex=0.6)

## ----blueGroupPickrell, fig.width=1.8, fig.height=2.3, fig.align='center', echo=FALSE----
par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5))

exprsBlue <- y[blueGroup,]

# Shift and scale epression levels.
exprsBlue <- t(scale(t(exprsBlue), TRUE, TRUE))

# Add canvas.
plot(0, 0, pch="", 
     xlim=c(2, 4.4), ylim=range(exprsBlue),
     ylab="Standard Expr. Level", xlab="Genes",
     cex.axis=0.7, cex.lab=0.7, xaxt="n")
axis(1, c(2.75, 4), 1:2, cex.axis=0.7) 
abline(h=0, col="gray", lty=3)

for(k in 1:nrow(exprsBlue)){
  if (k==1) j <- 2.5 else j <- 3.75
  
  points(jitter(rep(j, ncol(exprsBlue)), amount=0.2), 
         exprsBlue[k,], cex=0.35, col=gcol)
} 

# Add boxplots
boxplot(t(exprsBlue), at=c(3, 4.25), add=TRUE, 
        boxwex=0.4, cex=0.5, xlab="", ylab="", 
        xaxt="n", yaxt="n")
text(c(2.5, 3.75)-0.4, rep(0, nrow(exprsBlue)), 
     rownames(exprsBlue), srt=90, cex=0.6)

## ----volatilePickrell, fig.width=5, fig.height=4, fig.align='center'----
head(computeDvals(t3, t4, plot=TRUE))

## ----hammoudData-----------------------------------------------
hammoud <- examplesData$hammoud
rpkm <- hammoud$exprs
cond <- hammoud$cond
(tab <- table(cond))

## ----filterGenesHammoud, echo=FALSE----------------------------
minSamples <- min(tab)

# dim(rpkm) # Before filteration.
keep <- rowSums(rpkm > 1) >= minSamples
rpkm <- rpkm[keep, ]  
# dim(rpkm) # After filteration.

log2RPKM <- log2(rpkm + 1)

compute.residuals <- function (exprs, cond) {
    design <- model.matrix(~0 + as.factor(cond))
    beta.hat <- solve(t(design) %*% design) %*% t(design) %*% 
        t(as.matrix(exprs))
    res <- exprs - t(design %*% beta.hat)
    as.data.frame(res)
}

resids <- compute.residuals(log2RPKM, cond)

## ----Hammoud1, fig.width=7, fig.height=3, fig.align='center'----
par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5), mfrow=c(1, 2))

res1 <- fitShape(log2RPKM) # Without adjustment for cell type.
res2 <- fitShape(resids) # With adjustment for cell type.

plotSO(res1$lrats[, "t3"], res1$lrats[, "t4"], 
       dataName = "Hammoud (No adjustment)") 
plotSO(res2$lrats[, "t3"], res2$lrats[, "t4"], 
       dataName="Hammoud (Cell type adjusted)")       

## ----Hammoud2, fig.width=7, fig.height=3, fig.align='center', echo=FALSE----
par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5), mfrow=c(1, 2))
# Compute fold change
aveSpermatids <- rowMeans(log2RPKM[, cond=="Spermatids"])
aveSpermatocytes <- rowMeans(log2RPKM[, cond=="Spermatocytes"]) 
lfc <-  aveSpermatids - aveSpermatocytes

# Fold change vs. L-kurt (no adjustment for cell type)
rang <- range(c(res1$lrats[, "t4"], res2$lrats[, "t4"]))

plot(lfc, res1$lrats[, "t4"], pch=".", col="gray",
     ylab="L-kurt", xlab="log2 fold change",
     main="No adjustment", ylim=rang)
abline(h=0, lty=2, col="red")

# Fold change vs. L-kurt (adjusted for cell type)
plot(lfc, res2$lrats[, "t4"], pch=".", col="gray",
     ylab="L-kurt", xlab="log2 fold change",
     main="Resids", ylim=rang)
abline(h=0, lty=2, col="red")

## ----Hammoud3, fig.width=7, fig.height=2.5, fig.align='center', echo=FALSE----
set.seed(999)
par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5))

sel <- res1$lrats[, "t4"] < -0.2 # How many ?
subLog2RPKM <- log2RPKM[sample(which(sel), 100), ]
subLog2RPKM <- t(scale(t(subLog2RPKM), TRUE, TRUE))

matplot(subLog2RPKM[1:50,], pch=16, cex=0.5,
        col=ifelse(cond=="Spermatids", "blue", "red"),
        ylab="Standard Exprs. Level", 
        xlab="50 random genes (L-kurt < -0.2)",
        main="(1) Random 50")
abline(v=1:50, col="gray", lty=3)
abline(h=0, col="gray", lty=2)

## ----Hammoud4, fig.width=7, fig.height=2.5, fig.align='center', echo=FALSE----
par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5))

matplot(subLog2RPKM[51:100,], pch=16, cex=0.5,
        col=ifelse(cond=="Spermatids", "blue", "red"),
        ylab="Standard Exprs. Level", 
        xlab="50 random genes (L-kurt < -0.2)",
        main="(2) Random 50", xaxt="none")
axis(1, 1:50, labels = 51:100)
abline(v=1:50, col="gray", lty=3)
abline(h=0, col="gray", lty=2)

## ----bottomlyData----------------------------------------------
bottomly <- examplesData$bottomly
counts <- bottomly$exprs
cond <- bottomly$cond
(tab <- table(cond))

## ----filterGenesBottomly, echo=FALSE---------------------------
minSamples <- min(tab)

# dim(counts) # Before filteration.
counts <- filterCounts(counts, thresh=1, minSamples=minSamples)
# dim(counts)# After filteration.

pcounts <- counts + 1 # Pseudo-counts.

# Normalize (deseq method).
ref <- exp(rowMeans(log(pcounts)))
deseqScal <- apply(pcounts / ref, 2, median)
pcounts <- t(t(pcounts) / deseqScal)

log2pcounts <- log2(pcounts)
resids <- compute.residuals(log2pcounts, cond)

## ----Bottomly1, fig.width=7, fig.height=3, fig.align='center'----
par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5), mfrow=c(1, 2))

res1 <- fitShape(log2pcounts) # Without adjustment for strain.
res2 <- fitShape(resids) # With adjustment for strain.

plotSO(res1$lrats[, "t3"], res1$lrats[, "t4"], 
       dataName = "Bottomly (No adjustment)") 
plotSO(res2$lrats[, "t3"], res2$lrats[, "t4"], 
       dataName="Bottomly (Strain type adjusted)")

## ----bottomly2, fig.width=7, fig.height=3, fig.align='center', echo=FALSE----
par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0.5), mfrow=c(1, 2))
# Compute fold change
aveC57BL <- rowMeans(log2pcounts[, cond=="C57BL/6J"])
aveDBA <- rowMeans(log2pcounts[, cond=="DBA/2J"]) 
lfc <-  aveC57BL - aveDBA

# Fold change vs. L-kurt (no adjustment for strain type)
rang <- range(c(res1$lrats[, "t4"], res2$lrats[, "t4"]))

plot(lfc, res1$lrats[, "t4"], pch=".", col="gray",
     ylab="L-kurt", xlab="log2 fold change",
     main="No adjustment", ylim=rang)
abline(h=0, lty=2, col="red")

# Fold change vs. L-kurt (adjusted for strain type)
plot(lfc, res2$lrats[, "t4"], pch=".", col="gray",
     ylab="L-kurt", xlab="log2 fold change",
     main="Resids", ylim=c(rang))
abline(h=0, lty=2, col="red")

## ----sessionInfo-----------------------------------------------
sessionInfo()


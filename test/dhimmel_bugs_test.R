library(pheatmap)
pheatmap.blank <- function(...) pheatmap(cluster_rows = F, cluster_cols = F, show_rownames = F, ...)

source("GSC.R")
source("test/GSC2.R")
load("test/se.mat.RData")

UnderrepresentationWeight <- function(mat, functionName = "GSC") {
  # Returns underrepresentation weight for each column
  cat ("# Computing distance matrix...\n")
  col.dist <- stats::dist(t(mat), method = 'binary')
  cat ("# Hierarchical clustering...\n")
  col.clust <- hclust(col.dist, method = 'ward.D2') # change ward.D2 to ward if you;re running an earlier version than 3.0 of R
  cat ("# Computing the weights...\n")
  col.dendro <- as.dendrogram(col.clust)
  get(functionName)(col.dendro)
}



# Nan Problem -------------------------------------------------------------

UnderrepresentationWeight(se.mat[, 49:84]) # Minimum reproducible example
mat <- se.mat[,c(49:60,84)] # Let's stripdown the number of column
col.dist <- stats::dist(t(mat), method = 'binary')
col.clust <- hclust(col.dist, method = 'ward.D2')
col.dendro <- as.dendrogram(col.clust)
plot(col.dendro)
pheatmap(1-col.dist, clustering_distance_rows = col.dist, clustering_distance_cols = col.dist)
# Here is the sungularity that leads to NaN:
any(apply(se.mat[,c(49,84)], 1, diff) != 0) # Same vector... -> distance of 0 -> clustering with h = 0 for the two leaves -> division by 0 -> bug. (solved)

# Second "problem": 
# distance == 1 when the distance is computed with method = 'binary' with everybody else for these 3 guys (in a reduced matrix at least).
# This is gracefully handled by the algorithm, the weights are just equals. (because the nodes are at exactly the same heigths in the dendrogram)
mat3 <- se.mat[,colnames(se.mat) %in% paste0("C000", c(2880, 2624, 2631))]
pheatmap.blank(mat3) #It occurs in this case because each of these guys have too few, rare positives, so the distance == 1 with everybody else when using the binary distance. (no positive in common)

# Check the NaN fix:
mat1 <- mat
mat2 <- mat
# modify slightly the second matrix to make the first vector a bit different than the last:
mat2[2,1] <- 1
ttt1 <- UnderrepresentationWeight(mat1)
ttt2 <- UnderrepresentationWeight(mat2)
cbind(ttt1,ttt2) # The weigths at the singularity (ttt1) are close to the weights at the near-singularity (ttt2). The comportment at the limit is verified, all good.



# profiling -----------------------------------------------

library(microbenchmark)
library(ggplot2)

createRandomDendrogram <- function(n = 300) {
  rmat <- matrix(nc = n, nr = n, runif(n^2, 0, 1))
  col.dist <- as.dist(rmat + t(rmat))
  col.hc <- hclust(col.dist, method = "ward.D2")
  as.dendrogram(col.hc)
}

mbtimes <- sapply(1:15 * 200, function(x) {
  dd <- createRandomDendrogram(x)
  microbenchmark({ cat("# n = ", x, "\n"); GSC(dd)}, times = 3)
  })
mbtimes2 <- sapply(1:15 * 200, function(x) {
  dd <- createRandomDendrogram(x)
  microbenchmark({ cat("# n = ", x, "\n"); GSC(dd)}, times = 3)
})
dfbench <- rbind(
  data.frame(type = "normal" , 
             do.call(rbind, lapply(1:dim(mbtimes)[2], function(x) data.frame(n = x * 200, time = mbtimes[[2,x]]/10e6)))),
  data.frame(type = "verbose",
             do.call(rbind, lapply(1:dim(mbtimes)[2], function(x) data.frame(n = x * 200, time = mbtimes2[[2,x]]/10e6)))))
pdf("tests/profiling.pdf")
qplot(data = dfbench, x = n, y = time, color = type, shape = I(19), alpha = I(0.5), geom = c("point", "smooth"), 
      title = "profiling of the GSC algo, two versions", ylab = "time (ms)", xlab = "elements", method = "lm") +
  theme_bw()
dev.off()



# recursion check ----------------------------------------------------------

# Get the linear relationship between the size of the dendrogram and n, its number of leaves.
n <- 1:30 * 50
size <- sapply(n, function(x) object.size(createRandomDendrogram(x)))
qplot(n, size, geom = c("point", "smooth"))
lm(size ~ n)

# Test the new dude.
# Create a bid dendrogram
bigDendo <- do.call(merge, lapply(rep(1000, 5), createRandomDendrogram))
print(bigDendo)
# artificially put the number of allowed exressions low, retaining the initial pre-testing parameter.
print(init.expr <- options("expressions" = 50)[[1]])
# execute the GSC algo
ttt <- GSC(bigDendo)
str(ttt)
# reset the option while checking it was left untouched by the algo.
print(options("expressions" = init.expr))


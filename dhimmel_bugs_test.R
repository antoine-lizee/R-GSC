library(pheatmap)
pheatmap.blank <- function(...) pheatmap(cluster_rows = F, cluster_cols = F, show_rownames = F, ...)

source("GSC.R")
source("GSC2.R")
load("se.mat.RData")

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

mbtimes <- sapply(1:15 * 200, function(x) microbenchmark({ cat("# n = ", x*200, "\n"); UnderrepresentationWeight(se.mat[,1:x])}, times = 3))
mbtimes2 <- sapply(1:15 * 200, function(x) microbenchmark({ cat("# n = ", x*200, "\n"); UnderrepresentationWeight(se.mat[,1:x], "GSC2")}, times = 3))
dfbench <- rbind(
  data.frame(type = "normal" , 
             do.call(rbind, lapply(1:dim(mbtimes)[2], function(x) data.frame(n = x * 200, time = mbtimes[[2,x]]/10e9)))),
  data.frame(type = "verbose",
             do.call(rbind, lapply(1:dim(mbtimes)[2], function(x) data.frame(n = x * 200, time = mbtimes2[[2,x]]/10e9)))))
pdf("profiling.pdf")
qplot(data = dfbench, x = n, y = time / 10e9, color = type, shape = I(19), alpha = I(0.5), geom = c("point", "smooth"), 
      title = "profiling of the GSC algo, two versions", ylab = "time (s)", xlab = "elements") +
  theme_bw()
dev.off()

(tt.exp <- nls(time/10e9 ~ ca*exp(cb*n) + cc*n + cd, data = dfbench,
          start = list(ca = 0.1, cb = 0.001, cc = 0.1, cd = 0.1),
          trace = T))
(tt.square <- nls(time/10e9 ~ ca*n^2 + cc*n + cd, data = dfbench,
              start = list(ca = 0.1, cc = 0.1, cd = 0.1),
              trace = T))
(ttlm.square <- lm(time/10e9 ~ I(n^2) + n , data = dfbench))
dfbench$time.exp <- predict(tt.exp)
dfbench$time.square <- predict(ttlm.square)

pdf("profiling2.pdf")
qplot(data = dfbench, x = n, y = time/10e9, shape = I(19), alpha = I(0.5), geom = c("point"), 
      title = "profiling of the GSC algo with fit on top", ylab = "time (s)", xlab = "elements") +
  geom_line(aes(x = n, y = time.square), color = "red")
dev.off()



# recursion test ----------------------------------------------------------



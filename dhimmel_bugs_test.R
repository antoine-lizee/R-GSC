source ("GSC.R")
load("se.mat.RData")

UnderrepresentationWeight <- function(mat) {
  # Returns underrepresentation weight for each column
  col.dist <- stats::dist(t(mat), method = 'binary')
  col.clust <- hclust(col.dist, method = 'ward')
  col.dendro <- as.dendrogram(col.clust)
  plot(col.dendro)
  GSC(col.dendro)
}


# Nan Problem -------------------------------------------------------------

UnderrepresentationWeight(se.mat[, 49:84]) # Minimum reproducible example
mat <- se.mat[,c(49:60,84)] # Let's stripdown the number of column
col.dist <- stats::dist(t(mat), method = 'binary')
col.clust <- hclust(col.dist, method = 'ward')
col.dendro <- as.dendrogram(col.clust)
plot(col.dendro)
pheatmap(1-col.dist, clustering_distance_rows = col.dist, clustering_distance_cols = col.dist)
# Here is the sungularity that leads to NaN:
any(apply(se.mat[,c(49,84)], 1, diff) != 0) # Same vector... -> distance of 0 -> clustering with h = 0 for the two leaves -> division by 0 -> bug. (solved)

# Second "problem": 
# distance == 1 when the distance is computed with method = 'binary' with everybody else for these 3 guys (in a reduced matrix at least).
# This is gracefully handled by the algorithm, the weights are just equals. (because the nodes are at exactly the same heigths in the dendrogram)
mat3 <- se.mat[,colnames(se.mat) %in% paste0("C000", c(2880, 2624, 2631))]
pheatmap(mat3) #It occurs in this case because each of these guys have too few, rare positives, so the distance == 1 with everybody else when using the binary distance. (no positive in common)

# Check the NaN fix:
mat1 <- mat
mat2 <- mat
# modify slightly the second matrix to make the first vector a bit different than the last:
mat2[2,1] <- 1
ttt1 <- UnderrepresentationWeight(mat1)
ttt2 <- UnderrepresentationWeight(mat2)
cbind(ttt1,ttt2) # The weigths at the singularity (ttt1) are close to the weights at the near-singularity (ttt2). The comportment at the limit is verified, all good.



# recursivity limit problem -----------------------------------------------



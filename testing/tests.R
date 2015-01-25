# This is the actual example in the paper appendix [doi:10.1016/0022-2836(94)90012-4].
# The reader should be warned that the results in the paper are
# wrong, due to some imprecision.

testhc <- list("merge" = matrix(nrow = 3, byrow = T, c(-1,-2,-3,1, -4,2)),
               "height" = c(20, 50, 80),
               "order" = c(1,2,3,4),
               "labels" = LETTERS[1:4])
attr(testhc, "class") <- "hclust"
testdend <- as.dendrogram(testhc)
plot(testdend)
print(GSC(testdend))


# A more real-looking use-case:

hc <- hclust(dist(mtcars))
dd <- as.dendrogram(hc)
plot(as.dendrogram(hc))
ddGSCs <- GSC(dd) #compute the weights
hc$labels[hc$order] <- paste(names(ddGSCs), sprintf("%.1f", ddGSCs), sep = " - ") #add them to the label names
plot(as.dendrogram(hc)) #plot them to have a look.

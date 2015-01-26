# Less elegant version of GSC that divide by two the recursion depth, by avoiding calling the nicer helper function.
# Profiling show no improvement whatsoever, so the code is kept here for reference only.

GSCrec2 <- function(dd) {
  h0 <- attr(dd, "height")
  
  ddcAttrs <- attributes(dd[[1]])
  if (ddcAttrs$members == 1) { # we have a leaf.
    stopifnot(ddcAttrs$height == 0, ddcAttrs$leaf) # leaf should be at height 0, and have the leaf signature
    # Create the coeffs
    GSCs1 <- h0
    names(GSCs1) <- ddcAttrs$label
  } else { # we have a tree
    gscs <- GSCrec(dd[[1]]) + .Machine$double.eps
    GSCs1 <- gscs + gscs * (h0-ddcAttrs$height) / sum(gscs) - .Machine$double.eps
  }
  
  ddcAttrs <- attributes(dd[[2]])
  if (ddcAttrs$members == 1) { # we have a leaf.
    stopifnot(ddcAttrs$height == 0, ddcAttrs$leaf) # leaf should be at height 0, and have the leaf signature
    # Create the coeffs
    GSCs2 <- h0
    names(GSCs2) <- ddcAttrs$label
  } else { # we have a tree
    gscs <- GSCrec(dd[[2]]) + .Machine$double.eps
    GSCs2 <- gscs + gscs * (h0-ddcAttrs$height) / sum(gscs) - .Machine$double.eps
  }
  
  return(c(GSCs1, GSCs2))
}

# Main function to do some checks + normalize the weights at the end.
GSC2 <- function(dd) {
  if(class(dd) != "dendrogram") {
    stop("Argument should be a dendrogram object. Please use as.dendrogram()")
  }
  GSCs <- GSCrec2(dd)
  return(GSCs / sum(GSCs) * length(GSCs))
}


# test --------------------------------------------------------------------

if (test <- FALSE) { # Change to TRUE to launch testing code
  
  # This is the actual example of the paper. The reader should be warned that the results in the paper are
  # wrong, due to some imprecision.
  testhc <- list("merge" = matrix(nrow = 3, byrow = T, c(-1,-2,-3,1, -4,2)),
                 "height" = c(20, 50, 80),
                 "order" = c(1,2,3,4),
                 "labels" = LETTERS[1:4])
  attr(testhc, "class") <- "hclust"
  testdend <- as.dendrogram(testhc)
  plot(testdend)
  print(GSC2(testdend))
  
  # A more real-looking use-case:
  hc <- hclust(dist(mtcars))
  dd <- as.dendrogram(hc)
  plot(as.dendrogram(hc))
  ddGSCs <- GSC2(dd) #compute the weights
  hc$labels[hc$order] <- paste(names(ddGSCs), sprintf("%.1f", ddGSCs), sep = " - ") #add them to the label names
  plot(as.dendrogram(hc)) #plot them to have a look.
}

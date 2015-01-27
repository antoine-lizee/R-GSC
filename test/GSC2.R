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


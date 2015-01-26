# Implementation of the Gerstein-Sonnhammer-Chothia algorithm
#
# This script makes discoverable the funtion "GSC" that takes as argument a dendrogram and returns
# a named vector of weight for each of the leaves of the input dendrogram. These weigths are supposed
# to account for similarity between objects when considering them in a model. The more similar, the 
# lower the weights.
# The initial description of the implementation can be found here: https://pdf.yt/d/Sx3jMbr8vANgxAej/download
# As per request for: http://thinklab.org/p/rephetio/seeking-an-open-source-implementation-of-the-gerstein-sonnhammer-chothia-algorithm/25
#
# Copyright Antoine Lizee 2015/01 antoine.lizee@gmail.com


# Algorithm ---------------------------------------------------------------------

# getGSCs is a handler for the main recursve funtion below. It actually does all the work:
# It applies the main function to the child tree and does the main weighing computation,
# or create the right data structure (and check consistency) for a leaf.
getGSCs <- function(ddchild, parentHeight) {
  ddcAttrs <- attributes(ddchild)
  if (ddcAttrs$members == 1) { # we have a leaf.
    stopifnot(ddcAttrs$height == 0, ddcAttrs$leaf) # leaf should be at height 0, and have the leaf signature
    # Create the coeffs
    GSCs <- parentHeight
    names(GSCs) <- ddcAttrs$label
  } else { # we have a tree
    gscs <- GSCrec(ddchild) + .Machine$double.eps
    GSCs <- gscs + gscs * (parentHeight-ddcAttrs$height) / sum(gscs) - .Machine$double.eps
  }
  GSCs
}

# Main recursive function
GSCrec <- function(dd) {
  h0 <- attr(dd, "height")
  return(c(getGSCs(dd[[1]], h0), getGSCs(dd[[2]], h0)))
}

# Main function to do some checks + normalize the weights at the end.
GSC <- function(dd) {
  if(class(dd) != "dendrogram") {
    stop("Argument should be a dendrogram object. Please use as.dendrogram()")
  }
  GSCs <- GSCrec(dd)
  return(GSCs / sum(GSCs) * length(GSCs))
}


#####

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



# Tests -------------------------------------------------------------------

if (test <- TRUE) { # Change to TRUE to launch testing code
  
  # This is the actual example of the paper. The reader should be warned that the results in the paper are
  # wrong, due to some imprecision.
  testhc <- list("merge" = matrix(nrow = 3, byrow = T, c(-1,-2,-3,1, -4,2)),
                 "height" = c(20, 50, 80),
                 "order" = c(1,2,3,4),
                 "labels" = LETTERS[1:4])
  attr(testhc, "class") <- "hclust"
  testdend <- as.dendrogram(testhc)
  plot(testdend)
  print(GSC(testdend))
  print(GSC2(testdend))
  
  # A more real-looking use-case:
  hc <- hclust(dist(mtcars))
  dd <- as.dendrogram(hc)
  plot(as.dendrogram(hc))
  ddGSCs <- GSC(dd) #compute the weights
  hc$labels[hc$order] <- paste(names(ddGSCs), sprintf("%.1f", ddGSCs), sep = " - ") #add them to the label names
  plot(as.dendrogram(hc)) #plot them to have a look.
    
}

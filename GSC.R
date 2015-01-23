# Implementation of the Gerstein-Sonnhammer-Chothia algorithm
#
# This script makes discoverable the funtion "GSC" that takes as argument a dendrogram and returns
# a named vector of weight for each of the leaves of the input dendrogram. These weigths are supposed
# to account for similarity between objects when considering them in a model. The more similar, the 
# lower the weights.
# The initial description of the implementation can be found here: https://pdf.yt/d/Sx3jMbr8vANgxAej/download
# As per request for: http://thinklab.org/p/rephetio/seeking-an-open-source-implementation-of-the-gerstein-sonnhammer-chothia-algorithm/25


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
    gscs <- GSCrec(ddchild)
    GSCs <- gscs + gscs * (parentHeight-ddcAttrs$height) / sum(gscs)
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



# Tests -------------------------------------------------------------------

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
  print(GSC(testdend))
  
  # A more real-looking use-case:
  hc <- hclust(dist(mtcars))
  dd <- as.dendrogram(hc)
  plot(as.dendrogram(hc))
  ddGSCs <- GSC(dd) #compute the weights
  hc$labels[hc$order] <- paste(names(ddGSCs), sprintf("%.1f", ddGSCs), sep = " - ") #add them to the label names
  plot(as.dendrogram(hc)) #plot them to have a look.
   
}

# Implementation of the Gerstein-Sonnhammer-Chothia algorithm.
#
# This script makes discoverable the funtion "GSC" that takes as argument a dendrogram and returns
# a named vector of weight for each of the leaves of the input dendrogram. These weigths are supposed
# to account for similarity between objects when considering them in a model. The more similar, the 
# lower the weights.
# The algorithm was initially descriped in:
# > Gerstein M, Sonnhammer ELL, Chothian C (1994) Volume Changes in Protein Evolution.
# > Journal of Molecular Biology. doi:10.1016/0022-2836(94)90012-4
#
# Copyright Antoine Lizee 2015/01 antoine.lizee@gmail.com

getGSCs <- function(ddchild, parentHeight) {
  # getGSCs is a handler for the main recursve funtion below. It actually does all the work:
  # It applies the main function to the child tree and does the main weighing computation,
  # or create the right data structure (and check consistency) for a leaf.
  
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

GSCrec <- function(dd) {
  # Main recursive function
  
  h0 <- attr(dd, "height")
  return(c(getGSCs(dd[[1]], h0), getGSCs(dd[[2]], h0)))
}


GSC <- function(dd) {
  # Main function to do some checks + normalize the weights at the end.
  
  # Check the class
  if(class(dd) != "dendrogram") {
    stop("Argument should be a dendrogram object. Please use as.dendrogram()")
  }
  
  # Check the number of elements and adjust the max recursive depth accordingly.
  # Make sure that the dendrogram is not bigger than it says:
  nElem <- as.integer(attr(dd, "members"))
  nElem.size <- (object.size(dd) + 320) / 1224 # cf dhimmel_bugs_test for determining this.
  if ( length(nElem) == 0 || is.na(nElem) || nElem < nElem.size /2) { 
    warning("dendrogram seems bigger than the static 'members' attribute of the root node says. Guessing...")
    nElem <- c(round(nElem.size * 1.2))
  }
  
  # Set the options to increase the max depth
  opt.old <- options("expressions")[[1]]
  opt.new <- nElem * 2 + 500 # O(log(nElem)), worst case nElem -> nElem * 2 we're large.
  if (opt.new > opt.old)
    options("expressions" = opt.new)
  
  # Launch the rest inside a tryCatch to make sure we reset the options.
  tryCatch( {
    GSCs <- GSCrec(dd)
    return(GSCs / sum(GSCs) * length(GSCs))
  }, 
  finally = {options("expressions" = opt.old)}
  )
}

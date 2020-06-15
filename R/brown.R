

#' @title Brown Clustering of Words
#' @description Cluster words in groups by assuming that the usage and probabilities of occurrence of words 
#' are based on the cluster that the previous word belongs to: \eqn{P(w_i|w_{i-1}) = P(w_i|c_{i}) P(c_i|c_{i-1})}.\cr
#' Wrapper for the \url{https://github.com/percyliang/brown-cluster} C++ implementation.
#' @param x either the path to a file or a character vector of words where words are separated by spaces
#' @param initC the number of clusters
#' @param min_occur how many times a word should appear to be considered. Defaults to 5.
#' @param num_threads number of threads to use for parallel computation. Defaults to 1.
#' @return an object of class 'brown' which is a list with elements
#' \itemize{
#' \item{initC: the \code{initC} argument provided or if the number of phrases found in \code{x} is smaller than the amount provided, the lowest of these 2 values}
#' \item{min_occur: the \code{min_occur} argument provided}
#' \item{clusters: a data.frame with 1 row per word providing the path to the cluster, the frequency of occurrence of the word and the kullback-leiber divergence to the left and right cluster path}
#' \item{collocations: a data.frame with collocations found in \code{x}}
#' }
#' @references \url{https://github.com/percyliang/brown-cluster}, Brown, et al.: Class-Based n-gram Models of Natural Language
#' @export
#' @examples 
#' ## Take some data and standardise it a bit
#' library(udpipe)
#' data(brussels_reviews, package = "udpipe")
#' x <- subset(brussels_reviews, language == "nl")
#' x <- tolower(x$feedback)
#' 
#' ## Perform Brown clustering
#' wordgroups <- brown(x, initC = 50, min_occur = 3)
#' 
#' ## Inspect the output
#' head(wordgroups$clusters)
#' fat <- split(wordgroups$clusters$word, wordgroups$clusters$path)
#' fat <- lapply(fat, txt_collapse)
#' fat
brown <- function(x, initC = 100, min_occur = 5, num_threads = 1L){
  initC <- as.integer(initC)
  min_occur <- as.integer(min_occur)
  num_threads <- as.integer(num_threads)
  plen <- 1L
  x <- as.character(x)
  newfile <- tempfile(pattern = "brown_", fileext = ".txt")
  on.exit({
    if(file.exists(newfile)) file.remove(newfile)
  })
  if(length(x) == 1 && file.exists(x)){
    file.copy(from = x, to = newfile)
  }else{
    writeLines(x, con = newfile)
  }
  out <- cluster_brown(text_file = newfile, 
                       output_dir = tempdir(), 
                       min_occur = min_occur, initC = initC, plen = plen,
                       num_threads = num_threads)
  out$clusters <- out$clusters[order(out$clusters$path, -out$clusters$freq, decreasing = FALSE), ]
  out
}
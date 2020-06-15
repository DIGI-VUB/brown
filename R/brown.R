

#' @title Brown Clustering of Words
#' @description Cluster words in groups by assuming that the usage and probabilities of occurrence of words 
#' are based on the cluster that the previous word belongs to: \eqn{P(w_i|w_{i-1}) = P(w_i|c_{i}) P(c_i|c_{i-1})}.\cr
#' Wrapper for the \url{https://github.com/percyliang/brown-cluster} C++ implementation.
#' @param x either the path to a file or a character vector of words where words are separated by spaces
#' @param initC the number of clusters
#' @param min_occur how many times a word should appear to be considered. Defaults to 1.
#' @param vocabulary optionally, a character vector of words to use as the vocabulary
#' @param num_threads number of threads to use for parallel computation. Defaults to 1.
#' @param trace logical indicating to print the evolution of the algorithm. Defaults to FALSE.
#' @param ... passed on to \code{readLines} in case \code{x} is a file
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
#' ## Take some data and standardise it a bit, put it in a vector of words
#' library(udpipe)
#' data(brussels_reviews, package = "udpipe")
#' x <- subset(brussels_reviews, language == "nl")
#' x <- strsplit(x$feedback, "[[:space:][:punct:]]+")
#' x <- unlist(x)
#' x <- tolower(x)
#' 
#' ## Perform Brown clustering
#' wordgroups <- brown(x, initC = 50, min_occur = 3)
#' 
#' ## Inspect the output
#' head(wordgroups$clusters)
#' fat <- split(wordgroups$clusters$word, wordgroups$clusters$path)
#' fat <- lapply(fat, paste, collapse = " ")
#' fat
#' 
#' 
#' ## See the evolution of the algorithm by setting trace = TRUE
#' txt <- c("the cat chased the mouse", 
#'          "the dog chased the cat",
#'          "the mouse chased the dog")
#' txt <- strsplit(txt, " ")
#' txt <- unlist(txt)
#' wordgroups <- brown(txt, vocabulary = c("cat", "dog", "mouse", "chased"), 
#'                     initC = 4, min_occur = 1, trace = TRUE)
brown <- function(x, initC = 100, min_occur = 1, vocabulary = character(), num_threads = 1L, trace = FALSE, ...){
  initC <- as.integer(initC)
  min_occur <- as.integer(min_occur)
  num_threads <- as.integer(num_threads)
  plen <- 1L
  x <- as.character(x)
  if(length(x) == 1 && file.exists(x)){
    x <- readLines(x, ...)
    x <- strsplit(x, "[[:space:][:punct:]]+")
    x <- unlist(x)
  }
  if(trace){
    out <- cluster_brown(x = x, 
                         vocabulary = vocabulary, 
                         min_occur = min_occur, initC = initC, plen = plen,
                         num_threads = num_threads)
  }else{
    log <- utils::capture.output(out <- cluster_brown(x = x, 
                                                      vocabulary = vocabulary, 
                                                      min_occur = min_occur, initC = initC, plen = plen,
                                                      num_threads = num_threads))
  }
  out$clusters <- out$clusters[order(out$clusters$path, -out$clusters$freq, decreasing = FALSE), ]
  out
}
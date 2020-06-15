# brown

This repository contains an R package for performing clustering of words based on the [Brown algorithm](https://en.wikipedia.org/wiki/Brown_algorithm). 

- References:
    - Brown, et al.: Class-Based n-gram Models of Natural Language: http://acl.ldc.upenn.edu/J/J92/J92-4003.pdf
    - Liang: Semi-supervised learning for natural language processing: http://cs.stanford.edu/~pliang/papers/meng-thesis.pdf

### Installation

- For installing the development version of this package: `remotes::install_github("DIGIT-VUB/brown")`

Look to the documentation of the functions

```
help(package = "brown")
```

## Example

- Take some data and standardise it a bit

```{r}
library(udpipe)
data(brussels_reviews, package = "udpipe")
x <- subset(brussels_reviews, language == "nl")
x <- strsplit(x$feedback, "[[:space:][:punct:]]+")
x <- unlist(x)
x <- tolower(x)
```

- Perform Brown clustering and look at the clusters

```{r}
library(brown)
wordgroups <- brown(x, initC = 50, min_occur = 3)
str(wordgroups)
List of 4
 $ initC       : int 50
 $ min_occur   : int 3
 $ clusters    :'data.frame':	1002 obs. of  5 variables:
  ..$ word    : chr [1:1002] "is" "hadden" "beter" "deed" ...
  ..$ path    : chr [1:1002] "0000" "0000" "0000" "0000" ...
  ..$ freq    : int [1:1002] 509 54 6 6 6 6 5 5 5 3 ...
  ..$ kl_left : num [1:1002] 0.808 2.558 3.651 3.086 3.861 ...
  ..$ kl_right: num [1:1002] 0.685 1.494 3.778 2.893 2.7 ...
 $ collocations:'data.frame':	2500 obs. of  3 variables:
  ..$ term1      : chr [1:2500] "het" "het" "in" "in" ...
  ..$ term2      : chr [1:2500] "appartement" "centrum" "brussel" "de" ...
  ..$ collocation: num [1:2500] 0.0325 0.0183 0.015 0.0129 0.0107 ...
length(table(wordgroups$clusters$path))
[1] 50
```

```{r}
fat <- split(wordgroups$clusters$word, wordgroups$clusters$path)
fat <- lapply(fat, paste, collapse = " ")
fat
```


### DIGI

By DIGI: Brussels Platform for Digital Humanities: https://digi.research.vub.be

![](tools/logo.png)
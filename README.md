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
x <- tolower(x$feedback)
x <- gsub(pattern = "[[:space:][:punct:]]+", replacement = " ", x)
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
  ..$ word    : chr [1:1002] "doen" "gaan" "verblijven" "ondanks" ...
  ..$ path    : chr [1:1002] "000" "000" "000" "000" ...
  ..$ freq    : int [1:1002] 18 18 18 17 16 14 13 13 12 12 ...
  ..$ kl_left : num [1:1002] 1.71 1.84 1.81 2.01 2.81 ...
  ..$ kl_right: num [1:1002] 1.77 1.5 1.84 2.06 3.07 ...
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
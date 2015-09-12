HTShape: Shape analysis of high-throughput experiments data.
==========================================================

Distributional shape is often characterized by 
two features 
(1) skewness: a measure of how far the
shape of the distribution deviates from symmetry
around its location and 
(2) kurtosis: a measure of how much
weight is at the tails of the distribution relative 
to the weight around the location.

Similar to traditional moments,
the theory of L-moments
forms the basis of many statistical methods 
such as parameter estimation, hypothesis testing, 
and model selection. However, L-moments enjoy 
many theoretical and practical advantages over
traditional moments.
In this package we use L-moments ratios to
provide robust summaries of the shape of high-throughput
genomics data.

The first four L-moments L1, L2, L3 
and L4 measure location, variance,
skewness, and kurtosis of data respectively. 
Unit free measures of relative variance, skewness,
and kurtosis are defined as: 
(L-CV) L2 / L1, 
(L-skew) L3 / L2,
and 
(L-kurt) L4 / L2.

The purpose of this package is to compute
the shape (i.e. L-skew and L-kurt) statistics of each 
transcript (e.g. gene) or sample
in a high-throughput dataset (e.g. RNA-seq, microarry).
When put together these shape statistics give an overall
description of the entire high-throughput dataset.

The ability to describe the shape of high-throughput
genomics data is useful for two reasons: 
1. It provides a universal means of checking the distributional 
assumptions of statisical methods,
2. It provides a means for finding outlier genes, and
3. It provides a means for testing whether the empirical distribution
of samples differ across biological conditions.

There are four main functions in this package:

a. `fitShape()`

b. `computeDvals()`

c. `plotSO()`

d. `shapeManova()`

* Given a dataset such as a high-throughput expression matrix 
  (or just a vector of measurements) the function `fitShape()`
  will compute and return the L-CV, L-skew, and L-kurt estimates
  for each gene.

* Given the shape (i.e. L-skew and L-kurt) estimate of each gene, 
  the function `computeDvals()` computes a dissimilarity 
  distance (d-values) between each gene's shape estimate  
  and the typical gene's shape estimate.
  The d-values range from 0 to 1; where 1 is very close and 
  0 is very far.

* The function `plotSO()` shows 
  each gene's shape esitmate on a single plot. 

* The function `shapeManova()` summarizes the shape,
  (i.e. L-skew and L-kurt) estimate, of each sample 
  (column in expression matrix) and perfoms a one-way
  multivariate anova (MANOVA) to test whether sample shapes 
  are different across biological groups.

## Installation

Use [devtools](https://github.com/hadley/devtools) to install the latest
version of shape from Github:

```r
install.packages("Lmoments")
require(devtools)
install_github("kokrah/HTShape")
```

If all went well you should now be able to load shape:
```r
require(HTShape)
vignette("HTShape")
```
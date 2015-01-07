shape: Shape analysis of high-throughput experiments data.
==========================================================

Distributional shape is often characterized by 
two features 
(1) Skewness: a measure of how far the
shape of the distribution deviates from symmetry
around its location and 
(2) Kurtosis: a measure of how much
weight is at the tails of the distribution relative 
to the weight around the location.

Similar to traditional moments,
the theory of L-moments
forms the basis of many statistical methods 
such as parameter estimation, hypothesis testing, 
and model selection. However, L-moments enjoy 
many theoretical and practical advantages over
traditional moments.
In this package we focus on its ability to
provide robust statistics that summarize a given 
dataset.
The first four L-moments L1, L2, L3 
and L4 measure location, variance,
skewness, and kurtosis of data respectively. 
Unitless measures of relative variance, skewness,
and kurtosis are defined as: 
(L-CV) L2 / L1, 
(L-skew) L3 / L2,
and 
(L-kurt) L4 / L2.

The purpose of this package is to compute
the shape (L-skew, L-kurt) statistics of each 
transcript (eg. gene) 
in a high-throughput dataset (eg. RNA-seq, microarry).
Using these statistics we can find
genes within a dataset
whose sample shape is markedly different from 
the majority of gene in the same dataset. 
When put together these shape statistics give an overall
description of the entire high-throughput dataset.

The ability to describe the shape of high-throughput
genomics data is useful for two reasons: 
1. It enriches the exploratory data analysis process, and
2. It provides a means of checking the distributional 
assumptions of statisical methods.

There are three main functions in this package:

a. `fitShape()`

b. `computeDvals()`

c. `plotSO()`

* Given a dataset such as a high-throughput expression matrix 
  (or just a vector of measurements) the function `fitShape()`
  will compute and return the L-CV, L-skew, and L-kurt for each gene.

* Given the L-skew and L-kurt of each gene, the 
  function `computeDvals()` computes a dissimilarity 
  score (d-values) between each gene's (L-skew, L-shape) estimate 
  and the typical gene's (L-skew, L-shape) estimate. 
  The d-values range from 0 to 1; where 1 is very close and 
  0 is very far.

* The function `plotSO()` shows 
  each gene's (L-skew, L-shape) on a single plot (SO-plot). 


## Installation

Use [devtools](https://github.com/hadley/devtools) to install the latest
version of shape from Github:

```r
require(devtools)
install_github("kokrah/shape")
```

If all went well you should now be able to load shape:
```r
require(shape)
vignette("shape")
```
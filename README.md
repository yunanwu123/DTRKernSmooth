# Estimate and Make Inference about Optimal Treatment Regimes via Smoothed Methods

This R package probides functions to estimate the optimal treatment regime among all linear regimes via smoothed estimation methods, 
and construct element-wise confidence intervals for the optimal linear treatment regime vector, as well as the confidence interval 
for the optimal value via wild bootstrap procedures, if the population follows treatments recommended by the optimal linear regime. 

To install the package, please run the following codes in R:

```{r}
library(devtools)
install_github("yunanwu123/DTRKernSmooth")
```


## Reference

Wu, Y. and Wang, L. (2021), ***Resampling-based Confidence Intervals for Model-free Robust Inference on Optimal Treatment Regimes**, 
Biometrics, 77: 465– 476*, [doi:10.1111/biom.13337](https://doi.org/10.1111/biom.13337).

Nesterov, Y. (2007). ***Gradient methods for minimizing composite objective function. Core discussion papers, Université catholique 
de Louvain**, Center for Operations Research and Econometrics (CORE)*.

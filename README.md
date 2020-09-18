# snob
Flexible mixture models for automatic clustering

SNOB is a Matlab implementation of finite mixture models. Snob uses the minimum message length criterion to estimate the structure of the mixture model (i.e., the number of sub-populations; which sample belongs to which sub-population) and estimate all mixture model parameters. SNOB allows the user to specify the desired number of sub-populations, however if this is not specified, SNOB will automatically try to discover this information. Currently, SNOB supports mixtures of the following distributions:

-Univariate Gaussian distribution

-Multivariate Gaussian distribution

-Weibull distribution

-Exponential distribution

-Inverse Gaussian distribution

-Poisson distribution

-Negative binomial distribution

-Geometric distribution

-Multinomial distribution

-Gamma distribution

-Laplace distribution

-Single factor analysis model

-Gaussian linear regression

-Logistic regression

The program is easy to use and allows missing data - all missing data should be coded as NaN. Examples of how to use the program are provided; see data/mm_example?.m.

References:

Wallace, C. S. & Dowe, D. L. MML clustering of multi-state, Poisson, von Mises circular and Gaussian distributions. Statistics and Computing, 2000 , 10, pp. 73-83

Wallace, C. S. Intrinsic Classification of Spatially Correlated Data. The Computer Journal, 1998, 41, pp. 602-611

Wallace, C. S. Statistical and Inductive Inference by Minimum Message Length. Springer, 2005

Schmidt, D. F. & Makalic, E. Minimum Message Length Inference and Mixture Modelling of Inverse Gaussian Distributions. AI 2012: Advances in Artificial Intelligence, Springer Berlin Heidelberg, 2012, 7691, pp. 672-682

Edwards, R. T. & Dowe, D. L. Single factor analysis in MML mixture modelling. Research and Development in Knowledge Discovery and Data Mining, Second Pacific-Asia Conference (PAKDD-98), 1998, 1394

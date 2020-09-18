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

# DaME

The goal of `DaME` is to calculate a k-dimensional Dabrowska estimator using simulated or real survival data. 

# Installation

`DaME` can be installed from GitHub as follows:
```
# install.packages("devtools")
devtools::install_github("HeyItsKirill/DaME")
```

# Usage

Calculation of the k-dimensional Dabrowska estimator via `DaME` requires a data frame of k-dimensional survival data composed of survival times and censoring indicators. Users may use their own data, or simulate k-dimensional survival data using a built-in Clayton copula simulating process. 

```
x <- genClaytonk(n = 10,
                 theta = 5,
                 lambdaC = c(0.1, 0.2, 0.3))

estimate <- dabrowska(data = x,
                      k = 3)
```

# Additional Resources

For a more in-depth guide to using `DaME`, please see the attached vignette.

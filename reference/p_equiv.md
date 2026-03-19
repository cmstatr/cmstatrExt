# p-Value for one-sample equivalency

Calculates the p-Value for a one-sample acceptance test based on Vangel
(2002). This test considers the sample size of the acceptance sample
(`m`).

Two test statistics are required:

t1 = (mu - Y_min) / sigma

t2 = (mu - Y_mean) / sigma

Where:

- mu is the mean of the population

- sigma is the standard deviation of the population

- Y_min is the minimum from the acceptance sample

- Y_mean is the mean of the acceptance sample

## Usage

``` r
p_equiv(m, t1, t2)
```

## Arguments

- m:

  the size of the acceptance sample

- t1:

  the test statistic described above. May be a vector.

- t2:

  the test statistic described above. May be a vector.

## Value

a vector of p-Values of the same length as t1 and t2

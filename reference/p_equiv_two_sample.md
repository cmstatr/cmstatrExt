# p-Value for two-sample equivalency

Calculates the p-Value for a two-sample acceptance test. This test
considers the sample size of the qualification sample (`n`) and the
acceptance sample (`m`).

Two test statistics are required:

t1 = (X_mean - Y_min) / S

t2 = (X_mean - Y_mean) / S

Where:

- X_mean is the mean of the qualification sample

- S is the standard deviation of the qualification sample

- Y_min is the minimum from the acceptance sample

- Y_mean is the mean of the acceptance sample

## Usage

``` r
p_equiv_two_sample(n, m, t1, t2)
```

## Arguments

- n:

  the size of the qualification sample

- m:

  the size of the acceptance sample

- t1:

  the test statistic described above. May be a vector.

- t2:

  the test statistic described above. May be a vector.

## Value

a vector of p-Values of the same length as t1 and t2

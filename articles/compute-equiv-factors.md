# Equivalency Factor Calculator

This page provides an online calculator to determine two-sample
equivalency factors. This calculator gives the factors $k_{1}$ and
$k_{2}$ as well as determining the power of the test for detecting
reduction in mean. The basis of this method is the following paper. More
details are given at the bottom of this page.

S. Kloppenborg, “Lot acceptance testing using sample mean and extremum
with finite qualification samples,” Journal of Quality Technology, 2023.
[DOI:
10.1080/00224065.2022.2147884](https://doi.org/10.1080/00224065.2022.2147884)

*This calculator is provided as-is without any warranty. Users are
advised to review the code to verify correctness.*

## Calculator

### Input

Qualification Sample Size ($n$):  

Equivalency Sample Size ($m$):  

Significance ($\alpha$):  

**emscripten**

Downloading…

Compute Factors

  

### Factors

  

### Power for Reduction in Mean

Keep Power Curves

  

## Details

- Method Details
- Software Details

Based on a user selected qualification sample size ($n$), equivalency
sample size ($m$) and significance level ($\alpha$), the factors $k_{1}$
and $k_{2}$ are calculated. Equivalency limits are set as:

\$\$ W\_{min\\indiv} = \bar{x} - k_1 \cdot s \\ W\_{avg} = \bar{x} - k_2
\cdot s \$\$

The power of this equivalency criteria is investigated through
simulation. In this simulation, 2500 qualification samples are drawn
from a standard normal distribution ($N(\mu,\sigma)$) and equivalency
limits are computed based on each qualification sample. Next 2500
equivalency samples are drawn from a $N(\mu - \delta\sigma,\sigma)$
distribution. Each of the equivalency samples are compared against each
of the equivalency limits and the proportion of equivalency samples
rejected are reported. Thus, a total of 6,250,000 comparisons are made.
This is repeated for several values of $\delta$.

The functionality of this page is provided by the same C++ code that is
used by the `cmstatrExt` R package. This code is compiled to WebAssembly
so that it can run inside a web browser without the user installing any
special software. This software is licensed under the
[AGPL-3](https://www.r-project.org/Licenses/AGPL-3) license. Source code
is available [here](https://github.com/cmstatr/cmstatrExt).

Graphing is provided by the [Plotly JavaScript
library](https://plotly.com/javascript/), which is licensed under the
MIT license.

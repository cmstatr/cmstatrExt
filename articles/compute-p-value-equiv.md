# Equivalency p-Value Calculator

This page provides an online calculator to a p-value for two-sample
equivalency. Based on the qualification sample statistics and the
equivalency sample statistics, this calculator computes the p-value. The
basis of this method is the following paper. More details are given at
the bottom of this page.

S. Kloppenborg, “Lot acceptance testing using sample mean and extremum
with finite qualification samples,” Journal of Quality Technology, 2023.
[DOI:
10.1080/00224065.2022.2147884](https://doi.org/10.1080/00224065.2022.2147884)

*This calculator is provided as-is without any warranty. Users are
advised to review the code to verify correctness.*

## Calculator

### Input

Qualification Sample Mean (${\bar{x}}_{qual}$):  

Qualification Sample Standard Deviation ($s_{qual}$):  

Qualification Sample Size ($n$):  

Equivalency Sample Minimum ($x_{min\, equiv}$):  

Equivalency Sample Mean (${\bar{x}}_{equiv}$):  

Equivalency Sample Size ($m$):  

**emscripten**

Downloading…

Compute p-Value

  

### p-Value

  

## Details

- Method Details
- Software Details

Based on a user selected qualification sample statistics
(${\bar{x}}_{qual}$, $s_{qual}$, and $n$), equivalency sample statistics
($x_{min\, equiv}$, ${\bar{x}}_{equiv}$, and $m$), the following two
statistics are computed:

\$\$ t_1 = \frac{\bar{x}\_{qual} - x\_{min\\equiv}}{s\_{qual}} \\ t_2 =
\frac{\bar{x}\_{qual} - \bar{x}\_{equiv}}{s\_{qual}} \$\$

From these statistics, a p-Value is computed.

The functionality of this page is provided by the same C++ code that is
used by the `cmstatrExt` R package. This code is compiled to WebAssembly
so that it can run inside a web browser without the user installing any
special software. This software is licensed under the
[AGPL-3](https://www.r-project.org/Licenses/AGPL-3) license. Source code
is available [here](https://github.com/cmstatr/cmstatrExt).

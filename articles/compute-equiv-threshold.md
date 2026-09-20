# Equivalency Threshold Calculator

This page provides an online calculator to determine one- and two-sample
equivalency factors and thresholds. This calculator gives the factors
$`k_1`$ and $`k_2`$ as well as determining the power of the test for
detecting reduction in mean. The one-sample factors and thresholds are
based on the paper by Vangel below and the two-sample factors and
thresholds are based on the paper by Kloppenborg below. More details are
given at the bottom of this page.

M. Vangel, “Lot Acceptance and Compliance Testing Using the Sample Mean
and an Extremum,” Technometrics, vol. 44, no. 3. pp. 242–249, Aug-2002.

S. Kloppenborg, “Lot acceptance testing using sample mean and extremum
with finite qualification samples,” Journal of Quality Technology, 2023.
[DOI:
10.1080/00224065.2022.2147884](https://doi.org/10.1080/00224065.2022.2147884)

*This calculator is provided as-is without any warranty. Users are
advised to review the code to verify correctness.*

## Calculator

**emscripten**

Downloading…

[TABLE]

## Details

- Software Details

The functionality of this page is provided by the same C++ code that is
used by the `cmstatrExt` R package. This code is compiled to WebAssembly
so that it can run inside a web browser without the user installing any
special software. This software is licensed under the
[AGPL-3](https://www.r-project.org/Licenses/AGPL-3) license. Source code
is available [here](https://github.com/cmstatr/cmstatrExt).

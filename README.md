
<!-- README.md is generated from README.Rmd. Please edit that file -->

# What it Does

*NOTE: This package is under active development. Its API and the results
that it produces may change at any time without notice.*

`cmstatrExt` provides statistical methods intended for composite
material data that are not included in CMH-17-1G.

# Installation

This package is not yet ready to be submitted to CRAN. You can, however,
install it from GitHub using the following code:

``` r
devtools::install_github("cmstatr/cmstatrExt")
```

# Development

This package is under active development. If you have a suggestion or
question, please create an [issue on
GitHub](https://github.com/cmstatr/cmstatrExt/issues).

If you want to contribute to this package, Pull Requests are welcome!

## Debugging C++ Source

First, to add debugging symbols to the binaries, create the file
`src/Makevars` with the following content:

    ALL_CXXFLAGS = -ggdb -O0 -Wall
    ALL_CFLAGS = -ggdb -O0 -Wall
    CFLAGS =    -ggdb -O0 -Wall #-O3 -Wall -pipe -pedantic -std=gnu99

Next, launch R with a debugger. The following shows an example of doing
this, then adding a break point at line 157 of `acceptance.cpp` and
running the R function `runTests()`, which will presumably hit the break
point.

    R --debugger=gdb
    gdb> break acceptance.cpp:157
    gdb> run
    R> devtools::load_all()
    R> runTests()

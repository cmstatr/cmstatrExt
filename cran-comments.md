Submission of package `cmstatrExt` is a minor update in preparation for an
upcoming change to the Rcpp package. Calls to Rf_error have been replaced
with Rcpp::stop.

## Test environments
- win-builder (`devel`, `release`, `oldrelease`)
- local Ubuntu 24.04, R 4.5.2
- GitHub Action runners:
  - MacOS, R `release`
  - Windows, R `release`
  - Ubuntu, R `devel`
  - Ubuntu, R `release`
  - Ubuntu, R `oldrel`

## R CMD check results
There were no `ERRORs` or `WARNINGs`.

## Downstream dependencies
There are no downstream dependencies.


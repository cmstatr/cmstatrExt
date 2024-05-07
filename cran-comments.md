Re-submission of package `cmstatrExt` to address CRAN comments on submission
made 06-May-2024.

## Test environments
- win-builder (`devel`, `release`, `oldrelease`)
- local Ubuntu 22.04, R 4.4.0
- GitHub Action runners:
  - MacOS, R `release`
  - Windows, R `release`
  - Ubuntu, R `devel`
  - Ubuntu, R `release`
  - Ubuntu, R `oldrel`

## R CMD check results
There were no `ERRORs` or `WARNINGs`.

There was one expected `NOTE` on win-builder due to this being a new submission.

When built on Ubuntu, there is a `NOTE` related to the `libs` directory being
over 1 MB. I understand this to be a common `NOTE` for packages using `Rcpp`.

## Downstream dependencies
There are no downstream dependencies.

## Response to Initial Submission (2024-06-06)


>> Please add \value to .Rd files regarding exported methods and explain 
>> the functions results in the documentation. Please write about the 
>> structure of the output (class) and also what the output means. (If a 
>> function does not return a value, please document that too, e.g. 
>> \value{No return value, called for side effects} or similar)
>> Missing Rd-tags:
>>      power_sim_dual.Rd: \value

The \value tag was added to power_sim_dual.Rd. No other instances of missing
\value tag were found.

>> \dontrun{} should only be used if the example really cannot be executed 
>> (e.g. because of missing additional software, missing API keys, ...) by 
>> the user. That's why wrapping examples in \dontrun{} adds the comment 
>> ("# Not run:") as a warning for the user. Does not seem necessary. 
>> Please replace \dontrun with \donttest.
>> Please unwrap the examples if they are executable in < 5 sec, or replace 
>> dontrun{} with \donttest{}.
>> Please wrap examples that need packages in ‘Suggests’ in 
>> if(requireNamespace("pkgname")){} instead.
>> -> iso_equiv_two_sample.Rd

\dontrun replaced with \donttest since this example takes ~30 seconds to run
on my PC. The example wrapped in `if(requireNamespace("tidyverse")){}`

>> Please fix and resubmit.

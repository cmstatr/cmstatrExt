#' Example stress-strain data
#'
#' Example tension stress-strain data. This data was generated by tracing the
#' stress-strain graph for PA12 from the manuscript referenced below. The
#' non-linearity seen at low strain in the original data set was removed, then
#' the data was re-sampling to produce more tightly spaced strain values.
#' Normally-distributed error was added to the stress.
#' The code used to generate the data set can be found at
#' <https://github.com/cmstatr/cmstatrExt/blob/master/data-raw/pa12-tension.R>
#'
#' @format ## `pa12_tension`
#' A data frame with 3,212 rows and 3 columns:
#' \describe{
#'   \item{Coupon}{the coupon ID}
#'   \item{Strain}{the strain measurement `[mm/mm]`}
#'   \item{Stress}{the stress measurement `[MPa]`}
#' }
#' @source Alomarah, Amer & Ruan, Dong & Masood, S. & Gao, Zhanyuan. (2019).
#'         Compressive properties of a novel additively manufactured 3D auxetic
#'         structure. Smart Materials and Structures. 28.
#'         10.1088/1361-665X/ab0dd6.
"pa12_tension"

#' Example shear stress-shear strain data
#'
#' Example shear stress-strain data. This data was collected using a novel
#' shear test method. Coupons were made using a thermoset via
#' Fused Filament Fabrication (FFF).
#' This data requires some clean-up, including removal of the "toe", offsetting
#' the strain, and removal of the post-failure data points. These procedures
#' are demonstrated in the `stress-strain` vignette. See:
#' \code{vignette("stress-strain", package = "cmstatrExt")}
#'
#' @format ## `fff_shear`
#' A data frame with 2,316 rows and 3 columns:
#' \describe{
#'   \item{Coupon}{the coupon ID}
#'   \item{Stress}{the shear stress measurement `[psi]`}
#'   \item{Strain}{the shear strain measurement `[in/in]`}
#' }
"fff_shear"

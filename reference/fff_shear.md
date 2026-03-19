# Example shear stress-shear strain data

Example shear stress-strain data. This data was collected using a novel
shear test method. Coupons were made using a thermoset via Fused
Filament Fabrication (FFF). This data requires some clean-up, including
removal of the "toe", offsetting the strain, and removal of the
post-failure data points. These procedures are demonstrated in the
`stress-strain` vignette. See:
[`vignette("stress-strain", package = "cmstatrExt")`](https://cmstatrExt.cmstatr.net/articles/stress-strain.md)

## Usage

``` r
fff_shear
```

## Format

### `fff_shear`

A data frame with 2,316 rows and 3 columns:

- Coupon:

  the coupon ID

- Stress:

  the shear stress measurement `[psi]`

- Strain:

  the shear strain measurement `[in/in]`

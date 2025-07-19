# MoransI

Moran's I spatial autocorrelation index for images.

Moran's I is a measure of spatial autocorrelation that indicates whether nearby pixels in an image have similar values (positive autocorrelation),
dissimilar values (negative autocorrelation), or random spatial distribution.

Values range from -1 to 1:
- I > 0: Positive spatial autocorrelation (similar values cluster together)
- I = 0: Random spatial distribution
- I < 0: Negative spatial autocorrelation (dissimilar values cluster together)

## Installation

If [available in Hex](https://hex.pm/docs/publish), the package can be installed by adding `moransi` to your list of dependencies in `mix.exs`:

```elixir
def deps do
  [
    {:moransi, "~> 0.1.0"}
  ]
end
```

## Key Features

* Global Moran's I: Computes the overall spatial autocorrelation for the entire image. Provides an overall measure of spatial autocorrelation
* Local Moran's I (LISA): Calculates local indicators of spatial association for each pixel, Identifies specific locations of clusters and outliers
* Flexible connectivity: Supports both 4-connectivity (rook) and 8-connectivity (queen) neighborhoods
* Statistical testing: Includes z-scores and p-values for significance testing
* Cluster classification: Identifies spatial cluster types (High-High, Low-Low, High-Low, Low-High)
* Proper S0, S1, S2 calculations:
  * S0: Sum of all spatial weights
  * S1: Sum of squared weights considering symmetry
  * S2: Sum of squared row and column sums of the weights matrix
* Correct variance formula: Uses the standard variance formula for Moran's I under the normality assumption:
```
  Var(I) = [n((n²-3n+3)S1 - nS2 + 3S0²) - b2((n²-n)S1 - 2nS2 + 6S0²)] / [(n-1)(n-2)(n-3)S0²]
```


### Interpretation

* Moran's I Range: Values range from -1 to 1
  * Positive values: Indicate spatial clustering (similar values near each other)
  * Values near 0: Indicate random spatial distribution
  * Negative values: Indicate spatial dispersion (dissimilar values near each other)

Statistical Significance: Z-scores and p-values help determine if observed patterns are statistically significant.

### Cluster Types

* :hh - High-High (positive spatial autocorrelation)
* :ll - Low-Low (positive spatial autocorrelation)
* :hl - High-Low (negative spatial autocorrelation)
* :lh - Low-High (negative spatial autocorrelation)
* :ns - Not significant

## Main Functions

* global_morans_i/2: Computes global Moran's I with statistical significance
* local_morans_i/2: Computes local Moran's I for each pixel

## Usage

```elixir
# Compute global Moran's I
global_result = MoransI.global_morans_i(image)

# Compute local Moran's I
local_results = MoransI.local_morans_i(image)

# Use different connectivity
result = MoransI.global_morans_i(image, connectivity: :rook)
```

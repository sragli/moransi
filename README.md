# MoransI

Moran's I spatial autocorrelation index for images.

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

* Global Moran's I: Computes the overall spatial autocorrelation for the entire image
* Local Moran's I (LISA): Calculates local indicators of spatial association for each pixel
* Flexible connectivity: Supports both 4-connectivity (rook) and 8-connectivity (queen) neighborhoods
* Statistical testing: Includes z-scores and p-values for significance testing
* Cluster classification: Identifies spatial cluster types (High-High, Low-Low, High-Low, Low-High)

## Main Functions

* global_morans_i/2: Computes global Moran's I with statistical significance
* local_morans_i/2: Computes local Moran's I for each pixel
* create_test_image/2: Generates test images with different spatial patterns

## Usage

```elixir
# Create a test image
image = MoransI.create_test_image(:clustered, 5)

# Compute global Moran's I
global_result = MoransI.global_morans_i(image)

# Compute local Moran's I
local_results = MoransI.local_morans_i(image)

# Use different connectivity
result = MoransI.global_morans_i(image, connectivity: :rook)
```

## Cluster Types

* :hh - High-High (positive spatial autocorrelation)
* :ll - Low-Low (positive spatial autocorrelation)
* :hl - High-Low (negative spatial autocorrelation)
* :lh - Low-High (negative spatial autocorrelation)
* :ns - Not significant

The module handles edge cases and provides robust statistical calculations for analyzing spatial patterns in image data.
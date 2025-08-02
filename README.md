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
* Uses the standard variance formula for Moran's I under the normality assumption:

$$
  Var(I) = \frac{n((n²-3n+3)S1 - nS2 + 3S0²) - b2((n²-n)S1 - 2nS2 + 6S0²)}{(n-1)(n-2)(n-3)S0²}
$$

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

* `global_morans_i/2`: Computes global Moran's I with statistical significance
* `local_morans_i/2`: Computes local Moran's I for each pixel

## Usage

```elixir
# Basic usage
result = MoransI.global_morans_i(image)

# With options
local_results = MoransI.local_morans_i(image, 
  connectivity: :rook,
  parallel: true,
  chunk_size: 2000
)

# For very large images
huge_image_results = MoransI.local_morans_i(huge_image,
  parallel: true,
  chunk_size: 5000  # Adjust based on available memory
)
```

## How Moran's I Works

Moran's I measures the degree of spatial clustering or dispersion of these attribute values within the point cloud. Computing Moran's I for a point cloud involves assessing the spatial autocorrelation of attribute values associated with each point in the cloud.

### Steps to Compute Moran's I for a Point Cloud

* Define Attribute Values: Assign attribute values (e.g., population density, temperature, elevation) to each point in the point cloud dataset.
* Null Hypothesis: The null hypothesis for Moran's I is typically that the variable exhibits a random spatial pattern, meaning there is no spatial autocorrelation.
* Determine the Expected Value and Variance of Moran's I: Under the null hypothesis, calculate the expected value and variance of Moran's I. The expected value of Moran's I under the null hypothesis of spatial randomness is typically close to $-1/(n-1)$, where n is the number of spatial units.
* Compute Spatial Weights Matrix: Define spatial relationships between points. Common approaches include:
  * Distance-based weights (e.g., inverse distance, nearest neighbors).
  * K-nearest neighbors, distance bands, or other neighborhood definitions based on spatial proximity.
  * Create a spatial weights matrix (W) that quantifies the relationships between points based on the chosen method.
* Normalize Attribute Values: Standardize attribute values to have a mean of zero and a standard deviation of one. This normalization is crucial for Moran's I computation.
* Compute Moran's I: Use the formula for Moran's I to calculate the spatial autocorrelation:

$$
I = \frac{n}{\sum_{i=1}^{n} \sum_{j=1}^{n} w_{ij} \cdot \frac{ \sum_{i=1}^{n} \sum_{j=1}^{n} w_{ij} (x_i - \bar{x})(x_j - \bar{x}) }{ \sum_{i=1}^{n} (x_i - \bar{x})^2 }}
$$

  where:
  * n is the number of points
  * $x_i$​ and $x_j$​ are the standardized attribute values at points i and j respectively
  * $\bar{x}$ is the mean of standardized attribute values
  * $w_{ij}$​ represents the spatial weight between points i and j.

* Interpret Moran's I: Moran's I ranges from -1 (perfect dispersion) to +1 (perfect clustering). Values close to 0 indicate spatial randomness. The interpretation is context-dependent, the actual value is less important than the statistical significance of the result. Apart from statistical significance, consider the practical significance of your findings. In some cases, even a moderate level of spatial autocorrelation could have important implications for your research question.
* Test for Significance:
  * The significance of Moran's I is tested using its Z-score. If the Z-score falls within the critical range of a standard normal distribution (e.g., beyond ±$1.96 for a 95% confidence level), the null hypothesis of spatial randomness can be rejected.
  * Use of Permutation Tests: An alternative and often more robust approach is to use permutation tests. This involves randomly permuting the spatial locations of the values many times (e.g., 999 permutations) and recalculating Moran's I for each permutation. This creates a reference distribution against which the observed Moran's I can be compared.

  ### Reporting Moran's I

  Reporting Moran's I results in a scientific paper or report should be clear and comprehensive, providing all necessary details for readers to understand and evaluate your findings. Here's a guideline on how to report these results:

  * State the Moran's I Value: Begin by reporting the calculated Moran's I statistic. This value ranges from -1 to +1, where values close to +1 indicate positive spatial autocorrelation, values close to -1 indicate negative spatial autocorrelation, and values near 0 suggest no spatial autocorrelation. Example: "The calculated Moran's I for [variable of interest] was 0.45, suggesting a moderate positive spatial autocorrelation."
  * Provide the Statistical Significance: Include the p-value obtained from the significance test (e.g., permutation test) to indicate whether the observed spatial autocorrelation is statistically significant. Example: "This value was statistically significant (p = 0.01), indicating that the observed spatial pattern is unlikely to be due to random chance."
  * Report the Number of Permutations (if applicable): If a permutation test was used to assess significance, state the number of permutations. Example: "The significance of Moran's I was evaluated using 10,000 permutations."
  * Include the Expected Value of Moran's I (if relevant): Sometimes, reporting the expected value under the null hypothesis can provide context for the observed Moran's I.  Example: "Under the null hypothesis of no spatial autocorrelation, the expected Moran's I is -0.02."
  * Describe the Spatial Weights Matrix: Briefly describe the type of spatial weights matrix used, as this can affect the results. The choice depends on the nature of the spatial relationships in your data (e.g., distance-based, contiguity-based). Example: "The analysis utilized an inverse distance weighting scheme for the spatial weights matrix."
  * Provide Contextual Interpretation: Offer a brief interpretation of what the Moran's I value means in the context of your study. Example: "The positive Moran's I value indicates that areas with high values of _variable X_ are clustered together in the study region."
  * Mention Software or Tools Used: If relevant, state the software or statistical tools used to perform the Moran's I analysis. Example: "Moran's I was calculated using the PySAL library in Python."
  * Discuss Limitations or Assumptions (if necessary): If there were any notable limitations or assumptions in your spatial analysis, briefly discuss these. Example: "It should be noted that the analysis assumes stationarity in the spatial relationship across the study area."
  * Graphical Representation (optional): Sometimes, including a map or graphical representation of the spatial autocorrelation can help readers visualize the patterns. Example: "Figure 2 illustrates the spatial distribution of _variable X_, highlighting areas of significant clustering."
  
  In summary, a comprehensive report of Moran's I results includes the value of Moran's I, its statistical significance, details about the spatial weights matrix, and an interpretation of the findings in the context of the study. Providing these details ensures clarity and allows readers to fully understand the implications of your spatial analysis.
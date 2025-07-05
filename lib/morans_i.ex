defmodule MoransI do
  @moduledoc """
  Compute Moran's I spatial autocorrelation index for image data.

  Moran's I is a measure of spatial autocorrelation that indicates whether
  nearby pixels in an image have similar values (positive autocorrelation),
  dissimilar values (negative autocorrelation), or random spatial distribution.

  Values range from -1 to 1:
  - I > 0: Positive spatial autocorrelation (similar values cluster together)
  - I = 0: Random spatial distribution
  - I < 0: Negative spatial autocorrelation (dissimilar values cluster together)
  """

  @doc """
  Compute global Moran's I for an image using 8-connectivity (queen's case).

  ## Parameters
  - `image`: 2D list of numeric values representing the image
  - `options`: Keyword list of options
    - `:connectivity`: `:queen` (8-connectivity, default) or `:rook` (4-connectivity)

  ## Returns
  A map containing:
  - `:morans_i`: The Moran's I statistic
  - `:expected_i`: Expected value under null hypothesis
  - `:variance`: Variance of Moran's I
  - `:z_score`: Standardized z-score
  - `:p_value`: Approximate p-value (two-tailed)
  """
  def global_morans_i(image, options \\ []) do
    connectivity = Keyword.get(options, :connectivity, :queen)

    # Convert to flat lists with coordinates
    {values, coords} = image_to_flat_data(image)
    n = length(values)

    # Calculate mean
    mean = Enum.sum(values) / n

    # Calculate deviations from mean
    deviations = Enum.map(values, &(&1 - mean))

    # Build spatial weights matrix
    weights = build_weights_matrix(coords, connectivity)

    # Calculate Moran's I components
    numerator = calculate_numerator(deviations, weights)
    denominator = calculate_denominator(deviations)
    w_sum = weights |> List.flatten() |> Enum.sum()

    # Moran's I statistic
    morans_i = if w_sum > 0 and denominator > 0 do
      (n / w_sum) * (numerator / denominator)
    else
      0.0
    end

    # Expected value under null hypothesis
    expected_i = -1.0 / (n - 1)

    # Calculate variance and z-score
    variance = calculate_variance(n, weights)
    z_score = if variance > 0 do
      (morans_i - expected_i) / :math.sqrt(variance)
    else
      0.0
    end

    # Approximate p-value (two-tailed)
    p_value = 2 * (1 - standard_normal_cdf(abs(z_score)))

    %{
      morans_i: Float.round(morans_i, 6),
      expected_i: Float.round(expected_i, 6),
      variance: Float.round(variance, 6),
      z_score: Float.round(z_score, 6),
      p_value: Float.round(p_value, 6)
    }
  end

  @doc """
  Compute local Moran's I (LISA - Local Indicators of Spatial Association) for each pixel.

  ## Parameters
  - `image`: 2D list of numeric values representing the image
  - `options`: Keyword list of options
    - `:connectivity`: `:queen` (8-connectivity, default) or `:rook` (4-connectivity)

  ## Returns
  A 2D list of maps, each containing:
  - `:local_i`: Local Moran's I value for the pixel
  - `:z_score`: Standardized z-score
  - `:p_value`: Approximate p-value
  - `:cluster_type`: Type of spatial cluster (:hh, :ll, :hl, :lh, or :ns)
 """
  def local_morans_i(image, options \\ []) do
    connectivity = Keyword.get(options, :connectivity, :queen)

    # Convert to flat lists with coordinates
    {values, coords} = image_to_flat_data(image)
    n = length(values)

    # Calculate mean and variance
    mean = Enum.sum(values) / n
    variance = calculate_sample_variance(values, mean)

    # Calculate deviations from mean
    deviations = Enum.map(values, &(&1 - mean))

    # Build spatial weights matrix
    weights = build_weights_matrix(coords, connectivity)

    # Calculate local Moran's I for each pixel
    local_results =
      Enum.with_index(deviations)
      |> Enum.map(fn {dev_i, i} ->
        # Get neighbors for pixel i
        neighbors = Enum.at(weights, i)

        # Calculate local Moran's I
        weighted_sum =
          neighbors
          |> Enum.with_index()
          |> Enum.reduce(0, fn {w_ij, j}, acc ->
            acc + w_ij * Enum.at(deviations, j)
          end)

        local_i = if variance > 0 do
          dev_i * weighted_sum / variance
        else
          0.0
        end

        # Calculate z-score (simplified)
        expected_local = -1.0 / (n - 1)
        local_variance = calculate_local_variance(i, weights, variance)

        z_score = if local_variance > 0 do
          (local_i - expected_local) / :math.sqrt(local_variance)
        else
          0.0
        end

        # Approximate p-value
        p_value = 2 * (1 - standard_normal_cdf(abs(z_score)))

        # Determine cluster type
        value_i = Enum.at(values, i)
        neighbor_mean = calculate_neighbor_mean(i, values, weights)
        cluster_type = classify_cluster(value_i, neighbor_mean, mean, p_value < 0.05)

        %{
          local_i: Float.round(local_i, 6),
          z_score: Float.round(z_score, 6),
          p_value: Float.round(p_value, 6),
          cluster_type: cluster_type
        }
      end)

    # Convert back to 2D structure
    cols = length(hd(image))

    local_results
    |> Enum.chunk_every(cols)
  end

  @doc """
  Create a simple test image with spatial patterns for demonstration.

  ## Parameters
  - `type`: Type of pattern (`:clustered`, `:random`, `:dispersed`)
  - `size`: Size of the square image (default: 10)

  ## Returns
  A 2D list representing the test image.
  """
  def create_test_image(type \\ :clustered, size \\ 10) do
    case type do
      :clustered ->
        # Create clustered pattern
        for i <- 0..(size-1) do
          for j <- 0..(size-1) do
            cond do
              i < div(size, 2) and j < div(size, 2) -> 1
              i >= div(size, 2) and j >= div(size, 2) -> 1
              true -> 0
            end
          end
        end

      :random ->
        # Create random pattern
        :rand.seed(:exsplus, {1, 2, 3})
        for _i <- 0..(size-1) do
          for _j <- 0..(size-1) do
            :rand.uniform(2) - 1
          end
        end

      :dispersed ->
        # Create checkerboard pattern
        for i <- 0..(size-1) do
          for j <- 0..(size-1) do
            rem(i + j, 2)
          end
        end
    end
  end

  defp image_to_flat_data(image) do
    {values, coords} =
      image
      |> Enum.with_index()
      |> Enum.flat_map(fn {row, i} ->
        row
        |> Enum.with_index()
        |> Enum.map(fn {value, j} ->
          {value, {i, j}}
        end)
      end)
      |> Enum.unzip()

    {values, coords}
  end

  defp build_weights_matrix(coords, connectivity) do
    n = length(coords)

    for i <- 0..(n-1) do
      {row_i, col_i} = Enum.at(coords, i)

      for j <- 0..(n-1) do
        if i == j do
          0.0
        else
          {row_j, col_j} = Enum.at(coords, j)
          if are_neighbors({row_i, col_i}, {row_j, col_j}, connectivity) do
            1.0
          else
            0.0
          end
        end
      end
    end
  end

  defp are_neighbors({r1, c1}, {r2, c2}, connectivity) do
    row_diff = abs(r1 - r2)
    col_diff = abs(c1 - c2)

    case connectivity do
      :rook -> (row_diff == 1 and col_diff == 0) or (row_diff == 0 and col_diff == 1)
      :queen -> row_diff <= 1 and col_diff <= 1 and (row_diff + col_diff) > 0
    end
  end

  defp calculate_numerator(deviations, weights) do
    n = length(deviations)

    for i <- 0..(n-1), j <- 0..(n-1), reduce: 0.0 do
      acc ->
        w_ij = weights |> Enum.at(i) |> Enum.at(j)
        dev_i = Enum.at(deviations, i)
        dev_j = Enum.at(deviations, j)
        acc + w_ij * dev_i * dev_j
    end
  end

  defp calculate_denominator(deviations) do
    deviations
    |> Enum.map(&(&1 * &1))
    |> Enum.sum()
  end

  defp calculate_variance(n, weights) do
    # Simplified variance calculation
    w_sum = weights |> List.flatten() |> Enum.sum()

    if w_sum > 0 do
      2.0 / ((n - 1) * w_sum)
    else
      0.0
    end
  end

  defp calculate_local_variance(_i, _weights, variance) do
    # Simplified local variance calculation
    variance / 100
  end

  defp calculate_sample_variance(values, mean) do
    n = length(values)
    if n > 1 do
      sum_sq_dev = values |> Enum.map(&((&1 - mean) * (&1 - mean))) |> Enum.sum()
      sum_sq_dev / (n - 1)
    else
      0.0
    end
  end

  defp calculate_neighbor_mean(i, values, weights) do
    neighbors = Enum.at(weights, i)

    {sum, count} =
      neighbors
      |> Enum.with_index()
      |> Enum.reduce({0.0, 0}, fn {w_ij, j}, {acc_sum, acc_count} ->
        if w_ij > 0 do
          {acc_sum + Enum.at(values, j), acc_count + 1}
        else
          {acc_sum, acc_count}
        end
      end)

    if count > 0 do
      sum / count
    else
      0.0
    end
  end

  defp classify_cluster(value, neighbor_mean, global_mean, significant?) do
    if significant? do
      high_value = value > global_mean
      high_neighbors = neighbor_mean > global_mean

      case {high_value, high_neighbors} do
        {true, true} -> :hh    # High-High cluster
        {false, false} -> :ll  # Low-Low cluster
        {true, false} -> :hl   # High-Low outlier
        {false, true} -> :lh   # Low-High outlier
      end
    else
      :ns  # Not significant
    end
  end

  defp standard_normal_cdf(z) do
    # Approximation of standard normal CDF
    0.5 * (1 + erf(z / :math.sqrt(2)))
  end

  defp erf(x) do
    # Approximation of error function
    a = 0.147
    exp_term = :math.exp(-x * x * (4 / :math.pi() + a * x * x) / (1 + a * x * x))
    sign = if x >= 0, do: 1, else: -1
    sign * :math.sqrt(1 - exp_term)
  end
end

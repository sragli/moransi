defmodule MoransI do
  @moduledoc """
  Compute Moran's I spatial autocorrelation index for image data.
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

    mean = Enum.sum(values) / n

    deviations = Enum.map(values, &(&1 - mean))

    # Spatial weights matrix
    weights = build_weights_matrix(coords, connectivity)

    # Moran's I components
    numerator = calculate_numerator(deviations, weights)
    denominator = calculate_denominator(deviations)

    w_sum =
      weights
      |> List.flatten()
      |> Enum.sum()

    # Moran's I statistic
    morans_i =
      if w_sum > 0 and denominator > 0 do
        n / w_sum * (numerator / denominator)
      else
        0.0
      end

    # Expected value under null hypothesis
    expected_i = -1.0 / (n - 1)

    variance = calculate_variance(n, weights)

    z_score = z_score(variance, expected_i, morans_i)

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

    mean = Enum.sum(values) / n
    variance = calculate_sample_variance(values, mean)

    deviations = Enum.map(values, &(&1 - mean))

    # Spatial weights matrix
    weights = build_weights_matrix(coords, connectivity)

    # Local Moran's I for each pixel
    local_results =
      Enum.with_index(deviations)
      |> Enum.map(fn {dev_i, i} ->
        # Get neighbors for pixel i
        neighbors = Enum.at(weights, i)

        weighted_sum =
          neighbors
          |> Enum.with_index()
          |> Enum.reduce(0, fn {w_ij, j}, acc ->
            acc + w_ij * Enum.at(deviations, j)
          end)

        local_i =
          if variance > 0 do
            dev_i * weighted_sum / variance
          else
            0.0
          end

        expected_local = -1.0 / (n - 1)
        local_variance = calculate_local_variance(i, weights)

        z_score = z_score(local_variance, expected_local, local_i)

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

    for i <- 0..(n - 1) do
      {row_i, col_i} = Enum.at(coords, i)

      for j <- 0..(n - 1) do
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
      :queen -> row_diff <= 1 and col_diff <= 1 and row_diff + col_diff > 0
    end
  end

  defp calculate_numerator(deviations, weights) do
    n = length(deviations)

    for i <- 0..(n - 1), j <- 0..(n - 1), reduce: 0.0 do
      acc ->
        w_ij =
          weights
          |> Enum.at(i)
          |> Enum.at(j)

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

  #
  # Variance calculation using the proper mathematical formula for Moran's I variance instead
  # of the simplified approximation.
  #
  defp calculate_variance(n, weights) do
    w_sum = weights |> List.flatten() |> Enum.sum()

    if w_sum == 0 do
      0.0
    else
      # S0 (sum of all weights)
      s0 = w_sum

      # S1 (sum of squared weights, considering symmetry)
      s1 =
        for i <- 0..(n - 1), j <- 0..(n - 1), i != j, reduce: 0.0 do
          acc ->
            w_ij = weights |> Enum.at(i) |> Enum.at(j)
            w_ji = weights |> Enum.at(j) |> Enum.at(i)
            acc + (w_ij + w_ji) * (w_ij + w_ji)
        end

      s1 = s1 / 2.0

      # S2 (sum of squared row and column sums)
      row_sums = weights |> Enum.map(&Enum.sum/1)

      col_sums =
        for j <- 0..(n - 1) do
          for i <- 0..(n - 1) do
            weights |> Enum.at(i) |> Enum.at(j)
          end
          |> Enum.sum()
        end

      s2 =
        Enum.sum(Enum.map(row_sums, &(&1 * &1))) +
          Enum.sum(Enum.map(col_sums, &(&1 * &1)))

      # b2 (Kurtosis measure)
      # Assuming normal distribution (b2 = 3), but this should ideally be calculated from the actual data
      b2 = 3.0

      # Variance formula for Moran's I under normality assumption
      # Var(I) = [n((n²-3n+3)S1 - nS2 + 3S0²) - b2((n²-n)S1 - 2nS2 + 6S0²)] / [(n-1)(n-2)(n-3)S0²]
      numerator =
        n * ((n * n - 3 * n + 3) * s1 - n * s2 + 3 * s0 * s0) -
          b2 * ((n * n - n) * s1 - 2 * n * s2 + 6 * s0 * s0)

      denominator = (n - 1) * (n - 2) * (n - 3) * s0 * s0

      if denominator != 0 do
        numerator / denominator
      else
        # Fallback to simpler approximation
        2.0 / ((n - 1) * s0)
      end
    end
  end

  #
  # Local variance calculation using the actual mathematical formula for Local Indicators of Spatial
  # Association (LISA) variance, based on the theoretical derivation of local Moran's I variance
  # under the null hypothesis of no spatial autocorrelation.
  #
  defp calculate_local_variance(i, weights) do
    n = length(weights)

    # Weights for observation i
    w_i = Enum.at(weights, i)

    # Sum of weights for observation i
    w_i_sum = Enum.sum(w_i)

    if w_i_sum == 0 do
      0.0
    else
      # Sum of squared weights for observation i
      w_i_sq_sum = w_i |> Enum.map(&(&1 * &1)) |> Enum.sum()

      # b2 (Kurtosis) - assuming normal distribution for now
      b2 = 3.0

      # Local variance formula components
      # E[Ii²] under null hypothesis
      term1 = w_i_sq_sum * (n - b2) / (n - 1)
      term2 = w_i_sum * w_i_sum * (2 * b2 - n) / ((n - 1) * (n - 2))

      # Additional correction term for finite sample
      correction =
        if n > 3 do
          w_i_sum * w_i_sum / ((n - 1) * (n - 1))
        else
          0.0
        end

      # Combine terms for local variance
      local_var = term1 + term2 - correction

      # Ensure non-negative variance
      max(0.0, local_var)
    end
  end

  defp calculate_sample_variance(values, mean) when length(values) > 1 do
    sum_sq_dev =
      values
      |> Enum.map(&((&1 - mean) * (&1 - mean)))
      |> Enum.sum()

    sum_sq_dev / (length(values) - 1)
  end

  defp calculate_sample_variance(_values, _mean) do
    0.0
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

  defp classify_cluster(_value, _neighbor_mean, _global_mean, _significant = false) do
    # Not significant
    :ns
  end

  defp classify_cluster(value, neighbor_mean, global_mean, _significant = true) do
    high_value = value > global_mean
    high_neighbors = neighbor_mean > global_mean

    case {high_value, high_neighbors} do
      # High-High cluster
      {true, true} -> :hh
      # Low-Low cluster
      {false, false} -> :ll
      # High-Low outlier
      {true, false} -> :hl
      # Low-High outlier
      {false, true} -> :lh
    end
  end

  defp z_score(variance, expected_i, i) when variance > 0,
    do: (i - expected_i) / :math.sqrt(variance)

  defp z_score(_variance, _expected_i, _i), do: 0.0

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

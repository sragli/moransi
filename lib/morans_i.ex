defmodule MoransI do
  @moduledoc """
  Moran's I spatial autocorrelation index for image data.
  """

  @doc """
  Compute global Moran's I for an image using 8-connectivity (queen's case).

  ## Parameters
  - `image`: 2D list of numeric values representing the image
  - `options`: Keyword list of options
    - `:connectivity`: `:queen` (8-connectivity, default) or `:rook` (4-connectivity)
    - `:parallel`: boolean, whether to use parallel processing (default: true)

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

    # Convert to efficient data structures
    {values_array, coords_map, rows, cols} = prepare_image_data(image)
    n = :array.size(values_array)

    values_list = :array.to_list(values_array)
    mean = Enum.sum(values_list) / n
    deviations_array = :array.map(fn _i, val -> val - mean end, values_array)

    # Neighbor mapping (sparse representation)
    neighbors_map = build_neighbor_map(coords_map, rows, cols, connectivity)

    {numerator, w_sum} = calculate_global_components(deviations_array, neighbors_map)
    denominator = calculate_denominator(deviations_array)

    morans_i =
      if w_sum > 0 and denominator > 0 do
        n / w_sum * (numerator / denominator)
      else
        0.0
      end

    expected_i = -1.0 / (n - 1)
    variance = calculate_variance(n, neighbors_map, w_sum)

    z_score = z_score(variance, expected_i, morans_i)
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
    - `:parallel`: boolean, whether to use parallel processing (default: true)
    - `:chunk_size`: integer, size of chunks for parallel processing (default: 1000)

  ## Returns
  A 2D list of maps, each containing:
  - `:local_i`: Local Moran's I value for the pixel
  - `:z_score`: Standardized z-score
  - `:p_value`: Approximate p-value
  - `:cluster_type`: Type of spatial cluster (:hh, :ll, :hl, :lh, or :ns)
  """
  def local_morans_i(image, options \\ []) do
    connectivity = Keyword.get(options, :connectivity, :queen)
    parallel = Keyword.get(options, :parallel, true)
    chunk_size = Keyword.get(options, :chunk_size, 1000)

    # Convert to efficient data structures
    {values_array, coords_map, rows, cols} = prepare_image_data(image)
    n = :array.size(values_array)

    values_list = :array.to_list(values_array)
    mean = Enum.sum(values_list) / n
    variance = calculate_sample_variance(values_list, mean)

    deviations_array = :array.map(fn _i, val -> val - mean end, values_array)

    neighbors_map = build_neighbor_map(coords_map, rows, cols, connectivity)

    # Calculate local Moran's I for each pixel
    local_results =
      if parallel and n > chunk_size do
        calculate_local_parallel(
          values_array,
          deviations_array,
          neighbors_map,
          mean,
          variance,
          n,
          chunk_size
        )
      else
        calculate_local_sequential(
          values_array,
          deviations_array,
          neighbors_map,
          mean,
          variance,
          n
        )
      end

    # Convert back to 2D structure
    local_results
    |> Enum.chunk_every(cols)
  end

  defp prepare_image_data(image) do
    rows = length(image)
    cols = length(hd(image))

    # Flat arrays and coordinate mappings
    {values, coords} =
      image
      |> Enum.with_index()
      |> Enum.flat_map(fn {row, i} ->
        row
        |> Enum.with_index()
        |> Enum.map(fn {value, j} -> {value, {i, j}} end)
      end)
      |> Enum.unzip()

    values_array = :array.from_list(values)

    coords_map =
      coords
      |> Enum.with_index()
      |> Map.new(fn {coord, idx} -> {idx, coord} end)

    {values_array, coords_map, rows, cols}
  end

  # Build sparse neighbor representation - O(n) instead of O(n²)
  defp build_neighbor_map(coords_map, rows, cols, connectivity) do
    coords_map
    |> Enum.map(fn {idx, {row, col}} ->
      neighbors = get_neighbors(row, col, rows, cols, connectivity)

      neighbor_indices =
        neighbors
        |> Enum.map(fn coord ->
          Enum.find_value(coords_map, fn {i, c} -> if c == coord, do: i end)
        end)
        |> Enum.reject(&is_nil/1)

      {idx, neighbor_indices}
    end)
    |> Map.new()
  end

  # Direct coordinate calculation - much faster than distance checking
  defp get_neighbors(row, col, max_rows, max_cols, :queen) do
    for r <- max(0, row - 1)..min(max_rows - 1, row + 1),
        c <- max(0, col - 1)..min(max_cols - 1, col + 1),
        {r, c} != {row, col},
        do: {r, c}
  end

  defp get_neighbors(row, col, max_rows, max_cols, :rook) do
    [{row - 1, col}, {row + 1, col}, {row, col - 1}, {row, col + 1}]
    |> Enum.filter(fn {r, c} ->
      r >= 0 and r < max_rows and c >= 0 and c < max_cols
    end)
  end

  # Global calculation using sparse neighbors
  defp calculate_global_components(deviations_array, neighbors_map) do
    neighbors_map
    |> Enum.reduce({0.0, 0}, fn {i, neighbors}, {num_acc, w_acc} ->
      dev_i = :array.get(i, deviations_array)

      {neighbor_sum, neighbor_count} =
        neighbors
        |> Enum.reduce({0.0, 0}, fn j, {sum, count} ->
          dev_j = :array.get(j, deviations_array)
          {sum + dev_j, count + 1}
        end)

      numerator_contrib = dev_i * neighbor_sum
      {num_acc + numerator_contrib, w_acc + neighbor_count}
    end)
  end

  defp calculate_denominator(deviations_array) do
    deviations_array
    |> :array.to_list()
    |> Enum.reduce(0.0, fn dev, acc -> acc + dev * dev end)
  end

  defp calculate_variance(_n, _neighbors_map, 0), do: 0.0

  defp calculate_variance(n, neighbors_map, w_sum) do
    s0 = w_sum

    # S1: sum of squared weights (each neighbor relationship has weight 1)
    # Each edge counted twice in undirected graph
    s1 = 2 * w_sum

    # S2: sum of squared row and column sums
    row_sums = neighbors_map |> Enum.map(fn {_i, neighbors} -> length(neighbors) end)
    # Symmetric matrix
    s2 = 2 * Enum.sum(Enum.map(row_sums, &(&1 * &1)))

    # Assuming normal distribution
    b2 = 3.0

    numerator =
      n * ((n * n - 3 * n + 3) * s1 - n * s2 + 3 * s0 * s0) -
        b2 * ((n * n - n) * s1 - 2 * n * s2 + 6 * s0 * s0)

    denominator = (n - 1) * (n - 2) * (n - 3) * s0 * s0

    if denominator != 0, do: numerator / denominator, else: 2.0 / ((n - 1) * s0)
  end

  # Sequential local calculation
  defp calculate_local_sequential(values_array, deviations_array, neighbors_map, mean, variance, n) do
    0..(n - 1)
    |> Enum.map(fn i ->
      calculate_local_i(i, values_array, deviations_array, neighbors_map, mean, variance, n)
    end)
  end

  # Parallel local calculation
  defp calculate_local_parallel(
         values_array,
         deviations_array,
         neighbors_map,
         mean,
         variance,
         n,
         chunk_size
       ) do
    0..(n - 1)
    |> Enum.chunk_every(chunk_size)
    |> Task.async_stream(
      fn chunk ->
        Enum.map(chunk, fn i ->
          calculate_local_i(i, values_array, deviations_array, neighbors_map, mean, variance, n)
        end)
      end,
      max_concurrency: System.schedulers_online()
    )
    |> Enum.flat_map(fn {:ok, results} -> results end)
  end

  # Calculate local Moran's I for a single pixel
  defp calculate_local_i(i, values_array, deviations_array, neighbors_map, mean, variance, n) do
    dev_i = :array.get(i, deviations_array)
    neighbors = Map.get(neighbors_map, i, [])

    weighted_sum =
      neighbors
      |> Enum.reduce(0.0, fn j, acc ->
        dev_j = :array.get(j, deviations_array)
        # Weight is 1.0 for all neighbors
        acc + dev_j
      end)

    local_i =
      if variance > 0 do
        dev_i * weighted_sum / variance
      else
        0.0
      end

    expected_local = -1.0 / (n - 1)
    local_variance = calculate_local_variance(neighbors, n)

    z_score = z_score(local_variance, expected_local, local_i)
    p_value = 2 * (1 - standard_normal_cdf(abs(z_score)))

    # Determine cluster type
    value_i = :array.get(i, values_array)
    neighbor_mean = calculate_neighbor_mean(neighbors, values_array)
    cluster_type = classify_cluster(value_i, neighbor_mean, mean, p_value < 0.05)

    %{
      local_i: Float.round(local_i, 6),
      z_score: Float.round(z_score, 6),
      p_value: Float.round(p_value, 6),
      cluster_type: cluster_type
    }
  end

  defp calculate_local_variance([], _n), do: 0.0

  defp calculate_local_variance(neighbors, n) do
    w_i_sum = length(neighbors)

    # Each weight is 1.0, so sum of squares equals sum
    w_i_sq_sum = w_i_sum
    b2 = 3.0

    term1 = w_i_sq_sum * (n - b2) / (n - 1)
    term2 = w_i_sum * w_i_sum * (2 * b2 - n) / ((n - 1) * (n - 2))

    correction =
      if n > 3 do
        w_i_sum * w_i_sum / ((n - 1) * (n - 1))
      else
        0.0
      end

    local_var = term1 + term2 - correction
    max(0.0, local_var)
  end

  defp calculate_sample_variance(values, mean) when length(values) > 1 do
    sum_sq_dev =
      values
      |> Enum.reduce(0.0, fn val, acc ->
        dev = val - mean
        acc + dev * dev
      end)

    sum_sq_dev / (length(values) - 1)
  end

  defp calculate_sample_variance(_values, _mean), do: 0.0

  defp calculate_neighbor_mean([], _values_array), do: 0.0

  defp calculate_neighbor_mean(neighbors, values_array) do
    sum = neighbors |> Enum.reduce(0.0, fn j, acc -> acc + :array.get(j, values_array) end)
    sum / length(neighbors)
  end

  defp classify_cluster(_value, _neighbor_mean, _global_mean, false), do: :ns

  defp classify_cluster(value, neighbor_mean, global_mean, true) do
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

  defp z_score(variance, expected_i, i) when variance > 0 do
    (i - expected_i) / :math.sqrt(variance)
  end

  defp z_score(_variance, _expected_i, _i), do: 0.0

  defp standard_normal_cdf(z) do
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

defmodule MoransI.Utils.Renormalization do
  @moduledoc """
  Coarse-graining using the sum of pixels in partition_size blocks, mean thresholding and output to a binary matrix.
  """

  @spec coarse_grain(Nx.Tensor.t(), integer()) :: [list()]
  def coarse_grain(image_matrix, partition_size) do
    blocks =
      image_matrix
      |> partition(partition_size)
      |> Enum.map(fn m -> Nx.sum(m) end)

    min_max = Enum.min_max(blocks)
    threshold = (Nx.to_number(elem(min_max, 1)) - Nx.to_number(elem(min_max, 0))) / 2

    new_width = div(Nx.axis_size(image_matrix, 0), partition_size)

    blocks
    |> Enum.map(fn v -> if Nx.to_number(v) < threshold, do: 0, else: 1 end)
    |> Enum.chunk_every(new_width)
  end

  defp partition(tensor, partition_size) do
    case Nx.shape(tensor) do
      {} ->
        []

      {0} ->
        []

      {0, _} ->
        []

      {_, 0} ->
        []

      {rows, cols} ->
        num_row_blocks = div(rows, partition_size)
        num_col_blocks = div(cols, partition_size)

        for row_block <- 0..(num_row_blocks - 1),
            col_block <- 0..(num_col_blocks - 1) do
          start_row = row_block * partition_size
          start_col = col_block * partition_size
          Nx.slice(tensor, [start_row, start_col], [partition_size, partition_size])
        end
    end
  end
end

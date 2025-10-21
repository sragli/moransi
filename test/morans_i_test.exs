defmodule MoransITest do
  use ExUnit.Case

  test "calculates correct global result" do
    image = create_test_image(:clustered, 5)

    assert %{
             morans_i: 0.386663,
             expected_i: -0.041667,
             variance: 0.015006,
             z_score: 3.496612,
             p_value: 4.59e-4
           } = MoransI.global_morans_i(image)
  end

  test "calculates correct local result" do
    image = create_test_image(:clustered, 5)

    assert [
             [
               %{cluster_type: :ns, local_i: 2.658462, p_value: 0.082845, z_score: 1.734064},
               %{cluster_type: :ns, local_i: 0.738462, p_value: 0.684244, z_score: 0.406702},
               %{cluster_type: :ns, local_i: 1.2, p_value: 0.517478, z_score: 0.647314},
               %{cluster_type: :ll, local_i: 5.2, p_value: 0.006206, z_score: 2.732622},
               %{cluster_type: :ll, local_i: 3.12, p_value: 0.042193, z_score: 2.030471}
             ],
             [
               %{cluster_type: :ns, local_i: 0.738462, p_value: 0.684244, z_score: 0.406702},
               %{cluster_type: :ns, local_i: -0.295385, p_value: 0.909835, z_score: -0.113247},
               %{cluster_type: :ns, local_i: 0.32, p_value: 0.871755, z_score: 0.161431},
               %{cluster_type: :ns, local_i: 2.32, p_value: 0.291905, z_score: 1.054134},
               %{cluster_type: :ns, local_i: 1.2, p_value: 0.517478, z_score: 0.647314}
             ],
             [
               %{cluster_type: :ns, local_i: 1.2, p_value: 0.517478, z_score: 0.647314},
               %{cluster_type: :ns, local_i: 0.32, p_value: 0.871755, z_score: 0.161431},
               %{cluster_type: :ns, local_i: -0.295385, p_value: 0.909835, z_score: -0.113247},
               %{cluster_type: :ns, local_i: 1.550769, p_value: 0.477275, z_score: 0.710787},
               %{cluster_type: :ns, local_i: 0.738462, p_value: 0.684244, z_score: 0.406702}
             ],
             [
               %{cluster_type: :ll, local_i: 5.2, p_value: 0.006206, z_score: 2.732622},
               %{cluster_type: :ns, local_i: 2.32, p_value: 0.291905, z_score: 1.054134},
               %{cluster_type: :ns, local_i: 1.550769, p_value: 0.477275, z_score: 0.710787},
               %{cluster_type: :hh, local_i: 7.089231, p_value: 0.001428, z_score: 3.18289},
               %{cluster_type: :hh, local_i: 4.430769, p_value: 0.019601, z_score: 2.331601}
             ],
             [
               %{cluster_type: :ll, local_i: 3.12, p_value: 0.042193, z_score: 2.030471},
               %{cluster_type: :ns, local_i: 1.2, p_value: 0.517478, z_score: 0.647314},
               %{cluster_type: :ns, local_i: 0.738462, p_value: 0.684244, z_score: 0.406702},
               %{cluster_type: :hh, local_i: 4.430769, p_value: 0.019601, z_score: 2.331601},
               %{cluster_type: :ns, local_i: 2.658462, p_value: 0.082845, z_score: 1.734064}
             ]
           ] = MoransI.local_morans_i(image)
  end

  test "global_morans_i returns zero for uniform image and valid keys" do
    image = for _i <- 1..3, do: for(_j <- 1..3, do: 5)

    res = MoransI.global_morans_i(image)

    assert is_map(res)
    assert Map.has_key?(res, :morans_i)
    assert Map.has_key?(res, :p_value)
    assert res.morans_i == 0.0
    assert res.p_value >= 0.0 and res.p_value <= 1.0
  end

  test "global_morans_i raises on 1x1 image" do
    image = [[1]]

    assert_raise ArithmeticError, fn -> MoransI.global_morans_i(image) end
  end

  test "local_morans_i parallel and sequential produce same shape and valid entries" do
    image = for i <- 0..2, do: for(j <- 0..2, do: if(i == 1 and j == 1, do: 10, else: 0))

    seq = MoransI.local_morans_i(image, parallel: false)
    par = MoransI.local_morans_i(image, parallel: true, chunk_size: 2)

    assert length(seq) == length(par)
    assert length(hd(seq)) == length(hd(par))

    Enum.each(List.flatten(par), fn m ->
      assert is_map(m)
      assert Map.has_key?(m, :local_i)
      assert Map.has_key?(m, :cluster_type)
    end)
  end

  test "rook vs queen connectivity return morans_i floats" do
    image = for i <- 0..2, do: for(j <- 0..2, do: if(i + j < 2, do: 1, else: 0))

    r = MoransI.global_morans_i(image, connectivity: :rook)
    q = MoransI.global_morans_i(image, connectivity: :queen)

    assert is_float(r.morans_i)
    assert is_float(q.morans_i)
  end

  #
  # Create a simple test image with spatial patterns for demonstration.
  #
  ### Parameters
  # - Type of pattern (`:clustered`, `:random`, `:dispersed`)
  # - `size`: Size of the square image
  #
  ### Returns
  # A 2D list representing the test image.
  #
  defp create_test_image(:clustered, size) do
    for i <- 0..(size - 1) do
      for j <- 0..(size - 1) do
        cond do
          i < div(size, 2) and j < div(size, 2) -> 1
          i >= div(size, 2) and j >= div(size, 2) -> 1
          true -> 0
        end
      end
    end
  end

  defp create_test_image(:random, size) do
    :rand.seed(:exsplus, {1, 2, 3})

    for _i <- 0..(size - 1) do
      for _j <- 0..(size - 1) do
        :rand.uniform(2) - 1
      end
    end
  end

  defp create_test_image(:dispersed, size) do
    for i <- 0..(size - 1) do
      for j <- 0..(size - 1) do
        rem(i + j, 2)
      end
    end
  end
end

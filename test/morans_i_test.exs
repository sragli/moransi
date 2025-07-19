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

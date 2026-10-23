defmodule MoransITest do
  use ExUnit.Case

  test "calculates correct global result" do
    image = create_test_image(:clustered, 5)

    assert %{
             morans_i: 0.386663,
             expected_i: -0.041667,
             variance: 0.010738,
             z_score: 4.133454,
             p_value: 0.000035
           } = MoransI.global_morans_i(image)
  end

  test "calculates correct local result" do
    image = create_test_image(:clustered, 5)
    result = MoransI.local_morans_i(image)

    assert length(result) == 5
    assert length(hd(result)) == 5

    assert %{cluster_type: :ns, local_i: 2.769231, p_value: 0.073073, z_score: 1.791902} =
             result |> Enum.at(0) |> Enum.at(0)

    assert %{cluster_type: :hh, local_i: 7.384615, p_value: 0.000714, z_score: 3.377546} =
             result |> Enum.at(3) |> Enum.at(3)
  end

  test "local_morans_i uses population variance for scaling" do
    result = MoransI.local_morans_i([[1, 0], [0, 0]], parallel: false)

    assert %{local_i: -3.0} = result |> Enum.at(0) |> Enum.at(0)
  end

  test "local_morans_i handles two-cell images" do
    result = MoransI.local_morans_i([[1, 0]], parallel: false)

    assert [[first, second]] = result
    assert first.local_i == -1.0
    assert second.local_i == -1.0
    assert first.p_value == 1.0
    assert second.p_value == 1.0
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

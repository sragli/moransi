defmodule MoransITest do
  use ExUnit.Case

  test "calculates correct global result" do
    image = MoransI.create_test_image(:clustered, 5)

    assert %{
              morans_i: 0.386663,
              expected_i: -0.041667,
              variance: 5.79e-4,
              z_score: 17.805334,
              p_value: 0.0
            } = MoransI.global_morans_i(image)
  end

  test "calculates correct local result" do
    image = MoransI.create_test_image(:clustered, 5)

    assert [
              [
                %{
                  z_score: 52.953871,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 2.658462
                },
                %{
                  z_score: 15.299573,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 0.738462
                },
                %{
                  z_score: 24.351087,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 1.2
                },
                %{
                  z_score: 102.797541,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 5.2
                },
                %{
                  z_score: 62.005385,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 3.12
                }
              ],
              [
                %{
                  z_score: 15.299573,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 0.738462
                },
                %{
                  z_score: -4.975818,
                  p_value: 1.0e-6,
                  cluster_type: :hl,
                  local_i: -0.295385
                },
                %{
                  z_score: 7.092867,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 0.32
                },
                %{
                  z_score: 46.316094,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 2.32
                },
                %{
                  z_score: 24.351087,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 1.2
                }
              ],
              [
                %{
                  z_score: 24.351087,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 1.2
                },
                %{
                  z_score: 7.092867,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 0.32
                },
                %{
                  z_score: -4.975818,
                  p_value: 1.0e-6,
                  cluster_type: :hl,
                  local_i: -0.295385
                },
                %{
                  z_score: 31.230237,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 1.550769
                },
                %{
                  z_score: 15.299573,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 0.738462
                }
              ],
              [
                %{
                  z_score: 102.797541,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 5.2
                },
                %{
                  z_score: 46.316094,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 2.32
                },
                %{
                  z_score: 31.230237,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 1.550769
                },
                %{
                  z_score: 139.848405,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 7.089231
                },
                %{
                  z_score: 87.711684,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 4.430769
                }
              ],
              [
                %{
                  z_score: 62.005385,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 3.12
                },
                %{
                  z_score: 24.351087,
                  p_value: 0.0,
                  cluster_type: :ll,
                  local_i: 1.2
                },
                %{
                  z_score: 15.299573,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 0.738462
                },
                %{
                  z_score: 87.711684,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 4.430769
                },
                %{
                  z_score: 52.953871,
                  p_value: 0.0,
                  cluster_type: :hh,
                  local_i: 2.658462
                }
              ]
            ] = MoransI.local_morans_i(image)
  end
end

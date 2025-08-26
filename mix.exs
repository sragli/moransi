defmodule MoransI.MixProject do
  use Mix.Project

  def project do
    [
      app: :moransi,
      version: "0.1.1",
      elixir: "~> 1.17",
      start_permanent: Mix.env() == :prod,
      description: description(),
      package: package(),
      deps: deps(),
      name: "MoransI",
      source_url: "https://github.com/sragli/moransi",
      docs: docs()
    ]
  end

  def application do
    [
      extra_applications: [:logger]
    ]
  end

  defp description() do
    "Elixir module that implements Moran's I spatial autocorrelation index for images."
  end

  defp package() do
    [
      files: ~w(lib .formatter.exs mix.exs README.md LICENSE CHANGELOG),
      licenses: ["Apache-2.0"],
      links: %{"GitHub" => "https://github.com/sragli/moransi"}
    ]
  end

  defp docs() do
    [
      main: "MoransI",
      extras: ["README.md", "LICENSE", "examples.livemd", "CHANGELOG"]
    ]
  end

  defp deps do
    [
      {:ex_doc, ">= 0.0.0", only: :dev, runtime: false}
    ]
  end
end

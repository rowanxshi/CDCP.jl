using Documenter, CDCP

makedocs(
	sitename="Combinatorial Discrete Choice Problems",
	pages=[
		"Introduction" => "index.md",
		"Manual" => "manual.md",
		"Internals" => "internals.md",
		"Backwards compatibility" => "compat.md",
	],
	# remotes = nothing
)

deploydocs(repo = "github.com/rowanxshi/CDCP.jl.git")

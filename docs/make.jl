using Documenter, CDCP

makedocs(
	sitename="Combinatorial Discrete Choice Problems",
	pages= [
		"Introduction" => "index.md",
		"Manual" => "manual.md"
	]
)

deploydocs( repo = "github.com/rowanxshi/CDCP.jl.git")

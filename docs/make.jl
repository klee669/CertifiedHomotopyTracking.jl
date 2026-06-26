using Documenter
using DocumenterCitations
using CertifiedHomotopyTracking

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"), style=:alpha)

makedocs(
    sitename = "CertifiedHomotopyTracking.jl",
    modules = [CertifiedHomotopyTracking],
    plugins = [bib],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Homotopies and Systems" => "homotopies-and-systems.md",
        "Krawczyk Test" => "krawczyk-test.md",
        "A Posteriori Path Certification" => [
            "A Posteriori Certification" => "posteriori/certification.md",
            "Numerical Trace" => "posteriori/low-level.md",
            "Visualization" => "posteriori/visualization.md",
        ],
        "Tracking" => [
            "Path Tracking" => "tracking/path-tracking.md",
            "Tracking Results" => "tracking/results.md",
        ],
        "Monodromy" => [
            "Solving Monodromy" => "monodromy/solving.md",
            "Homotopy Graph" => "monodromy/graph.md",
            "Results" => "monodromy/results.md",
            "GAP Computations" => "monodromy/gap.md",
        ],
        "Variety Approximations" => "variety-approximations.md",
        "Visualization" => "visualization.md",
    ],
    checkdocs = :none,
)

deploydocs(
    repo = "github.com/klee669/CertifiedHomotopyTracking.jl.git",
    devbranch = "dev/v1.0.0",
)

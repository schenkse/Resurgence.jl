using Documenter
using DocumenterCitations
using Resurgence

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

makedocs(;
    modules  = [Resurgence],
    sitename = "Resurgence.jl",
    authors  = "Sebastian Schenk",
    plugins  = [bib],
    format   = Documenter.HTML(;
        prettyurls   = get(ENV, "CI", "false") == "true",
        canonical    = "https://schenkse.github.io/Resurgence.jl",
        edit_link    = "dev",
        assets       = String[],
        mathengine   = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Tutorials" => [
            "Stieltjes / Euler series"      => "tutorials/stieltjes.md",
            "Lateral and median sums"       => "tutorials/lateral_sums.md",
            "Stokes / large-order diagnostics" => "tutorials/stokes_diagnostics.md",
            "Trans-series"                  => "tutorials/transseries.md",
        ],
        "Methods guide" => "methods_guide.md",
        "API reference" => [
            "Sequence acceleration" => "api/series_acceleration.md",
            "Padé"                  => "api/pade.md",
            "Borel and conformal"   => "api/borel.md",
            "Borel–Padé"            => "api/borel_pade.md",
            "Meijer-G"              => "api/meijerg.md",
            "Truncation"            => "api/truncation.md",
            "Stokes diagnostics"    => "api/stokes.md",
            "Trans-series"          => "api/transseries.md",
            "Unified `resum` API"   => "api/api.md",
        ],
        "Roadmap"    => "roadmap.md",
        "References" => "references.md",
    ],
    checkdocs = :exports,
    warnonly  = [:missing_docs, :cross_references],
)

deploydocs(;
    repo      = "github.com/schenkse/Resurgence.jl.git",
    devbranch = "dev",
    push_preview = true,
)

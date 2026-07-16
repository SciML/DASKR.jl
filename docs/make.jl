using DASKR
using Documenter

makedocs(;
    modules = [DASKR],
    checkdocs = :exports,
    sitename = "DASKR.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://docs.sciml.ai/DASKR/stable/",
        edit_link = "master",
    ),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/SciML/DASKR.jl",
    devbranch = "master",
)

# Build the manual:   julia --project=docs docs/make.jl
# (first time: julia --project=docs -e 'using Pkg; Pkg.instantiate()')
# The examples are not run: the example pages show the stored output (examples/output/<name>/output.txt) and the
# figures of the last run, which are copied to docs/src/assets/figures/<name>/.
using Documenter, DocumenterCitations, Lightcurves

const ROOT = dirname(@__DIR__)
for name in readdir(joinpath(ROOT, "examples", "output"))
    src = joinpath(ROOT, "examples", "output", name)
    isdir(src) || continue
    dst = mkpath(joinpath(@__DIR__, "src", "assets", "figures", name))
    for f in readdir(src)
        endswith(f, ".png") && cp(joinpath(src, f), joinpath(dst, f); force=true)
    end
end
example_output(name) = read(joinpath(ROOT, "examples", "output", name, "output.txt"), String)

# example pages: template docs/src/examples/<name>.md.in with the placeholder {{OUTPUT}}
for f in readdir(joinpath(@__DIR__, "src", "examples"))
    endswith(f, ".md.in") || continue
    name = replace(f, ".md.in" => "")
    txt = replace(read(joinpath(@__DIR__, "src", "examples", f), String), "{{OUTPUT}}" => example_output(name))
    write(joinpath(@__DIR__, "src", "examples", name * ".md"), txt)
end

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)
makedocs(; sitename="Lightcurves.jl", modules=[Lightcurves], authors="Mike Kretlow",
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true", edit_link="main", size_threshold=nothing),
         pages=["Home" => "index.md", "Methods" => "methods.md", "Data access" => "data.md",
                "Examples" => ["(1486) Marilyn (ALCDEF)" => "examples/1486_marilyn_alcdef.md"],
                "API" => "api.md", "References" => "references.md"],
         plugins=[bib], checkdocs=:exports, warnonly=[:missing_docs])

# deploydocs(; repo="github.com/mkretlow/Lightcurves.jl.git", devbranch="main")   # when the repository is public

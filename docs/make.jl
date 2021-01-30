

using Documenter, Mag2Dpoly

#modules = [EikonalSolvers],
makedocs(modules = [Mag2Dpoly],
         sitename="Mag2Dpoly.jl",
         authors = "Andrea Zunino, Alessandro Ghirotto",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true"),
         pages = [
             "Home" => "index.md",
         ]
         )

deploydocs(
    repo="github.com/inverseproblem/Mag2Dpoly.jl.git",
    devbranch = "main"
)

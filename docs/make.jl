

using Documenter, Mag2Dpoly

makedocs(sitename="Mag2Dpoly.jl",
         modules = [Mag2Dpoly],
         authors = "Andrea Zunino, Alessandro Ghirotto",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true"),
         pages = [
             "Home" => "index.md",
         ]
         )

deploydocs(
    repo="github.com/inverseproblem/Mag2Dpoly.jl.git",
    devbranch = "main",
    branch = "gh-pages" 
)

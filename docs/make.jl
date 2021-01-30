

using Documenter, Mag2Dpoly

if haskey(ENV, "GITHUB_ACTIONS")
         println("################################################")
         println("###  haskey(ENV, GITHUB_ACTIONS) == true ")
         println("################################################")
end

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



using Documenter, Mag2Dpoly, PyPlot

if haskey(ENV, "GITHUB_ACTIONS")
         println("################################################")
         println("###  haskey(ENV, GITHUB_ACTIONS) == true ")
         println("################################################")
end

adds = Documenter.auto_detect_deploy_system();
println("################################")
println("### Auto detect deploy system: ###")
println("###      $adds             ###")
println("################################")


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

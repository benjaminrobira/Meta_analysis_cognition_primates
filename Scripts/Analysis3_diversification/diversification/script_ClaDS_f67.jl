##/usr/local/src/julia-1.4.2/bin/julia

cd("/users/biodiv/bperez/appli/ClaDS_Julia/")

import Pkg
Pkg.add("RCall") 
Pkg.add("RDatasets")
Pkg.add("StatsBase")
Pkg.add("Distributions")
Pkg.add("BenchmarkTools")
Pkg.add("JLD2")
Pkg.add("DataFrames")
Pkg.add("DataFramesMeta")
Pkg.add("CSV")

using DataFrames
using RCall
using DataFramesMeta
using CSV

cd("/users/biodiv/bperez/appli/ClaDS_Julia/")

include("/users/biodiv/bperez/appli/ClaDS_Julia/load_ClaDS2_functions.jl")

cd("/users/biodiv/bperez/data/others/Benji/diversification/")

R"library(RPANDA)"


name="primate_complete"

@rput name


# assume 504 species based on https://advances.sciencemag.org/content/3/1/e1600946/tab-pdf
# sampling fraction 367/504


R"tree <- read.tree(paste0('tree_',name,'.tre'))"


@rget tree     #load it in Julia
tree = ape2Tree(tree)    # and make it a Julia Tree

Random.seed!(813)
n_iter = 400

sample_fraction = 0.67
sampler = run_ClaDS2(tree, n_iter,plot_tree = 0,
    print_state = Int64(n_iter/4), plot_chain = false,
    f=sample_fraction, ltt_steps = 250)

sampler_to_Rdata(tree, sampler, join(["ClaDS2_tree_", name, "_f67.Rdata"]) ; sample_fraction = sample_fraction, max_it_number = 2_500)
result = (tree, sampler, sample_fraction)
@save join(["ClaDS2_tree_", name, "_f67.jld2"]) result


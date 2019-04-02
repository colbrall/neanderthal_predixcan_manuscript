# human_dr.jl
# @author Laura Colbran
# 
# given a set of IDs, will pull them out of a population and call DR status for all predixcan models
#
# julia 1.1

using DataFrames
using CSV
using ArgParse
using GZip
using Statistics
using HypothesisTests
using GLM
using StatsPlots
pyplot(grid = false)

# parses command-line arguments
function parseCommandLine()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--predictions","-p"
            help = "path to directory containing population predixcan predictions"
            arg_type = String
            required = true
        "--ids", "-i"
            help = "path to file containing IDs for people we're calling DR status in"
            arg_type = String
            required = true
        "--comp_file","-c"
            help = "file of DR genes to compare to set we're calculating if you want"
            arg_type = String
        "--write","-w"
            help = "if you want to save list of DR genes for IDs"
            action = :store_true
        "--superpop"
            help = "if you want to compare across superpopulations"
            action = :store_true
    end
    return parse_args(s)
end

# returns tissue name from prediction file
function getTissue(file::String)::String
    return join(split(file,"")[1:end-24])
end

# my tissue names -> Eric's tissue names
function mapNames()
  return Dict{String,String}(
    "adipose_subcutaneous" => "Adipose-Subcutaneous",
    "adipose_visceral_omentum" => "adipose_visceral_omentum",
    "brain_putamen_basal_ganglia" => "Brain-Putamen-basalganglia",
    "pancreas" => "Pancreas", "breast_mammary_tissue" => "Breast-MammaryTissue",
    "pituitary" => "Pituitary", "adrenal_gland" => "AdrenalGland",
    "cells_ebv_transformed_lymphocytes" => "cells_ebv-transformed_lymphocytes",
    "anterior_cingulate_cortex" => "Brain-Anteriorcingulatecortex-BA24",
    "cells_transformed_fibroblasts" => "Cells-Transformedfibroblasts",
    "skin_nosun_suprapubic" => "Skin-NotSunExposed-Suprapubic",
    "artery_aorta" => "Artery-Aorta", "colon_sigmoid" => "Colon-Sigmoid",
    "skin_sun_lower_leg" => "Skin-SunExposed-Lowerleg",
    "artery_coronary" => "Artery-Coronary", "colon_transverse" => "Colon-Transverse",
    "small_intestine_terminal_ileum" => "SmallIntestine-TerminalIleum",
    "artery_tibial" => "Artery-Tibial", "spleen" => "Spleen",
    "esophagus_gastroesophageal_junction" => "Esophagus-GastroesophagealJunction",
    "brain_caudate_basal_ganglia" => "Brain-Caudate-basalganglia",
    "esophagus_mucosa" => "Esophagus-Mucosa", "stomach" => "Stomach",
    "brain_cerebellar_hemisphere" => "Brain-CerebellarHemisphere",
    "esophagus_muscularis" => "Esophagus-Muscularis",
    "brain_cerebellum" => "Brain-Cerebellum", "nerve_tibial" => "Nerve-Tibial",
    "heart_atrial_appendage" => "Heart-AtrialAppendage",
    "brain_cortex" => "Brain-Cortex", "liver" => "Liver", "lung" => "Lung",
    "brain_frontal_cortex" => "Brain-FrontalCortex-BA9",
    "brain_hippocampus" => "Brain-Hippocampus", "muscle_skeletal" => "Muscle-Skeletal",
    "whole_blood" => "WholeBlood", "brain_hypothalamus" => "Brain-Hypothalamus",
    "brain_nucleus_accumbens_basal_ganglia" => "Brain-Nucleusaccumbens-basalganglia",
    "prostate" => "prostate","ovary" => "Ovary", "vagina" => "vagina", "testis" => "testis",
    "uterus" => "uterus", "thyroid" => "thyroid","left_ventricle" => "left_ventricle")
end

# reads population expression values into a dictionary, and splice out individual of interest
function readDict(id::String,pop_path::String)
    pop_dict = Dict{String,Array{Float64,1}}()
    id_dict = Dict{String,Float64}()
    GZip.open(pop_path) do f
        e = Array{String,1}()
        ind = 0
        for line in eachline(f)
            if startswith(line,"gene")
                ind = findfirst(x-> x==id,split(chomp(line),"\t")[2:end])
                continue
            end
            if startswith(line,"ENSG")
                e = [parse(Float64,d) for d=split(chomp(line),"\t")[2:end]]
                id_dict[split(chomp(line),"\t")[1]] = splice!(e,ind)
                pop_dict[split(chomp(line),"\t")[1]] = e
            end
        end
    end
    return id_dict,pop_dict
end

# calculates empirical 2-sided p-value based on population
function popP(id::String,pop_path::String)
    #read into dictionary {gene -> pop}
    person,pop = readDict(id,pop_path)
    genes = []
    for gene in keys(pop) #iterate through genes
        stat = person[gene]
        p = 0.0::Float64
        m = median(pop[gene])
        for pers in pop[gene]
            if abs(m - pers) >= abs(m - stat)
                p += 1
            end
        end
        p = p/(length(pop[gene]))
        if p == 0 
            push!(genes,gene)
        end
    end
    return genes #list of dr genes
end

# returns list of DR genes for a person
function callDR(pred_dir::String,id_path::String)
    ids = CSV.read(id_path; delim='\t',header=1)
    dr_genes = DataFrames.DataFrame(tissue = String[],person = String[],gene=String[])
    empty_ids = DataFrames.DataFrame(tissue=String[],id=String[])
    for item in readdir(pred_dir)
        if !endswith(item,".full.gz") continue end
        for person in ids[:,1]
            empty = true
            for gene in popP(person,"$(realpath(pred_dir))/$(item)")
                empty=false
                push!(dr_genes, [getTissue(item),person,gene])
            end
            if empty
                push!(empty_ids, [getTissue(item),person])
            end
        end
    end
    return dr_genes, empty_ids
end

# compares DR genes between random set and another individual
function compareDR(dr_genes::DataFrames.DataFrame,empty_ids::DataFrames.DataFrame,comp_file::String)
    comp = CSV.read(comp_file; delim='\t',header=[:tissue,:p,:gene])
    deletecols!(comp, [:p])
    for tiss in 1:nrow(comp)
        comp[tiss,:tissue] = mapNames()[comp[tiss,:tissue]]
    end
    println("Overall: $(length(unique(comp[:gene])))")
    num_tiss = length(unique(dr_genes[:tissue]))
    people = unique(dr_genes[:,[:person,:gene]])
    people = by(people,:person,dr_count = :gene => y->count(!ismissing,y))
    summ_ids = by(empty_ids,:id, missing_count = :tissue => y -> count(!ismissing,y))
    for id in summ_ids[summ_ids[:missing_count].==num_tiss,:id]
        push!(people,[id,0])
    end
    describe(people[:dr_count])
    println(OneSampleTTest(mean(people[:dr_count]),std(people[:dr_count]),nrow(people),length(unique(comp[:gene]))))
    println("Tissue-Wise:")
    comp_tissue_wise = by(comp,:tissue,dr_count = :gene => y->count(!ismissing,y))
    for tiss in comp_tissue_wise[:tissue]
        if !(tiss in unique(dr_genes[:tissue])) continue end
        println("$tiss: $(comp_tissue_wise[comp_tissue_wise[:tissue].== tiss,:dr_count][1])")
        people = by(dr_genes[dr_genes[:tissue].== tiss,:],:person,dr_count = :gene => y->count(!ismissing,y))
        for per in empty_ids[empty_ids[:tissue].== tiss,:id]
            push!(people,[per,0])
        end
        describe(people[:dr_count])
        println(OneSampleTTest(mean(people[:dr_count]),std(people[:dr_count]),nrow(people),comp_tissue_wise[comp_tissue_wise[:tissue].== tiss,:dr_count][1]))
    end
end

# compares across superpops
function superPop(dr_genes::DataFrames.DataFrame,empty_ids::DataFrames.DataFrame,id_path::String)
    println("Overall:")
    num_tiss = length(unique(dr_genes[:tissue]))
    people = unique(dr_genes[:,[:person,:gene]])
    people = by(people,:person,dr_count = :gene => y->count(!ismissing,y))
    summ_ids = by(empty_ids,:id, missing_count = :tissue => y -> count(!ismissing,y))
    for id in summ_ids[summ_ids[:missing_count].==num_tiss,:id]
        push!(people,[id,0])
    end
    ids = CSV.read(id_path,delim='\t')[:,[:ID,:SuperPop]]
    rename!(ids,:ID => :person)
    people = join(people, ids, on=:person,kind=:left)
    for pop in unique(people[:SuperPop])
        println("\n$pop:")
        describe(people[people[:SuperPop].== pop,:dr_count])
    end
    model = lm(@formula(dr_count ~ 1 + SuperPop),people)
    nullmodel = lm(@formula(dr_count ~ 1),people)
    println("\nF Test for effect of Superpopulations on DR Count:")
    println(ftest(model.model,nullmodel.model))
    @df people density(:dr_count,
                    xlabel = "DR GeneCount",
                    ylabel = "Density",
                    group = :SuperPop, 
                    legend = :topleft)
    savefig("superpop_distribution.pdf")
end

function main()
    parsed_args = parseCommandLine()
    dr_genes, empty_ids = callDR(parsed_args["predictions"],parsed_args["ids"])
    if parsed_args["write"]
        CSV.write("dr_genes.txt",dr_genes,delim='\t')
    end
    if typeof(parsed_args["comp_file"]) != Nothing
        compareDR(dr_genes,empty_ids,parsed_args["comp_file"])
    end
    if parsed_args["superpop"]
        superPop(dr_genes,empty_ids,parsed_args["ids"])
    end
end

main()

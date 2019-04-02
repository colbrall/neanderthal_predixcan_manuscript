# comp_pop.jl
# @author Laura Colbran
# functions to compare PrediXcan results in archaics to Eric's bioVU results or 1kG
#
# contains functions for calculating empirical p-values, calling DR genes, direction bias, PCA,
# and computing a distance matrix for a tree
# 
# runs with Julia version 1.1 and requires R3.3 or greater to be present in the environment with package dendextend

using ArgParse
using GZip
using DataFrames
using Gadfly
using Cairo
using StatsBase
using MultivariateStats
using Distances
#using Clustering
using HypothesisTests
using CSV
using RCall

const DENDRO_SOURCE = "/dors/capra_lab/projects/neanderthal_predixcan/bin/WriteDendrogram.R"
#const N = 18621 #number people in population for empirical p-value calculation
const N = 2504

################################################################################
# USEFUL FUNCTIONS-- called by other functions
################################################################################

# parses command-line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--expr_mat","-e"
            help = "path to tissue x gene predicted expression file[s]"
            nargs = '*'
            arg_type = String
        "--pop_dir","-d"
            help = "path to directory containing population predictions (1 file per tissue)"
            arg_type = String
        "--emp_p"
            help = "if we're calculating an empirical p-value. requires --expr_mat and --pop_dir"
            action = :store_true
        "--p_mat","-p"
            help = "path to tissue x gene p-value file[s]"
            nargs = '*'
            arg_type = String
        "--manhattan" 
            help = "if you want to calculate statistics based on the p-values and make a manhattan plot. requires --p_mat"
            action = :store_true
        "--avg_p"
            help = "to summarize statistics and p-values across tissues. requires --p_mat"
            action = :store_true
        "--fdr"
            help = "returns how many genes pass FDR significance in each tissue. requires --p_mat"
            action = :store_true     
        "--pop_ids","-i"
            help = "path to file with IDs for all individuals being considered"
            arg_type = String
        "--pca"
            help = "to run a PCA with population and given number of individuals. requires --expr_mat --pop_dir and --pop_ids"
            action = :store_true    
        "--filter","-f"
            help = "path to gene list to filter on for tree-building. optional."
            arg_type = String
        "--exclude"
            help = "whether to exclude genes (as opposed to include only those). Default to false."
            action = :store_true 
        "--spearman"
            help = "whether to use spearman correlation for trees (rather than pearson). Default to false."
            action = :store_true                 
        "--tree"
            help = "hierarchical clustering to return newick tree. requires --expr_mat --pop_dir --pop_ids and WriteDendrogram.R"
            action = :store_true       
        "--bias"
            help = "calculate direction bias. requires --expr_mat --pop_dir --filter"
            action = :store_true  
    end
    return parse_args(s)
end

# tracks Julia memory use. N.B.: will crash if there's multiple instances of Julia running
function memUse()
  pid = parse(Int,chomp(readstring(pipeline(`ps axc`,`awk "{if (\$5==\"julia\") print \$1}"`))))
  return string(round(Int,parse(Int,chomp(readstring(`ps -p $pid -o rss=`)))/1024),"M")
end

# returns bonferroni multiple testing correction significance threshold
function bonferroni(m::Int64)
  return 0.05/m
end

# returns Benjamini Hochberg FDR multiple testing correction significance threshold
function FDR(a::Array{Float64,1})
  fdr = 0.0
  f = DataFrame(values = sort!(a), rank = collect(1:length(a)))
  for i in 1:nrow(f)
    if (0.05*f[i,2]/nrow(f)) < f[i,1]
      fdr = f[i,1]
      break
    end
  end
  return fdr
end

# converts p-value into one with pseudocount of 1
function pseudocount(d::Array{Float64,1},N::Int64)
  return ((d*N)+1)/(N+1)
end

################################################################################
# INPUT-- all called by other functions
################################################################################

# reads gene x tiss matrix into DataFrame
function readGeneDF(f_path::String)
  println("Reading file $f_path as data frame.......")
  if ispath(f_path)
    df = readtable(f_path,separator='\t')
    try rename!(df, :gene=> :gene_id) finally return df end
  else
    df = "NA"
    return df
  end
end

# reads file into DataFrame
function readDF(f_path::String)
  println("Reading file $f_path as data frame.......")
  if ispath(f_path)
    df = readtable(f_path,separator='\t', normalizenames=false)
  else
    df = "NA"
  end
  return df
end

# reads population expression values into a dictionary
function readDict(path::String)
  dict = Dict{String,Array{SubString{String},1}}()
  try
    gzopen(path) do f
      e = Array(String,1)
      for line in eachline(f)
        if startswith(line,"ENSG")
          e = split(chomp(line),"\t")
          dict[ascii(e[1])] = e[2:end]
        end
      end
    end
    return dict
  catch
    return "NA"
  end
end

# reads population and ids into a dictionary
function popDict(path::String)
  dict = Dict{SubString{String},SubString{String}}()
  try
    gzopen(path) do f
      for line in eachline(f)
        e = split(chomp(line),"\t")
        dict[e[1]] = e[2]
      end
    end
    return dict
  catch
    return "NA"
  end
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

# reverses key:value dictionary; naive-- doesn't check for unique values
function flipDict(dict::Dict{String,String})
  new_dict = Dict{String,String}()
  for key in keys(dict)
    new_dict[dict[key]] = key
  end
  return new_dict
end

#1kG pop -> superpop
function thousGenSuperPopDict()
  return Dict{String,String}(
    "GWD" => "AFR", "MSL" => "AFR","ESN" => "AFR", "MXL" => "AMR", "CLM" => "AMR",
    "PEL" => "AMR", "TSI" => "EUR", "IBS" => "EUR", "PJL" => "SAS", "STU" => "SAS",
    "ITU" => "SAS", "GBR" => "EUR", "CHB" => "EAS", "JPT" => "EAS", "CDX" => "EAS",
    "YRI" => "AFR", "LWK" => "AFR", "CEU" => "EUR", "GIH" => "SAS", "ASW" => "AFR",
    "CHS" => "EAS", "KHV" => "EAS", "ACB" => "AFR", "PUR" => "AMR", "BEB" => "SAS",
    "FIN" => "EUR")
end


################################################################################
# OPERATIONS
################################################################################

# calculates empirical 2-sided p-value based on population, saves new file: genes x tissues
function popP(mat::DataFrames.DataFrame,pop_path::String,tiss::Dict{String,String},inf::String,side::Int64)
  #const inds = collect(1:length(mat[1]))
  for j in keys(tiss) #iterate through tissues
    #read into dictionary {gene -> pop}
    println("Calculating p-value for $(j)")
    pop = readDict("$(realpath(pop_path))/$(tiss[j])_elasticNet0_0.5.full.gz")
    if pop == "NA"
      println("$(j) not in Population Directory")
      try deletecols!(mat,Symbol(j)) finally continue end #remove tissue from mat if it's there
    end
    for i in 1:nrow(mat) #iterate through genes
      stat = mat[Symbol(j)][i]
      try
        p = 0.0::Float64
        m = median([parse(Float64,d) for d=pop[ascii(mat[:gene_id][i])]])
        for person in pop[ascii(mat[:gene_id][i])]
          #println("med: $m, stat: $stat, person: $person")
          if side == 2 && abs(m - parse(Float64,person)) >= abs(m - stat)
            p += 1
          elseif side == 1 && ((stat > m && parse(Float64,person) >= stat) || (stat < m && parse(Float64,person) <= stat))
            p += 1
          end
        end
        mat[Symbol(j)][i] = (p/(length(pop[ascii(mat[:gene_id][i])]))) #replace stat with p
    catch
        mat[Symbol(j)][i] = NA #if gene id doesn't exist in population file, or doesn't have expression in that tissue
      end
    end
  end
  # write new file
  if side == 2
    out_f = join(vcat([splitext(inf)[1]],["_p_values_2sided.txt"]),"")
  else
    out_f = join(vcat([splitext(inf)[1]],["_p_values_1sided.txt"]),"")
  end
  CSV.write(out_f, mat; delim='\t')
  return mat
end

# for each gene, averages the p-values in all tissues it was modeled
function avgP(df::DataFrames.DataFrame)
  df[:mean] = 0.0
  for row in 1:nrow(df)
    num_val = 0.0
    has_val = false
    for col in 2:(ncol(df)-1)
      if isna(df[row,col])
        continue
      else
        has_val = true
        num_val += 1
        df[:mean][row] += df[row,col]
      end
    end
    if has_val
      df[:mean][row] = df[:mean][row]/num_val
    else
      df[:mean][row] = NA
    end
  end
  return DataFrame(gene_id = df[:gene_id], mean_p = df[:mean])
end

# does PCA for population plus individuals
function runPCA(expr_mat::Array{String,1},id_path::String,pop_dir::String)
  println("Running PCA.....")
  map_dict = mapNames()::Dict{String,String}
  for key in keys(map_dict)
    pop = readGeneDF("$(realpath(pop_dir))/$(map_dict[key])_elasticNet0_0.5.full.gz")
    if pop == "NA"
      println("$key not in pop directory")
      continue
    end
    df = DataFrame()
    for file in expr_mat #get values for individuals with gen_mat file
      tmp = readtable(file, separator='\t')
      if !(Symbol(key) in names(tmp))
        println("$key not in target file(s)")
      end
      if nrow(df) == 0
        df[:gene_id] = tmp[:gene_id]
        df[Symbol(splitext(basename(file))[1])] = tmp[Symbol(key)]
        continue
      end
      t = DataFrame()
      t[:gene_id] = tmp[:gene_id]
      t[Symbol(splitext(basename(file))[1])] = tmp[Symbol(key)]
      df = join(df, t, on = :gene_id) #default inner-- only want genes that exist in both
    end
    df = join(df, pop, on = :gene_id) # gene x inds
    dropmissing!(df)
    deletecols!(df, :gene_id)
    df_array = Array(df)
#    for row in 1:size(df_array)[1]
#      df_array[row,:] = zscore(df_array[row,:]) #doesn't actually do much- already ~zscores
      #df_array[row,:] = zscore(abs(df_array[row,:]))
      #df_array[row,:] = zscore(df_array[row,:]+ -minimum(df_array[row,:]) + 1)
      #df_array[row,:] = zscore(log2(df_array[row,:]+ -minimum(df_array[row,:]) + 1))
#    end
    pca = fit(PCA,df_array)
    pcaPlot(DataFrame(transpose(MultivariateStats.transform(pca,df_array))),names(df),id_path,key)
  end
end

# hierarchical clustering for population plus individuals
function hierClust(expr_mat::Array{String,1},id_path::String,pop_dir::String,filter,to_excl::Bool,spear::Bool)
  println("Clustering.....")
  excl = ["MXL", "CLM", "PUR", "ACB","ASW", "PJL", "PEL"]
  map_dict = mapNames()::Dict{String,String}
  for key in keys(map_dict)
    num_ind = length(expr_mat)
    if !isdirpath(pop_dir)
      pop = "none"
      println("No Population Directory")
    else
      pop = readGeneDF("$(realpath(pop_dir))/$(map_dict[key])_elasticNet0_0.5.full.gz")
    end
    if pop == "NA"
      println("$key not in pop directory")
      continue
    end
    df = DataFrame()
    for file in expr_mat #get values for individuals with gen_mat file
      tmp = readtable(file, separator='\t')::DataFrames.DataFrame
      if !(Symbol(key) in names(tmp))
        println("$key not in target file(s)")
        continue
      end
      if nrow(df) == 0
        df[:gene_id] = tmp[:gene_id]
        df[Symbol(splitext(basename(file))[1])] = tmp[Symbol(key)]
        continue
      end
      t = DataFrame()
      t[:gene_id] = tmp[:gene_id]
      t[Symbol(splitext(basename(file))[1])] = tmp[Symbol(key)]
      df = join(df, t, on = :gene_id) #default inner-- only want genes that exist in both
    end
    dropmissing!(df)
        if pop != "NA" df = join(df, pop, on = :gene_id) end # gene x inds, if there's a pop file
    pop_dict = popDict(id_path)
    # remove people to exclude (admixed populations)
    for ind in keys(pop_dict)
      if in(pop_dict[ind], excl)
          try deletecols!(df,Symbol(ind)) finally continue end
      end
    end
    # filter gene list
    if typeof(filter) != Nothing
        filter_genes = CSV.read(filter, delim='\t',header=["gene_id"])
        if to_excl #ie you want to exclude the genes in your filtered list
            df = join(df,filter_genes,on=:gene_id,kind=:anti)
        else
            df = join(df,filter_genes,on=:gene_id,kind=:inner)
        end
    end
    deletecols!(df, :gene_id)    
    if spear
        for col in names(df)
            df[col] = tiedrank(df[col]) #convert to ranks so CorrDist is actually Spearman correlation.
        end
    end
    #writedlm("$(key)_distance.txt",pairwise(CorrDist(),Matrix(df)))
    makeTree(pairwise(CorrDist(),Matrix(df),dims=2),names(df),pop_dict,key)
  end
end

# saves newick tree, has commented-out code to plot the tree as well.
# called by hierClust()
function makeTree(dist::Array{Float64,2},ids::Array{Symbol,1},pop_dict::Dict{SubString{String},SubString{String}},tiss::String)
    out_tree = "$(tiss)_dendrogram.newick"
    R"""
        suppressPackageStartupMessages(library(dendextend))
        source($DENDRO_SOURCE)
        dend <- as.dendrogram(hclust(as.dist($dist)), hang = -1)
        l <- order.dendrogram(dend)
        WriteDendrogram(dend,file=$out_tree,quoteLabels = F)
    """
## the rest plots a very ugly tree. in practice, I open the newick tree in FigTree to make it pretty.
  # out_f = "$(tiss)_dendrogram.pdf"
  #id = []
#  grp = []
#  for i in 1:length(ids)
#    #id = [id,string(ids[i])]
#    grp = [grp,pop_dict[string(ids[i])]]
#  end
#  @rget l
#  poplabels = Array{String,1}()
#  superpop = thousGenSuperPopDict()
#  for i in 1:length(l)
#    try
#      #poplabels = [poplabels,superpop[grp[l[i]]]] #if dividing into superpops
#      poplabels = [poplabels,grp[l[i]]] #if leaving as subpops
#    catch
#      poplabels = [poplabels,grp[l[i]]]
#    end
#  end
#  R"""
#  suppressPackageStartupMessages(library(ggplot2))
#  mycols <- c("darkgrey", "tan3", "red4", "deeppink", "dodgerblue", "yellow",
#    "navy", "turquoise", "chartreuse3", "grey28", "orangered", "gold3",
#    "darkgoldenrod", "darkorchid", "burlywood4", "orange", "hotpink", "steelblue",
#    "blue", "red", "tomato3", "darkolivegreen", "black", "palegreen3", "violetred",
#    "forestgreen", "purple4", "darkorange1", "sienna4") # for all subpops
#  #mycols <- c("grey69", "red4","forestgreen", "grey28", "navy", "yellow1",
#    #"black", "darkorchid") # All superpops
#  #mycols <- c("grey69", "red4","forestgreen", "navy", "yellow1",
#    #"black", "darkorchid") # if excluding DEN
#  palette(mycols)
#  """
#  R"""
#  dend <- set(dend, "labels", $(poplabels))
#  l <- sort(unique(labels(dend)))
#  dend <-set(dend,"branches_lwd",0.3)
#  dend <-sort(dend,type = "nodes")
#  dend <-set(dend,"labels_cex",0.3)
#  col <- 2
#  for (pop in l) {
#    dend <-set(dend,"by_labels_branches_col", value = c(pop), TF_values = c(col,Inf))
#    col = col+1
#  }
#  """
#  R"""
#  ggd <- as.ggdend(dend)
#  #print(ggd$labels[4])
#  p <- ggplot(ggd,labels=TRUE) + ylim(0, max(get_branches_heights(dend))) #+
#  #geom_point(data=ggd$labels,aes(x=x,y=y,color=col)) +
#  #geom_text(data=ggd$labels,aes(x=x,y=y,label=label))
#  ggsave($out_f,plot=p,width = 15, height = 7, units = c("in"))
#  """
end

# calculates direction bias of
function dirBias(gene_list::DataFrames.DataFrame,mat::DataFrames.DataFrame,pop_path::String,tiss::Dict{String,String})
println("Calculating direction bias")
  mat = join(mat, gene_list,on=:gene_id,kind=:inner)
  props = Dict{String,Array{Int64,1}}()
  total_up = 0
  total = 0
  for t in keys(tiss) #iterate through tissues
    #read into dictionary {gene -> pop}
    pop = readDict("$(realpath(pop_path))/$(tiss[t])_elasticNet0_0.5.full.gz")
    if pop == "NA"
      println("$(t) not in Population Directory")
      try deletecols!(mat,Symbol(t)) finally continue end #remove tissue from mat if it's there
    end
    num_genes = 0
    num_up = 0
    for gene in 1:nrow(mat) #iterate through genes
      sample = mat[gene,Symbol(t)]
      if isna(sample) continue end #skip if not in this tissue
      exp_dist = [parse(Float64,d) for d=pop[ascii(mat[:gene_id][gene])]]
      if (sample >= minimum(exp_dist) && sample <= maximum(exp_dist)) continue end #skip if not DR in this tissue
      m = median(exp_dist)
      total += 1
      num_genes += 1
      if sample > m
        total_up += 1
        num_up += 1
      end
    end
    props[tiss[t]] = [num_up,num_genes]
  end
  exp_prop = total_up/total #calculate proportion across all models
  println("Overall Prop. Up = $exp_prop")
  println(BinomialTest(total_up,total,0.5))
  for t in keys(props)
    println("$t:")
    println(BinomialTest(props[t][1],props[t][2],exp_prop))
  end
end

################################################################################
# OUTPUT
################################################################################

# plots p-value distribution for all tissues
function manhattanPlot(df::DataFrames.DataFrame,out_f::String)
  println("Making manhattan plot......")
  gap = 10 #for gap between tissues in plot
  s = melt(df,[1])::DataFrames.DataFrame
  dropmissing!(s) #remove NA
  s[:log] = -log10(pseudocount(s[:value],N))
  s[:x] = collect(1:nrow(s)) #set spacing for genes, tissues
  add = 0
  prev = s[1,1]
  for i in 1:nrow(s)
    if (s[i,1] != prev)
      add += gap
      prev = s[i,1]
    end
    s[i,5] += add
  end
  bonf = -log10(bonferroni(nrow(s)))
   fdr = -log10(FDR(deepcopy(s[:value])))
  f = DataFrame(x = s[:x], fdr = fdr, bonf = bonf)
  theme = Theme(highlight_width=0pt,key_position = :bottom,key_max_columns=1,background_color=colorant"white",
        major_label_font_size=20pt,minor_label_font_size=16pt,key_label_font_size=12pt,key_title_font_size=14pt)
  draw(PDF(out_f,40cm,30cm),plot(layer(f,x="x",y="bonf",Geom.line,Theme(default_color=colorant"red")),
      layer(f,x="x",y="fdr",Geom.line,Theme(default_color=colorant"purple")),layer(s, color="variable", x="x", y="log", Geom.point),
      Coord.Cartesian(xmin=0,xmax=maximum(s[:x])),theme,Guide.xlabel("Tissues"),Guide.xticks(label=false),Guide.ylabel("-log10(P)")))
end

#save PC1-PC2, PC2-PC3, PC1-PC3 plots
# called by runPCA()
function pcaPlot(pca::DataFrames.DataFrame,ids::Array{Symbol,1},pop_f::String,tiss::String)
  println("making PCA plots....")
  #id = []
  grp = []
  superpop = thousGenSuperPopDict()
  for i in 1:length(ids)
    try
      grp = [grp,superpop[pop_dict[string(ids[i])]]] #if dividing into superpops
    catch
      grp = [grp,pop_dict[string(ids[i])]]
    end
  end
  #for i in 1:length(ids)
  #  id = [id,string(ids[i])]
  #  grp = [grp,pop_dict[string(ids[i])]]
  #end
  #pca[:ids] = id
  pca[:Population] = grp
  theme = Theme(highlight_width=0pt,background_color=colorant"white",
    major_label_font_size=20pt,minor_label_font_size=16pt,key_label_font_size=12pt,key_title_font_size=14pt)
  if ncol(pca) > 5
    draw(PDF("$(tiss)_pc1pc2.pdf",40cm,30cm),plot(layer(pca,color ="Population",x="x1",y="x2",Geom.point,
      Theme(default_point_size=4pt,default_color=colorant"blue")),
      theme,Guide.xlabel("PC1"),Guide.ylabel("PC2")))
    draw(PDF("$(tiss)_pc1pc3.pdf",40cm,30cm),plot(layer(pca,color ="Population",x="x1",y="x3",Geom.point,
      Theme(default_point_size=4pt,default_color=colorant"blue")),
      theme,Guide.xlabel("PC1"),Guide.ylabel("PC3")))
    draw(PDF("$(tiss)_pc2pc3.pdf",40cm,30cm),plot(layer(pca,color ="Population",x="x2",y="x3",Geom.point,
      Theme(default_point_size=4pt,default_color=colorant"blue")),
      theme,Guide.xlabel("PC2"),Guide.ylabel("PC3")))
  elseif ncol(pca) > 4
    draw(PDF("$(tiss)_pc1pc2.pdf",40cm,30cm),plot(layer(pca,color ="Population",x="x1",y="x2",Geom.point,
      Theme(default_point_size=4pt,default_color=colorant"blue")),
      theme,Guide.xlabel("PC1"),Guide.ylabel("PC2")))
  else
    println("Not enough dimensions...")
  end
end

# writes file with summary info about a matrix of gene expression
# ARGS: PATH/TO/PVALUE/FILE
function vitalStats(df::DataFrames.DataFrame,out_path::String)
  println("Calculating stats....")
  out_f = join(vcat(out_path,["_stats.txt"]),"")
  f = open(out_f,"w")
  write(f,"Summary Stats for $(splitdir(out_path)[2]):\n")
  s = melt(df,[1])::DataFrames.DataFrame
  dropmissing!(s)
  fdr = FDR(deepcopy(s[:value]))
  write(f,"$(nrow(s)) Total tests\n")
  write(f,"$(nrow(s[s[:value] .< fdr, :])) passed FDR overall\n\n")
  groups = unique(s[:variable])
  for group in groups
    write(f, "$group\n")
    write(f,"$(nrow(s[s[:variable] .== group, :])) tests\n")
    write(f,"$(nrow(s[(s[:value] .< fdr) & (s[:variable] .== group), :])) passed FDR\n\n")
  end
  close(f)
end

# writes file with all genes <= p
# ARGS: PATH/TO/PVALUE/FILE
function topGenes(df::DataFrames.DataFrame,out_path::String,p::Float64)
  println("Assembling list of top genes.....")
  out_f = join(vcat(out_path,["_top_genes.txt"]),"")
  s = melt(df,[1])::DataFrames.DataFrame
  dropmissing!(s)
  open(out_f, "w") do f
    writedlm(f,"tissue\tpvalue\tgene_id")
    writedlm(f, convert(Array,s[s[:value] .<= p,:]), '\t')
  end
end

# calculates FDR threshold for each tissue, and prints how many passed it
# ARGS: PATH/TO/PVALUE/FILE
function fdrPerTissue(df::DataFrames.DataFrame)
  try deletecols!(df, :gene_id)
  catch; deletecols!(df, :gene) end
  for tiss in names(df)
    values = dropna(df[tiss])
    thresh = FDR(values)
    println("$tiss:\nThreshold: $thresh")
    println("$(length(values[values .< thresh, :]))/$(length(values)) passed")
  end
end

################################################################################

function main()
    parsed_args = parse_commandline()
#### calculate empirical p-values
    if parsed_args["emp_p"]
        map_dict = mapNames()::Dict{String,String}
        for in_file in parsed_args["expr_mat"]
            gen_mat = readGeneDF(in_file)::DataFrames.DataFrame
            popP(gen_mat,parsed_args["pop_dir"],map_dict,in_file,2)::DataFrames.DataFrame
            #popP(gen_mat,map_dict,in_file,1)::DataFrames.DataFrame
        end
    end
    
#### manhattan plots, multiple testing
    if parsed_args["manhattan"]
        for in_file in parsed_args["p_mat"]
            p_mat = readGeneDF(in_file)::DataFrames.DataFrame
            manhattanPlot(p_mat,join(vcat([splitext(in_file)[1]],["_tissues_manhattan.pdf"]),""))
            vitalStats(p_mat,join(vcat([splitext(in_file)[1]],["_tissues"]),""))
            topGenes(p_mat,join(vcat([splitext(in_file)[1]],["_tissues"]),""),0.0)
        end
    end

    if parsed_args["avg_p"]
        for in_file in parsed_args["p_mat"]
            p_mat = readGeneDF(in_file)::DataFrames.DataFrame
            avg_mat = avgP(p_mat)::DataFrames.DataFrame
            manhattanPlot(avg_mat,join(vcat([splitext(in_file)[1]],["_avg_manhattan.pdf"]),""))
            vitalStats(avg_mat,join(vcat([splitext(in_file)[1]],["_avg"]),""))
            topGenes(avg_mat,join(vcat([splitext(in_file)[1]],["_avg"]),""),0.0)
        end
    end

#### Calculating FDR thresholds per tissue
    if parsed_args["fdr"]
        for in_file in parsed_args["p_mat"]
            fdrPerTissue(readGeneDF(in_file))
        end
    end

#### PCA
    if parsed_args["pca"]
        in_files = ARGS[:,1]::Array{String,1}
        runPCA(parsed_args["expr_mat"],parsed_args["pop_ids"],parsed_args["pop_dir"])
    end
    
#### Clustering
    if parsed_args["tree"]
        in_files = ARGS[:,1]::Array{String,1}
        hierClust(parsed_args["expr_mat"],parsed_args["pop_ids"],parsed_args["pop_dir"],parsed_args["filter"],parsed_args["exclude"],parsed_args["spearman"])
    end

#### Direction Bias
    if parsed_args["bias"]
        map_dict = mapNames()::Dict{String,String}
        for in_file in parsed_args["expr_mat"]
            gen_mat = readGeneDF(in_file)::DataFrames.DataFrame
            dirBias(readDF(parsed_args["filter"]),gen_mat,parsed_args["pop_dir"],map_dict)
        end
    end
end

main()


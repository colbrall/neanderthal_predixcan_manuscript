# comp_pop.jl
# @author Laura Colbran
# functions to compare PrediXcan results in archaics to Eric's bioVU results or 1kG
#
# USAGE: julia comp_pop.jl ARGS
#
# ARGS vary by analysis
# look in main() or above individual functions for correct calls to run certain analyses
#
# runs with Julia version 0.6.1



using GZip
using DataFrames
using Gadfly
using Cairo
using StatsBase
using MultivariateStats
using Distances
using RCall

################################################################################
# USEFUL FUNCTIONS-- called by other functions
################################################################################

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
function FDR(a::DataArrays.DataArray{Float64,1})
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
function pseudocount(d::DataArrays.DataArray{Float64,1},N::Int64)
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
    try rename!(df, :gene, :gene_id) end
  else
    df = "NA"
  end
  return df
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
  const dict = Dict{String,Array{SubString{String},1}}()
  try
    gzopen(path) do f
      const e = Array(String,1)
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
  const dict = Dict{SubString{String},SubString{String}}()
  try
    gzopen(path) do f
      const e = Array(SubString{String},1)
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
    "cells_ebv_transformed_lymphocytes" => "Cells-EBV-transformedlymphocytes",
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
# ARGS: PATH/TO/EXPR/FILE PATH/TO/POP/DIR
function popP(mat::DataFrames.DataFrame,tiss::Dict{String,String},inf::String,side::Int64)
  const pop_path = ARGS[2]::String
  #const inds = collect(1:length(mat[1]))
  for j in keys(tiss) #iterate through tissues
    #read into dictionary {gene -> pop}
    println("Calculating p-value for $(j)")
    const pop = readDict("$(realpath(pop_path))/$(tiss[j])_elasticNet0_0.5.full.gz")
    if pop == "NA"
      println("$(j) not in Population Directory")
      try delete!(mat,Symbol(j)) end #remove tissue from mat if it's there
      continue
    end
    gc() #can comment out if time becomes an issue
    for i in 1:nrow(mat) #iterate through genes
      const stat = mat[Symbol(j)][i]
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
  writetable(out_f, mat, separator='\t')
  return mat
end

# for each gene, averages the p-values in all tissues it was modeled
# ARGS: PATH/TO/PVALUE/FILE
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
# ARGS: PATHS/TO/EXPR/FILES PATH/TO/POP/IDs PATH/TO/POP/DIR
function runPCA(paths::Array{String,1})
  println("Running PCA.....")
  const map_dict = mapNames()::Dict{String,String}
  for key in keys(map_dict)
    pop = readGeneDF("$(realpath(last(paths)))/$(map_dict[key])_elasticNet0_0.5.full.gz")
    if pop == "NA"
      println("$key not in pop directory")
      continue
    end
    const df = DataFrame()
    for file in paths[1:(length(paths)-2)] #get values for individuals with gen_mat file
      tmp = readtable(file, separator='\t')
      if !(Symbol(key) in names(tmp))
        println("$key not in target file(s)")
      end
      if nrow(df) == 0
        df[:gene_id] = tmp[:gene_id]
        df[Symbol(splitext(basename(file))[1])] = tmp[Symbol(key)]
        continue
      end
      const t = DataFrame()
      t[:gene_id] = tmp[:gene_id]
      t[Symbol(splitext(basename(file))[1])] = tmp[Symbol(key)]
      df = join(df, t, on = :gene_id) #default inner-- only want genes that exist in both
    end
    df = join(df, pop, on = :gene_id) # gene x inds
    complete_cases!(df)
    gc()
    delete!(df, :gene_id)
    df_array = Array(df)
#    for row in 1:size(df_array)[1]
#      df_array[row,:] = zscore(df_array[row,:]) #doesn't actually do much- already ~zscores
      #df_array[row,:] = zscore(abs(df_array[row,:]))
      #df_array[row,:] = zscore(df_array[row,:]+ -minimum(df_array[row,:]) + 1)
      #df_array[row,:] = zscore(log2(df_array[row,:]+ -minimum(df_array[row,:]) + 1))
#    end
    gc()
    pca = fit(PCA,df_array)
    pcaPlot(DataFrame(transpose(transform(pca,df_array))),names(df),paths[length(paths)-1],key)
  end
  gc()
end

# hierarchical clustering for population plus individuals
# ARGS: PATHS/TO/EXPR/FILES PATH/TO/POP/IDs PATH/TO/POP/DIR
function hierClust(paths::Array{String,1})
  println("Clustering.....")
  const excl = ["MXL", "CLM", "PUR", "ACB","ASW", "PJL", "PEL"]
  #const excl = ["MXL", "CLM", "PUR", "DEN"]
  #const excl = []
  const map_dict = mapNames()::Dict{String,String}
  for key in keys(map_dict)
    if !isdirpath(last(paths))
      pop = "none"
      println("No Population Directory")
      num_ind = length(paths) -1
    else
      pop = readGeneDF("$(realpath(last(paths)))/$(map_dict[key])_elasticNet0_0.5.full.gz")
      num_ind = length(paths) -2
    end
    if pop == "NA"
      println("$key not in pop directory")
      continue
    end
    const df = DataFrame()
    for file in paths[1:num_ind] #get values for individuals with gen_mat file
      const tmp = readtable(file, separator='\t')::DataFrames.DataFrame
      if !(Symbol(key) in names(tmp))
        println("$key not in target file(s)")
        continue
      end
      if nrow(df) == 0
        df[:gene_id] = tmp[:gene_id]
        df[Symbol(splitext(basename(file))[1])] = tmp[Symbol(key)]
        continue
      end
      const t = DataFrame()
      t[:gene_id] = tmp[:gene_id]
      t[Symbol(splitext(basename(file))[1])] = tmp[Symbol(key)]
      df = join(df, t, on = :gene_id) #default inner-- only want genes that exist in both
    end
    complete_cases!(df)
    try df = join(df, pop, on = :gene_id) end # gene x inds, if there's a pop file
    gc()
    pop_dict = popDict(paths[num_ind + 1])
    for ind in keys(pop_dict)
      if in(pop_dict[ind], excl)
        try delete!(df,Symbol(ind)) end
      end
    end
    delete!(df, :gene_id)
    makeTree(pairwise(CorrDist(),Array(df)),names(df),pop_dict,key)
  end
end
################################################################################
# OUTPUT
################################################################################

# plots p-value distribution for all tissues
# ARGS: PATH/TO/PVALUE/FILE
function manhattanPlot(df::DataFrames.DataFrame,out_f::String)
  println("Making manhattan plot......")
  #const N = 18621 #number people empirical P was calculated from
  const N = 2504
  const gap = 10 #for gap between tissues in plot
  s = melt(df,[1])::DataFrames.DataFrame
  complete_cases!(s) #remove NA
  s[:log] = -log10(pseudocount(s[:value],N))
  s[:x] = collect(1:nrow(s)) #set spacing for genes, tissues
  const add = 0
  const prev = s[1,1]
  for i in 1:nrow(s)
    if (s[i,1] != prev)
      add += gap
      prev = s[i,1]
    end
    s[i,5] += add
  end
  const bonf = -log10(bonferroni(nrow(s)))
  const fdr = -log10(FDR(deepcopy(s[:value])))
  f = DataFrame(x = s[:x], fdr = fdr, bonf = bonf)
  const theme = Theme(highlight_width=0pt,key_position = :bottom,key_max_columns=1,background_color=colorant"white",
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
  const theme = Theme(highlight_width=0pt,background_color=colorant"white",
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

# draws and saves tree
# called by hierClust()
function makeTree(dist::Array{Float64,2},ids::Array{Symbol,1},pop_dict::Dict{SubString{String},SubString{String}},tiss::String)
  println("Making dendrogram...")
  const out_f = "$(tiss)_dendrogram.pdf"
  const out_tree = "$(tiss)_dendrogram.newick"
  #id = []
  const grp = []
  for i in 1:length(ids)
    #id = [id,string(ids[i])]
    grp = [grp,pop_dict[string(ids[i])]]
  end
  R"""
  suppressPackageStartupMessages(library(dendextend))
  suppressPackageStartupMessages(library(ggplot2))
  source("/dors/capra_lab/projects/predixcan/bin/WriteDendrogram.R")
  dend <- as.dendrogram(hclust(as.dist($(dist))), hang = -1)
  l <- order.dendrogram(dend)
  """
  R"""
  WriteDendrogram(dend,file=$out_tree,quoteLabels = F)
  """
  @rget l
  const poplabels = Array{String,1}()
  superpop = thousGenSuperPopDict()
  for i in 1:length(l)
    try
      #poplabels = [poplabels,superpop[grp[l[i]]]] #if dividing into superpops
      poplabels = [poplabels,grp[l[i]]] #if leaving as subpops
    catch
      poplabels = [poplabels,grp[l[i]]]
    end
  end
  R"""
  mycols <- c("darkgrey", "tan3", "red4", "deeppink", "dodgerblue", "yellow",
    "navy", "turquoise", "chartreuse3", "grey28", "orangered", "gold3",
    "darkgoldenrod", "darkorchid", "burlywood4", "orange", "hotpink", "steelblue",
    "blue", "red", "tomato3", "darkolivegreen", "black", "palegreen3", "violetred",
    "forestgreen", "purple4", "darkorange1", "sienna4") # for all subpops
  #mycols <- c("grey69", "red4","forestgreen", "grey28", "navy", "yellow1",
    #"black", "darkorchid") # All superpops
  #mycols <- c("grey69", "red4","forestgreen", "navy", "yellow1",
    #"black", "darkorchid") # if excluding DEN
  palette(mycols)
  """
  R"""
  dend <- set(dend, "labels", $(poplabels))
  l <- sort(unique(labels(dend)))
  dend <-set(dend,"branches_lwd",0.3)
  dend <-sort(dend,type = "nodes")
  dend <-set(dend,"labels_cex",0.3)
  col <- 2
  for (pop in l) {
    dend <-set(dend,"by_labels_branches_col", value = c(pop), TF_values = c(col,Inf))
    col = col+1
  }
  """
  R"""
  ggd <- as.ggdend(dend)
  #print(ggd$labels[4])
  p <- ggplot(ggd,labels=TRUE) + ylim(0, max(get_branches_heights(dend))) #+
  #geom_point(data=ggd$labels,aes(x=x,y=y,color=col)) +
  #geom_text(data=ggd$labels,aes(x=x,y=y,label=label))
  ggsave($out_f,plot=p,width = 15, height = 7, units = c("in"))
  """
end

# writes file with summary info about a matrix of gene expression
# ARGS: PATH/TO/PVALUE/FILE
function vitalStats(df::DataFrames.DataFrame,out_path::String)
  println("Calculating stats....")
  #const N = 18621 #number people empirical P was calculated from
  const N = 2504
  out_f = join(vcat(out_path,["_stats.txt"]),"")
  f = open(out_f,"w")
  write(f,"Summary Stats for $(splitdir(out_path)[2]):\n")
  s = melt(df,[1])::DataFrames.DataFrame
  complete_cases!(s)
  const fdr = FDR(deepcopy(s[:value]))
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
  complete_cases!(s)
  open(out_f, "w") do f
    writedlm(f, Array{String,1}(names(s)), '\t')
    writedlm(f, convert(Array,s[s[:value] .<= p,:]), '\t')
  end
end

# calculates FDR threshold for each tissue, and prints how many passed it
# ARGS: PATH/TO/PVALUE/FILE
function fdrPerTissue(df::DataFrames.DataFrame)
  try delete!(df, :gene_id)
  catch delete!(df, :gene) end
  for tiss in names(df)
    const values = dropna(df[tiss])
    const thresh = FDR(DataArray(values))
    println("$tiss:\nThreshold: $thresh")
    println("$(length(values[values .< thresh, :]))/$(length(values)) passed")
  end
end

################################################################################

function main()
#### calculate empirical p-values; ARGS: PATH/TO/EXPR/FILE PATH/TO/POP/DIR
  in_file = ARGS[1]::String
  const gen_mat = readGeneDF(in_file)::DataFrames.DataFrame
  const map_dict = mapNames()::Dict{String,String}
  popP(gen_mat,map_dict,in_file,2)::DataFrames.DataFrame
#  popP(gen_mat,map_dict,in_file,1)::DataFrames.DataFrame

#### manhattan plots; ARGS: PATH/TO/PVALUE/FILE
#  in_file = ARGS[1]::String
#  p_mat = readGeneDF(in_file)::DataFrames.DataFrame
#  manhattanPlot(p_mat,join(vcat([splitext(in_file)[1]],["_tissues_manhattan.pdf"]),""))
#  vitalStats(p_mat,join(vcat([splitext(in_file)[1]],["_tissues"]),""))
#  topGenes(p_mat,join(vcat([splitext(in_file)[1]],["_tissues"]),""),0.0)
#  avg_mat = avgP(p_mat)::DataFrames.DataFrame
#  manhattanPlot(avg_mat,join(vcat([splitext(in_file)[1]],["_avg_manhattan.pdf"]),""))
#  vitalStats(avg_mat,join(vcat([splitext(in_file)[1]],["_avg"]),""))
#  topGenes(avg_mat,join(vcat([splitext(in_file)[1]],["_avg"]),""),0.0)

#### Calculating FDR thresholds per tissue; ARGS: PATH/TO/PVALUE/FILE
#  in_file = ARGS[1]::String
#  fdrPerTissue(readGeneDF(in_file))

#### PCA; ARGS: PATHS/TO/EXPR/FILES PATH/TO/POP/IDs PATH/TO/POP/DIR
#  in_files = ARGS[:,1]::Array{String,1}
#  runPCA(in_files)

#### Clustering; ARGS: PATHS/TO/EXPR/FILES PATH/TO/POP/IDs PATH/TO/POP/DIR
#  in_files = ARGS[:,1]::Array{String,1}
#  hierClust(in_files)
end

main()

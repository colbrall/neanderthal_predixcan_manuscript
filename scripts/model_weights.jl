# model_weights.jl
# @author Laura Colbran
# 2017-11-28
#
# compares models weights between two sets of SNPs of interest

using DataFrames
using GZip
using StatsBase
using SQLite
using Gadfly
using Cairo
using HypothesisTests

# tracks Julia memory use. N.B.: will crash if there's multiple instances of Julia running
# in that case, error is "LoadError: ArgumentError: extra characters after whitespace"
function memUse()
  pid = parse(Int,chomp(readstring(pipeline(`ps axc`,`awk "{if (\$5==\"julia\") print \$1}"`))))
  return string(round(Int,parse(Int,chomp(readstring(`ps -p $pid -o rss=`)))/1024),"M")
end

function standardTheme()
  return Theme(highlight_width=0pt,key_position = :bottom,key_max_columns=1,default_color=colorant"steelblue",
        major_label_font_size=20pt,minor_label_font_size=16pt,key_label_font_size=12pt,key_title_font_size=14pt,
        boxplot_spacing=96pt)
end

function smallTheme()
  return Theme(highlight_width=0pt,key_position = :bottom,key_max_columns=1,default_color=colorant"steelblue",
        major_label_font_size=12pt,minor_label_font_size=10pt,key_label_font_size=10pt,key_title_font_size=12pt)
end

function readDF(f_path::String)
  println("Reading file $f_path as DataFrame.......")
  if ispath(f_path)
    df = readtable(f_path,separator='\t', normalizenames=false, allowcomments=true,header=false)
  else
    df = "NA"
  end
  return df
end

# called by compweights()
function getWeights(snp_df::DataFrames.DataFrame,db_dir::String,tiss_stat::String)
  n_db = length(filter(x -> endswith(x,".db"),readdir(db_dir)))
  n_done = 0
  for item in readdir(db_dir)
    if endswith(item, ".db")
      println("Pulling weights from $item")
      n_done += 1
      db = SQLite.DB("$(db_dir)/$item")
      snp_df[Symbol("$(item[4:(length(item)-11)])_stat")] =  -100.0
      snp_df[Symbol("$(item[4:(length(item)-11)])_Nmodels")] =  0.0
      if n_done == n_db
        snp_df[Symbol("Nmodel_max")] = 0.0
        snp_df[Symbol("Nmodel_min")] = 0.0
        snp_df[Symbol("Nmodel_mean")] = 0.0
        snp_df[Symbol("weights_max")] = -100.0
        snp_df[Symbol("weights_min")] = -100.0
        snp_df[Symbol("weights_mean")] = -100.0
      end
      for i in 1:nrow(snp_df)
        out = SQLite.query(db,"SELECT weight FROM 'weights' WHERE rsid='$(snp_df[i,:x4])'")
        if nrow(out) != 0
          if tiss_stat == "max"
            snp_df[i,Symbol("$(item[4:(length(item)-11)])_stat")] = maximum(abs(Array(out[:weight])))
          elseif tiss_stat == "min"
            snp_df[i,Symbol("$(item[4:(length(item)-11)])_stat")] = minimum(abs(Array(out[:weight])))
          else
            snp_df[i,Symbol("$(item[4:(length(item)-11)])_stat")] = mean(abs(Array(out[:weight])))
          end
          snp_df[i,Symbol("$(item[4:(length(item)-11)])_Nmodels")] = nrow(out)
        end
        if n_done == n_db
          snp_df[i,Symbol("Nmodel_max")] = maximum(Array(snp_df[i,5:(ncol(snp_df)-6)][:,2:2:end]))
          snp_df[i,Symbol("Nmodel_min")] = minimum(Array(snp_df[i,5:(ncol(snp_df)-6)][:,2:2:end]))
          snp_df[i,Symbol("Nmodel_mean")] = mean(Array(snp_df[i,5:(ncol(snp_df)-6)][:,2:2:end]))
          #if snp_df[i,Symbol("Nmodel_min")] < 0 ## neither nean_snps or not_nean_snps ever triggered this
          #  println("$(snp_df[i,:x4]), $(snp_df[i,Symbol("Nmodel_min")])")
          #end
          if snp_df[i,Symbol("Nmodel_max")] == 0
            snp_df[i,Symbol("Nmodel_max")] = NA
          end
          try #fails if snp isn't in any models
            snp_df[i,Symbol("weights_max")] = maximum(filter(x -> x!= -100,Array(snp_df[i,5:(ncol(snp_df)-6)][:,1:2:end])))
            snp_df[i,Symbol("weights_min")] = minimum(filter(x -> x!= -100,Array(snp_df[i,5:(ncol(snp_df)-6)][:,1:2:end])))
            snp_df[i,Symbol("weights_mean")] = mean(filter(x -> x!= -100,Array(snp_df[i,5:(ncol(snp_df)-6)][:,1:2:end])))
          end
        end
      end
    end
  end
  #println(nrow(snp_df))
  completecases!(snp_df) #remove all lines where Nmodel_max = 0; ie SNPs that weren't in any models
  #println(nrow(snp_df))
end

# stats to compare weights of SNPs
# does within-tiss and overall
# ARGS: positive_set neg_set stat db_dir
# assumes BED format for both- chr start end rsid
function compWeights(args::Array{String,1})
  positives = readDF("$(realpath(args[1]))")::DataFrames.DataFrame
  negatives = readDF("$(realpath(args[2]))")::DataFrames.DataFrame
  getWeights(positives, "$(realpath(args[4]))",args[3])
  writedlm("positives.txt",convert(Array,positives),'\t')
  #writedlm("positives.txt", convert(Array,positives), '\t')
  getWeights(negatives, "$(realpath(args[4]))",args[3])
  writedlm("negatives.txt",convert(Array,negatives),'\t')
  num_tiss= length(names(positives)[5:(ncol(positives)-6)])
  tmp = DataFrames.DataFrame(Tissue = vcat(names(positives)[5:(ncol(positives)-6)],names(positives)[5:(ncol(positives)-6)]),
    set = vcat(fill("Positives",num_tiss),fill("Negatives",num_tiss)),
    stat = 0.0,min=0.0,max=0.0)
  println("\nTissue-Specific MWU Tests (with $(args[3]) of SNP weights across models):")
  for i in 1:num_tiss# p-values for tissue distributions
    pos = filter(x -> x != -100,positives[:,i+4])
    neg = filter(x -> x != -100,negatives[:,i+4])
    try
      println("$(names(positives)[i+4]) p = $(pvalue(MannWhitneyUTest(pos,neg)))")
    catch e
      println("$(names(positives)[i+4]) p = NA")
    end
    # if length(filter(x -> x < 0,pos)) > 0
    #   println("pos: $(length(filter(x -> x < 0,pos)))")
    # end
    # if length(filter(x -> x < 0,neg)) > 0
    #   println("neg: $(length(filter(x -> x < 0,neg)))")
    # end
    if length(pos) != 0
      tmp[i,:stat] = mean(pos)
      try
        tmp[i,:min] = confint(OneSampleTTest(pos))[1]
        tmp[i,:max] = confint(OneSampleTTest(pos))[2]
      catch
        tmp[i,:min] = tmp[i,:stat]
        tmp[i,:max] = tmp[i,:stat]
      end
    else
      tmp[i,:stat] = -1
      tmp[i,:min] = -1
      tmp[i,:max] = -1
    end
    if length(neg) != 0
      tmp[i+num_tiss,:stat] = mean(neg)
      try
        tmp[i+num_tiss,:min] = confint(OneSampleTTest(neg))[1]
        tmp[i+num_tiss,:max] = confint(OneSampleTTest(neg))[2]
      catch
        tmp[i+num_tiss,:min] = tmp[i+num_tiss,:stat]
        tmp[i+num_tiss,:max] = tmp[i+num_tiss,:stat]
      end
    else
      tmp[i+num_tiss,:stat] = -1
      tmp[i+num_tiss,:min] = -1
      tmp[i+num_tiss,:max] = -1
    end
  end
  writedlm("tiss_sum_stats.txt", convert(Array,tmp), '\t')
  out_f = join(vcat([splitext("$(realpath(args[1]))")[1]],["_tissues_$(args[3]).pdf"]),"")
  #p = plot(tmp[1:2:end,:], x="Tissue", y="stat", color=:set,Geom.point,
  p = plot(tmp[1:2:end,:], x="Tissue", y="stat", ymin = :min, ymax = :max,color=:set,Geom.point, Geom.errorbar,
      smallTheme(),Guide.xlabel("Tissues"),Guide.ylabel("$(args[3])"),Guide.XTicks(orientation=:vertical))
  draw(PDF(out_f,16cm,20cm),p)
  out_f = join(vcat([splitext("$(realpath(args[1]))")[1]],["_tissues_Nmodels.pdf"]),"")
  #p = plot(tmp[2:2:end,:], x="Tissue", y="stat",color=:set,Geom.point,
  p = plot(tmp[2:2:end,:], x="Tissue", y="stat", ymin = :min, ymax = :max,color=:set,Geom.point, Geom.errorbar,
      smallTheme(),Guide.xlabel("Tissues"),Guide.ylabel("Number of Models"),Guide.XTicks(orientation=:vertical),Coord.Cartesian(ymin=0))
  draw(PDF(out_f,16cm,20cm),p)

  println("\nOverall Summary Stats:")
  for i in (ncol(positives)-5):ncol(positives)
    pos = filter(x -> x != -100,positives[:,i])
    neg = filter(x -> x != -100,negatives[:,i])
    out_f = join(vcat([splitext("$(realpath(args[1]))")[1]],["_$(names(positives)[i]).txt"]),"")
    out_pdf = join(vcat([splitext("$(realpath(args[1]))")[1]],["_$(names(positives)[i]).pdf"]),"")
    println("$(names(positives)[i]) p = $(pvalue(MannWhitneyUTest(pos,neg)))")
    if length(filter(x -> x < 0,pos)) > 0
      println("pos: $(length(filter(x -> x < 0,pos)))")
    end
    if length(filter(x -> x < 0,neg)) > 0
      println("neg: $(length(filter(x -> x < 0,neg)))")
    end
    tp = DataFrames.DataFrame(set = vcat(fill("Positives",length(pos)),fill("Negatives",length(neg))),
      value = vcat(pos,neg))
    writedlm(out_f,convert(Array,tp),'\t')
    p = plot(tp, x="set", y="value", Geom.violin, standardTheme(),
          Guide.xlabel(""),Guide.ylabel("$(names(positives)[i])"))
    draw(PDF(out_pdf,10cm,8cm),p)
  end
end

function main()
  infiles = ARGS[:,1]::Array{String,1}
  compWeights(infiles) #stat can be mean, max, min- determines how snps are summarized within a tissue
end

main()

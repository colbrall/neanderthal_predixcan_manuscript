# comp_inds.jl
# @author Laura Colbran 01-18-18
#
# comparisons of predicted expression between individuals
# usage: julia comp_inds.jl PATH/TO/IND1 PATH/TO/IND2 PATH/TO/IND/3

using DataFrames
using Gadfly
using Cairo

function readDF(f_path::String)
  println("Reading file $f_path as DataFrame.......")
  if ispath(f_path)
    df = readtable(f_path,separator='\t', normalizenames=false, allowcomments=true,header=true)
  else
    df = "NA"
  end
  return df
end

# standard theme for plots
function standardTheme()
  return Theme(highlight_width=0pt,key_position = :bottom,key_max_columns=1,default_color=colorant"steelblue",
        major_label_font_size=14pt,minor_label_font_size=12pt,key_label_font_size=10pt,key_title_font_size=12pt,
        boxplot_spacing=96pt)
end

# calculates and plots the difference of expression between all pairs of individuals given in ARGS
# assumes they're all gene x tissue matrices in the same order
function diffExpHist(args::Array{String,1})
  try run(`mkdir species_diffExpHist`) end
  for i in 1:(length(args)-1)
    const ind1 = readDF(realpath(args[i]))::DataFrames.DataFrame
    const name1 = splitext(splitdir("$(realpath(args[i]))")[2])[1]
    for j in (i+1):length(args)
      const ind2 = readDF(realpath(args[j]))::DataFrames.DataFrame
      const name2 = splitext(splitdir("$(realpath(args[j]))")[2])[1]
      for k in 2:length(names(ind1))
        const tmp = completecases!(DataFrame(x1 = ind1[names(ind1)[k]],x2 = ind2[names(ind1)[k]]))::DataFrames.DataFrame
        out_f = "species_diffExpHist/$(name1)_vs_$(name2)_diffExp_$(names(ind1)[k]).pdf"
        p = plot(DataFrame(Diff = abs(tmp[:x1]-tmp[:x2])),x = "Diff",Geom.histogram,standardTheme(),Guide.xlabel("abs(diff in exp)"),Guide.ylabel("# models"))
        draw(PDF(out_f,10cm,8cm),p)
      end
    end
  end
end

# calculates the predicted difference in expression between pairs of individuals
# output file is gene x tissue matrix containing those values
function diffExpGenes(args::Array{String,1})
  for i in 1:(length(args)-1)
    const ind1 = readDF(realpath(args[i]))::DataFrames.DataFrame
    const name1 = splitext(splitdir("$(realpath(args[i]))")[2])[1]
    for j in (i+1):length(args)
      const ind2 = readDF(realpath(args[j]))::DataFrames.DataFrame
      const name2 = splitext(splitdir("$(realpath(args[j]))")[2])[1]
      const out_df = DataFrame(gene=[],x1=[], x2=[], difference=[],tissue=[])
      for k in 2:length(names(ind1))
        const tmp = completecases!(DataFrame(gene = ind1[names(ind1)[1]],x1 = ind1[names(ind1)[k]],x2 = ind2[names(ind1)[k]]))::DataFrames.DataFrame
        tmp[:difference] = abs(tmp[:x1]-tmp[:x2])
        tmp[:tissue] = "$(names(ind1)[k])"
        if nrow(out_df) == 0
          out_df = tmp
        else
          out_df = append!(out_df, tmp)
        end
      end
      out_f = "$(name1)_vs_$(name2).txt"
      names!(out_df.colindex, map(parse, ["gene",name1, name2, "difference", "tissue"]))
      writedlm(out_f, convert(Array,out_df), '\t')
    end
  end
end

function main()
  args = ARGS[:,1]::Array{String,1}
  #diffExpHist(args)
  diffExpGenes(args)
end

main()

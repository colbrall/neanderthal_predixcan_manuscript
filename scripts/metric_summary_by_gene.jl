# @author Laura Colbran
#
# runs stats on gene metric in Gene regions 
# requires R3.3 or greater and packages dplyr, dunn.test
# julia 1.1

using ArgParse
using CSV
using DataFrames
using StatsBase
using Statistics
using HypothesisTests
using StatsPlots
using RCall
pyplot(grid = false,leg = false)

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--metric","-m"
            help = "path to file with metric of interest. format gene id, metric, tab-delim"
            arg_type = String
            required = true
        "--label","-l"
            help = "mabel of metric to use for axes (e.g. Gene Density)"
            arg_type = String
            required = true
        "--regions", "-r"
            nargs='*'
            help = "list of paths to files containing gene lists for each category to compare"
            arg_type = String
        "--bins", "-b"
            help = "parameters for binning in proportion or OR analyses. format: MAX_BIN,NUM_BINS"
            arg_type = String
        "--targetlist", "-t"
            help = "list of gene IDs for your target set (usually DR genes)"
            arg_type = String
	"--only"
	    help="if you Only want the regions specifies by -r"
	    action = :store_true
        "--exclude", "-e"
            nargs='*'
            help = "file(s) for genes to exclude"
            arg_type = String
        "--summary"
            help = "if you want summary stats and a density plot of the gene metric in your metric file"
            action = :store_true
        "--comparison"
            help = "if you want compare stats in certain genes. requires --regions"
            action = :store_true
        "--proportion"
            help = "if you are comparing proportions of a gene category across gene metric. requires --comparison, --regions, --targetlist"
            action = :store_true
        "--odds"
            help = "to calculate odds ratio of target set by genes and matched by gene metric. requires --comparison, --regions, --targetlist"
            action = :store_true
    end
    return parse_args(s)
end

# pulls label for set out of file name
function labelSet(file_name::String)
    return split(basename(file_name),'.')[1]
end

# removes version info from array of Ensembl IDs
function idOnly(ids::Array{String,1})
    for i in 1:length(ids)
        ids[i] = split(ids[i],'.')[1]
    end
    return ids
end

# bins a distribution of gene metric, returns array of lower bounds of bins
function bins(dist::Array{Float64,1},max::Float64,num_bins::Int64)
    cut_points = collect(range(minimum(dist),max,step=max/num_bins))
    bins = zeros(length(dist))
    for i in 1:length(dist)
        j = 1
        while j <= length(cut_points) && dist[i] >= cut_points[j]
            j+=1
        end
        bins[i] = cut_points[j-1]
    end
    return bins
end

# returns a set of regions matched to target by metric. using bins. should end up with same distribution, not same number
function matchRates(targets::DataFrames.DataFrame, candidates::DataFrames.DataFrame)
    bin_counts = by(targets,:bins, gene_count = :target => y-> count(!ismissing,y))
    matched_set = DataFrames.DataFrame(gene = String[], score = Float64[], set = String[],target=Int64[],bins=Float64[])
    # calculate number of times I can sample before I empty a bin
    n = nrow(candidates)
    for i in 1:nrow(bin_counts)
        reps = div(nrow(candidates[candidates[:bins] .== bin_counts[i,:bins],:]),bin_counts[i,:gene_count])
        if reps < n n = reps end
    end
    if n > 0
        for i in 1:nrow(bin_counts)
            bin = bin_counts[i,:bins]
            count = bin_counts[i,:gene_count]
            bin_pick = candidates[candidates[:bins] .== bin,:]
            for match in sample(collect(1:nrow(bin_pick)), n*count; replace = false)
                push!(matched_set,bin_pick[match,:])
            end
        end
        else #triggered if the target set isn't consistently smaller than the candidate set across bins
        for i in 1:nrow(bin_counts)
            bin = bin_counts[i,:bins]
            count = bin_counts[i,:gene_count]
            bin_pick = candidates[candidates[:bins] .== bin,:]
            if count <= nrow(bin_pick)
                for match in sample(collect(1:nrow(bin_pick)), count; replace = false)
                    push!(matched_set,bin_pick[match,:])
                end
            else
                count = nrow(bin_pick)
                bin_pick = targets[targets[:bins] .== bin,:]
                for match in sample(collect(1:nrow(bin_pick)), count; replace = false)
                    push!(matched_set,bin_pick[match,:])
                end
            end
        end
    end
    return matched_set
end

# calculates and prints summary stats of a distribution, and plots the distribution
function summaryStats(metric_file::String,label::String,excl)
    metric = CSV.read(metric_file; delim='\t',allowmissing=:none,header=[:gene,:score],comment="#")    
    metric[:gene] = idOnly(metric[:gene])
    metric[:score] = metric[:score]*1. 
    for file_path in excl
        tmp = CSV.read(file_path; delim='\t',allowmissing=:none,header=[:gene])
        tmp[:gene] = idOnly(tmp[:gene])
        metric = join(metric,tmp,on=:gene,kind=:anti)
    end
    describe(metric[:score])
    println("Std. Dev: $(mean_and_std(metric[:score])[2])")
    dist_plot = histogram(metric[:score], nbins = 100,
                    xlabel = "$label",
                    ylabel = "Count",
                    color = [:black])
    savefig(dist_plot, "$(label)_distribution.pdf")
end

# runs stats between all sets of regions given
function genomeWide(metric::DataFrames.DataFrame,regions::Array{String,1},label::String)    
    @df metric boxplot(:set,:score, 
                    xlabel = "Region Set",
                    ylabel = "$label",
                    color=[:white],
                    outliers = false
                    )
    savefig("$(label)_boxplot.pdf")
    if length(unique(metric[:set])) > 2
        R"""
        suppressPackageStartupMessages(library(dplyr))
        suppressPackageStartupMessages(library(dunn.test))
        df = $metric
        group_by(df,set) %>% summarise(
            count = n(),
            mean = mean(score, na.rm = TRUE),
            sd = sd(score, na.rm = TRUE),
            median = median(score, na.rm = TRUE),
            IQR = IQR(score, na.rm = TRUE)
        )
        dunn.test(df$score,df$set,method="bh")
        """
    elseif length(unique(metric[:set])) == 2
        println(unique(metric[:set])[1])
        println(median(metric[metric[:set] .== unique(metric[:set])[1],:score]))
        println(unique(metric[:set])[2])
        println(median(metric[metric[:set] .== unique(metric[:set])[2],:score]))
        println(MannWhitneyUTest(metric[metric[:set] .== unique(metric[:set])[1],:score],metric[metric[:set] .== unique(metric[:set])[2],:score]))
    end
end

# labels windows from map according  to overlap with regions, runs stats and plots on sets of windows
function prop(metric::DataFrames.DataFrame, target_file::String,label::String,bin_params::String)    
    proplabel = rename!(CSV.read(target_file; delim='\t',allowmissing=:none,header=0),:Column1 => :gene)
    proplabel[:gene] = idOnly(proplabel[:gene])
    tmp = join(metric,proplabel,on=:gene,kind = :inner)
    tmp[:target] = 1
    metric = join(metric,proplabel,on=:gene,kind=:anti)
    metric[:target] = 0
    metric = vcat(metric,tmp)
    metric[:bins] = bins(metric[:score],parse(Float64,split(bin_params,",")[1]),parse(Int64,split(bin_params,",")[2]))    
    for region in unique(metric[:set])
        binned_rates = by(metric[metric[:set] .== region,:],:bins,
                targetprop = :target => x-> sum(x)/count(!ismissing,x),
                gene_count = :target => y-> count(!ismissing,y))
        println("$(region):")
        println("Rho: $(corspearman(binned_rates[:bins],binned_rates[:targetprop]))")
        println(OneSampleZTest(atanh(corspearman(binned_rates[:bins],binned_rates[:targetprop])),
                             1, nrow(binned_rates)))
        @df binned_rates scatter!(:bins,:targetprop, 
                    xlabel = "$label",
                    ylabel = "Proportion in Target Set",
                    leg = true,
                    label = region,
                    smooth= true,
                    #marker = 10,
                    )
        savefig("$(label)_proportion_scatter.pdf")
        println("$region:")
        println(binned_rates)
    end
end

# calculates weighted mean and max for a set of regions, plots and stats by bins
function oddsRatio(metric::DataFrames.DataFrame,regions::Array{String,1},targets::String,label::String,bin_params::String)
    #label for odds ratio calculating
    proplabel = rename!(CSV.read(targets; delim='\t',allowmissing=:none,header=0),:Column1 => :gene)
    proplabel[:gene] = idOnly(proplabel[:gene])
    tmp = join(metric,proplabel,on=:gene,kind = :inner)
    tmp[:target] = 1
    metric = join(metric,proplabel,on=:gene,kind=:anti)
    metric[:target] = 0
    metric = vcat(metric,tmp)
    metric[:bins] = bins(metric[:score],parse(Float64,split(bin_params,",")[1]),parse(Int64,split(bin_params,",")[2]))  
    
    # calulate Odds Ratio of Target for each set vs. everything else
    r = metric[metric[:set] .== labelSet(regions[1]),:]
    CSV.write("dr_desert_genes.txt", r[r[:target] .== 1,:];delim='\t')
    r_target = nrow(r[r[:target] .== 1,:])
    matched = matchRates(r,metric[metric[:set] .!= labelSet(regions[1]),:])
    match_target = nrow(matched[matched[:target] .== 1,:])
    println("$(labelSet(regions[1])):")
    println(FisherExactTest(r_target,nrow(r)-r_target,match_target,nrow(matched)-match_target))
end

# splits genes in metric file into sets of regions, calls functions for various sets of stats
function comparison(metric_file::String,regions::Array{String,1},label::String,bin_params,targetlist,proportion,odds,excl,only)    
    metric = CSV.read(metric_file; delim='\t',allowmissing=:none,header=[:gene,:score],comment="#")    
    metric[:gene] = idOnly(metric[:gene])
    metric[:score] = metric[:score]*1. #make sure it's a float
    metric[:set] = "Other"
    for file_path in excl
        tmp = CSV.read(file_path; delim='\t',allowmissing=:none,header=[:gene])
        tmp[:gene] = idOnly(tmp[:gene])
        metric = join(metric,tmp,on=:gene,kind=:anti)
    end
    for file_path in regions
        for gene in CSV.read(file_path; delim='\t',header=[:gene])[:gene]
            metric[metric[:gene] .== split(gene,'.')[1],:set] = labelSet(file_path)
        end
    end
    if only
	metric = metric[metric[:set] .!= "Other",:]
    end
    if proportion
        prop(metric,targetlist,label,bin_params)
    elseif odds
        oddsRatio(metric,regions,targetlist,label,bin_params)
    else        
        genomeWide(metric,regions,label)
    end
end

function main()
    parsed_args = parse_commandline()
    if parsed_args["summary"]
        summaryStats(parsed_args["metric"], parsed_args["label"],parsed_args["exclude"])
    end
    if parsed_args["comparison"]
        comparison(parsed_args["metric"], parsed_args["regions"], parsed_args["label"], parsed_args["bins"], parsed_args["targetlist"], parsed_args["proportion"],
                parsed_args["odds"],parsed_args["exclude"],parsed_args["only"])
    end
end

main()

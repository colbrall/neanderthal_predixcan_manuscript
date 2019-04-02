# @author Laura Colbran
#
# calculates recombination rate for a given window and sets of files 
# assumes, if it's genome wide, that you have a file per chromosome (otherwise locations aren't mapped to a chromosome)
#
# require system-installed intersectBed
# julia 1.1

using ArgParse
#using Libz
using CSV
using DataFrames
using StatsBase
using HypothesisTests
using StatsPlots
using GLM
using RCall
pyplot(grid = false,leg = false)


# parses command-line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--distances","-d"
            help = "path to directory (if split by chromosome) or file (if combined) with genetic distances"
            arg_type = String
            required = true
        "--columns", "-c"
            help = "columns to consider from distance file; format: bp, centimorgan"
            arg_type = String
            required = true
        "--noheader" 
            help = "if there is no header in your distance file (default assumes yes)"
            action = :store_true
        "--summary"
            help = "if you want summary stats and a density plot of the recombination rate in your distance file"
            action = :store_true
        "--window", "-w"
            help = "window size (kb) for calculating recomb rate. Default 1. Required for summary."
            arg_type = Int
            default = 1
        "--write"
            help = "if you want to save a table of the windows or regions with calculated recomb rate"
            action = :store_true
        "--comparison"
            help = "if you want to compare recombination rates of regions"
            action = :store_true
        "--genomewide"
            help = "if your comparison whould be window-based and split the entire input map into categories"
            action = :store_true
        "--proportion"
            help = "if you are comparing proportions of a gene category across recombination rates"
            action = :store_true
        "--odds"
            help = "to calculate odds ratio of target set by region and matched by recombination rate"
            action = :store_true
        "--regions", "-r"
            nargs='*'
            help = "list of paths to bedfiles for regions you want to compare. required for comparison. Assumes they have a 4th column for the gene ID"
            arg_type = String
        "--exclude", "-e"
            help = "bed file for regions to exclude"
            arg_type = String
        "--targetlist", "-l"
            help = "list of gene IDs for your target set (usually DR genes)"
        "--model"
            help = "compares model stats and recombination rate and DR status"
            action = :store_true
        "--modfile","-m"
            help = "model summary file output from predixcan_model_summary.jl"
            arg_type = String
    end
    return parse_args(s)
end


# pulls label for set out of file name
function label(file_name::String)
    return split(basename(file_name),'.')[1]
end


# reads file for distances, returns user-specified columns as dataframe
function readDist(file_path::String,cols::Array{SubString{String},1},noheader)
    if endswith(file_path, ".gz")
        println("Woops: gzipped file handling not implemented yet")
        exit()
        #df = CSV.read(ZlibInflateInputStream(open(file_path));delim='\t',comment='#', 
        #            strings=:raw)
        #rename!(df,Symbol("#IDs")=>:gene_id)
    else
        if noheader
            df = CSV.read(file_path; delim=' ',header=0)
        else
            df = CSV.read(file_path; delim=' ')
        end
    end
    try 
        df = df[:,[parse(Int64,cols[1]),parse(Int64,cols[2])]]::DataFrames.DataFrame #filter df to columns in order specified by user
    catch crash
        if isa(crash, BoundsError)
            println("Woops: Nonexistent column specified!")
            exit()
        end
    end
    # any tweaks to formatting
    return df
end


# returns chromosome designation pulled from a file name, or returns an empty string if there's no designation
function chrom(file_name::String)
    if typeof(findfirst("chr",file_name)) == Nothing
        return ""
    end
    chr_pieces = "0123456789chrX" #define legal characters for chromosome designation
    chr_start = collect(findfirst("chr",file_name))[1]#identify place in file name where chr[] starts
    if chr_start+6 <= length(file_name)
        chr = file_name[chr_start:chr_start+6] ## some files do chr.[] for some reason
    else 
        chr = file_name[chr_start:length(file_name)]
    end
    chr = join(filter(x -> occursin(x,chr_pieces), collect(chr)),"")
    return chr
end


# count and return recombination rates over window
function recombRates(dist::DataFrames.DataFrame,window::Int64,chr::String)
    df = DataFrame(chr=String[], start=Int64[], stop=Int64[], rate = Float64[])
    first = 1::Int64
    last = 2::Int64
    while last <= nrow(dist)
        if dist[last,1]-dist[first,1] >= window
            if (abs((dist[last,1]-dist[first,1])-window) > abs((dist[last-1,1]-dist[first,1])-window)) && (last-1 != first) # if previous marker was closer to window 
                last -= 1
            end
            rate = (dist[last,2]-dist[first,2])/((dist[last,1]-dist[first,1])/1000000)::Float64 
            push!(df, [chr,dist[first,1], dist[last,1], rate])
            first = last
            last += 1
        else
            last += 1 
        end
    end
    return df
end


# calculates recombination rates across the genome for a given window size
# may write table of windows and rates
function calculateRates(dist_path::String,cols::Array{SubString{String},1},noheader,window::Int,write)
    if isfile(dist_path)
        dist = readDist(dist_path,cols,noheader)::DataFrames.DataFrame
        rates = recombRates(dist,window,chrom(dist_path))::DataFrames.DataFrame
    elseif isdir(dist_path)
        rates = DataFrame(chr=String[],start=Int64[],stop=Int64[],rate=Float64[])::DataFrames.DataFrame
        for item in readdir(dist_path)
            dist = readDist("$(realpath(dist_path))/$item",cols,noheader)::DataFrames.DataFrame
            rates = vcat(rates,recombRates(dist,window,chrom(item)))::DataFrames.DataFrame
        end
    end
    if write 
        CSV.write("recomb_rates_$(div(window,1000))kb.bed", rates, delim = '\t')
    end
    return rates
end


# bins a distribution of recombination rate, returns array of lower bounds of bins
function bins(dist::Array{Float64,1},max::Float64,num_bins::Int64)
    cut_points = collect(range(0,max,step=max/num_bins))
    bins = zeros(length(dist))
    for i in 1:length(dist)
        j = 1
        while j <= length(cut_points) && dist[i] > cut_points[j]
            j+=1
        end
        bins[i] = cut_points[j-1]
    end
    #println(Dict([(i,count(x->x==i,bins)) for i in unique(bins)]))
    return bins
end

# calculates and prints summary stats of a distribution, and plots the distribution
function summaryStats(distribution::Array{Float64,1},window::Int)
    describe(distribution)
    println("Std. Dev: $(mean_and_std(distribution)[2])")
    dist_plot = histogram(distribution, nbins = 100,
                    xlabel = "Recombination Rate",
                    ylabel = "Count",
                    color = [:black])
    savefig(dist_plot, "rates_distribution_$(div(window,1000))kb.pdf")
end


# calculates weighted mean and max for a set of regions, plots and stats
function meanMax(regions::Array{String,1},excl)
    region_rates = DataFrames.DataFrame(gene = String[], w_mean = Float64[], max = Float64[], set = String[])
    for file_path in regions
        tmp = DataFrames.DataFrame(gene = String[], length = Int64[], rate = Float64[], overlap = Int64[])
        if typeof(excl) != Nothing
            intersect = read(pipeline(`intersectBed -a $file_path -b $excl -v`,`intersectBed -a stdin -b tmp.bed -wo`),String)
        else
            intersect = read(`intersectBed -a $file_path -b tmp.bed -wo`,String)
        end
        for line in split(chomp(intersect),'\n') # intersect w/ map, summarize to one rate per row
            if line != ""
                entries = split(line,'\t')
                push!(tmp,[entries[4],parse(Int64,entries[3])-parse(Int64,entries[2]),parse(Float64,entries[8]),parse(Int64,entries[9])])
            end
        end
        for gene in unique(tmp[:gene])
            w_mean = 0.0
            t = tmp[tmp[:gene] .== gene,:]
            for row in 1:nrow(t)
                w_mean += t[row,:rate]*(t[row,:overlap]/(t[row,:length]))
            end
            push!(region_rates, [gene,w_mean, maximum(t[:,:rate]), label(file_path)])
        end
    end
    @df region_rates boxplot(:set,:w_mean, 
                    ylabel = "Recombination Rate (cM/Mb)",
                    color = [:blue],
                    #outliers = false
                    )
    savefig("rates_boxplot_w_mean.pdf")
    @df region_rates boxplot(:set,:max, 
                    xlabel = "Region Set",
                    ylabel = "Recombination Rate (cM/Mb)",
                    color = [:blue],
                    #outliers = false,
                    )
    savefig("rates_boxplot_max.pdf")
    if length(regions) > 2
        R"""
        suppressPackageStartupMessages(library(dplyr))
        suppressPackageStartupMessages(library(dunn.test))
        df = $region_rates
        group_by(df,set) %>% summarise(
            count = n(),
            mean = mean(w_mean, na.rm = TRUE),
            sd = sd(w_mean, na.rm = TRUE),
            median = median(w_mean, na.rm = TRUE),
            IQR = IQR(w_mean, na.rm = TRUE)
        )
        dunn.test(df$region_rates,df$set,method="bh")
        """
    elseif length(regions) == 2
        println("Weighted Mean:")
        println(MannWhitneyUTest(region_rates[region_rates[:set] .== label(regions[1]),:w_mean],region_rates[region_rates[:set] .== label(regions[2]),:w_mean]))
        println("Max:")
        println(MannWhitneyUTest(region_rates[region_rates[:set] .== label(regions[1]),:max],region_rates[region_rates[:set] .== label(regions[2]),:max]))
    end
end


# labels windows from map according  to overlap with regions, runs stats and plots on sets of windows
function genomeWide(regions::Array{String,1},excl)
    region_rates = DataFrames.DataFrame(rate = Float64[], set = String[])
    if typeof(excl) != Nothing
        run(pipeline(`intersectBed -a tmp.bed -b $excl -v`, stdout="t.bed"))
        run(`mv t.bed tmp.bed`)
    end
    for file_path in regions
        intersect = chomp(read(`intersectBed -a tmp.bed -b $file_path -u`,String))
        for line in split(intersect,'\n') # intersect w/ map, summarize to one rate per row
            if line != ""
                entries = split(line,'\t')
                push!(region_rates,[parse(Float64, entries[4]),label(file_path)])
            end
        end
        run(pipeline(`intersectBed -a tmp.bed -b $file_path -v`, stdout="t.bed"))
        run(`mv t.bed tmp.bed`)
    end
    for line in split(chomp(read(`cut -f 4 tmp.bed`,String)),'\n') # intersect w/ map, summarize to one rate per row
        push!(region_rates,[parse(Float64,line),"Other"])
    end
    @df region_rates boxplot(:set,:rate, 
                    xlabel = "Region Set",
                    ylabel = "Recombination Rate (cM/Mb)",
                    color = [:blue],
                    outliers = false
                    )
    savefig("rates_boxplot_genomewide.pdf")
    if length(regions) > 1
        println("kruskal wallis not implemented. run me in R.")
        CSV.write("region_rates_for_kw.txt",region_rates,delim = '\t')
        #KruskalWallisTest()
    elseif length(regions) == 1
        println(MannWhitneyUTest(region_rates[region_rates[:set] .== label(regions),:rate],regions_rates[region_rates[:set] .== "Other",:rate]))
    end
    # Kruskal-Wallis if >2 sets, MWU if 2
end


# pulls RR for a set of regions by intersecting with the map file generated by calculateRates()
function intersectRegions(file_path::String)
    df = DataFrames.DataFrame(gene = String[], length = Int64[], rate = Float64[], overlap = Int64[])
    intersect = read(`intersectBed -a $file_path -b tmp.bed -wo`,String)
    for line in split(chomp(intersect),'\n') # intersect w/ map, summarize to one rate per row
        if line != ""
            entries = split(line,'\t')
            push!(df,[entries[4],parse(Int64,entries[3])-parse(Int64,entries[2]),parse(Float64,entries[8]),parse(Int64,entries[9])])
        end
    end
    return df
end


# calculates weighted mean and max RR for windows across a gene
function rateSummary(df::DataFrames.DataFrame)
    w_mean = 0.0
    for row in 1:nrow(df)
        w_mean += df[row,:rate]*(df[row,:overlap]/(df[row,:length]))
    end
    return [w_mean, maximum(df[:,:rate])]
end


# returns a set of regions matched to target by r.r. using bins. should end up with same distribution, not same number
function matchRates(targets::DataFrames.DataFrame, candidates::DataFrames.DataFrame)
    bin_counts = sort!(by(targets,:mean_bins, gene_count = :target => y-> count(!ismissing,y)),:mean_bins)
    matched_set = DataFrames.DataFrame(gene = String[], w_mean = Float64[], max = Float64[], set = String[],target=Int64[],mean_bins=Float64[])
    # calculate number of times I can sample before I empty a bin
    n = nrow(candidates)
    for i in 1:nrow(bin_counts)
        reps = div(nrow(candidates[candidates[:mean_bins] .== bin_counts[i,:mean_bins],:]),bin_counts[i,:gene_count])
        if reps < n n = reps end
    end
    if n > 0
        for i in 1:nrow(bin_counts)
            bin = bin_counts[i,:mean_bins]
            count = bin_counts[i,:gene_count]
            bin_pick = candidates[candidates[:mean_bins] .== bin,:]
            for match in sample(collect(1:nrow(bin_pick)), n*count; replace = false)
                push!(matched_set,bin_pick[match,:])
            end
        end
        else #triggered if the target set isn't consistently smaller than the candidate set across bins
        for i in 1:nrow(bin_counts)
            bin = bin_counts[i,:mean_bins]
            count = bin_counts[i,:gene_count]
            bin_pick = candidates[candidates[:mean_bins] .== bin,:]
            if count <= nrow(bin_pick)
                for match in sample(collect(1:nrow(bin_pick)), count; replace = false)
                    push!(matched_set,bin_pick[match,:])
                end
            else
                count = nrow(bin_pick)
                bin_pick = targets[targets[:mean_bins] .== bin,:]
                for match in sample(collect(1:nrow(bin_pick)), count; replace = false)
                    push!(matched_set,bin_pick[match,:])
                end
            end
        end
    end
    return matched_set
end


# calculates weighted mean and max for a set of regions, plots and stats by bins
function prop(regions::Array{String,1},proplist::String,excl)
    region_rates = DataFrames.DataFrame(gene = String[], w_mean = Float64[], max = Float64[], set = String[])
    for file_path in regions
        tmp = intersectRegions(file_path,excl)
        for gene in unique(tmp[:gene])
            r = rateSummary(tmp[tmp[:gene] .== gene,:])
            push!(region_rates, [gene, r[1],r[2], label(file_path)])
        end 
    end
    #label for proportion calculating
    proplabel = rename!(CSV.read(proplist; delim='\t',header=0),:Column1 => :gene)
    tmp = join(region_rates,proplabel,on=:gene,kind = :inner)
    tmp[:target] = 1
    region_rates = join(region_rates,proplabel,on=:gene,kind=:anti)
    region_rates[:target] = 0
    region_rates = vcat(region_rates,tmp)
    #describe(region_rates[region_rates[:set] .== "desert_2Mb_ldexpand",:w_mean])
    region_rates[:mean_bins] = bins(region_rates[:w_mean],2.4,30)    
    #region_rates[:max_bins] = bins(region_rates[:max])
    for region in unique(region_rates[:set])
        binned_rates = by(region_rates[region_rates[:set] .== region,:],:mean_bins,
                targetprop = :target => x-> sum(x)/count(!ismissing,x),
                gene_count = :target => y-> count(!ismissing,y))
        println("$(region):")
        println("Spearman rho: $(corspearman(binned_rates[:mean_bins],binned_rates[:targetprop]))")
        println(OneSampleZTest(atanh(corspearman(binned_rates[:mean_bins],binned_rates[:targetprop])),
                             1, nrow(binned_rates)))
        @df binned_rates scatter!(:mean_bins,:targetprop, 
                    xlabel = "Recombination Rate (cM/Mb)",
                    ylabel = "Proportion in Target Set",
                    leg = true,
                    label = region,
                    smooth= true,
                    )
        savefig("rates_proportion_scatter.pdf")
        println("$region:")
        println(binned_rates)
    end
end


# calculates weighted mean and max for a set of regions, plots and stats by bins
function oddsRatio(regions::Array{String,1},targets::String,excl)
    region_rates = DataFrames.DataFrame(gene = String[], w_mean = Float64[], max = Float64[], set = String[])
    for file_path in regions
        tmp = intersectRegions(file_path)
        for gene in unique(tmp[:gene])
            r = rateSummary(tmp[tmp[:gene] .== gene,:])
            push!(region_rates, [gene, r[1],r[2], label(file_path)])
        end 
    end
    if typeof(excl) != Nothing
        region_rates = join(region_rates,CSV.read(excl; delim='\t',allowmissing=:none,header=[:gene]),on=:gene,kind=:anti)        
    end
    #label for odds ratio calculating
    proplabel = rename!(CSV.read(targets; delim='\t',header=0),:Column1 => :gene)
    tmp = join(region_rates,proplabel,on=:gene,kind = :inner)
    tmp[:target] = 1
    region_rates = join(region_rates,proplabel,on=:gene,kind=:anti)
    region_rates[:target] = 0
    region_rates = vcat(region_rates,tmp)
    region_rates[:mean_bins] = bins(region_rates[:w_mean],2.4,30)   
    
    # calulate Odds Ratio of Target for each set vs. everything else
    r = region_rates[region_rates[:set] .== label(regions[1]),:]
    r_target = nrow(r[r[:target] .== 1,:])
    matched = matchRates(r,region_rates[region_rates[:set] .!= label(regions[1]),:])
    match_target = nrow(matched[matched[:target] .== 1,:])
    println("$(label(regions[1])):")
    println(FisherExactTest(r_target,nrow(r)-r_target,match_target,nrow(matched)-match_target))
end


# plots and runs stats to compare recombination rates between sets of genomic regions
function comparison(map::DataFrames.DataFrame, regions::Array{String,1},genomewide,proportion,targetlist,odds,excl)
    if genomewide
        genomeWide(regions,excl)
    elseif proportion
        prop(regions,targetlist,excl)
    elseif odds
        oddsRatio(regions,targetlist,excl)
    else 
        meanMax(regions,excl)
    end
end

        
# plots and compares RR to model performance and #SNPs, and proportion of DR
function modelSummary(region_path::Array{String,1},model_file::String,excl)
    region_rates = DataFrames.DataFrame(gene_id = String[], w_meanRR = Float64[], maxRR = Float64[], set = String[])
    for file_path in region_path
        tmp = intersectRegions(file_path,excl)
        for gene in unique(tmp[:gene])
            r = rateSummary(tmp[tmp[:gene] .== gene,:])
            push!(region_rates, [gene, r[1],r[2], label(file_path)])
        end 
    end
    models = CSV.read(model_file; delim='\t')
    for i in 1:nrow(models)
        models[i,:gene_id] = split(models[i,:gene_id],".")[1]
    end
    models = by(models, :gene_id, meanR2 = :R2 => x-> sum(x)/count(!ismissing,x), maxR2 = :R2 => maximum,
                meanSNPs = :NumSNPs => x-> sum(x)/count(!ismissing,x), maxSNPs = :NumSNPs => maximum)
    region_rates = join(region_rates,models,on=:gene_id,kind=:inner)
    println("mean RR vs mean R2")
    println("Spearman rho: $(corspearman(region_rates[:w_meanRR],region_rates[:meanR2]))")
        println(OneSampleZTest(atanh(corspearman(region_rates[:w_meanRR],region_rates[:meanR2])),
                             1, nrow(region_rates)))
    println("mean RR vs max R2")
    println("Spearman rho: $(corspearman(region_rates[:w_meanRR],region_rates[:maxR2]))")
        println(OneSampleZTest(atanh(corspearman(region_rates[:w_meanRR],region_rates[:maxR2])),
                             1, nrow(region_rates)))
    println("mean RR vs mean Num SNPs")
    println("Spearman rho: $(corspearman(region_rates[:w_meanRR],region_rates[:meanSNPs]))")
        println(OneSampleZTest(atanh(corspearman(region_rates[:maxRR],region_rates[:meanSNPs])),
                             1, nrow(region_rates)))
    println("mean RR vs max Num SNPs")
    println("Spearman rho: $(corspearman(region_rates[:w_meanRR],region_rates[:maxSNPs]))")
        println(OneSampleZTest(atanh(corspearman(region_rates[:maxRR],region_rates[:maxSNPs])),
                             1, nrow(region_rates)))
end


function main()
    parsed_args = parse_commandline()
    map = calculateRates(parsed_args["distances"],split(parsed_args["columns"],','),parsed_args["noheader"],
        parsed_args["window"]*1000,parsed_args["write"])::DataFrames.DataFrame
    if parsed_args["summary"]
        summaryStats(map[:,4],parsed_args["window"]*1000)
    elseif parsed_args["comparison"]
        CSV.write("tmp.bed", map, delim = '\t', writeheader=false)
        comparison(map,parsed_args["regions"],parsed_args["genomewide"],parsed_args["proportion"],
            parsed_args["targetlist"],parsed_args["odds"],parsed_args["exclude"])
        run(`rm tmp.bed`)
    elseif parsed_args["model"]
        CSV.write("tmp.bed", map, delim = '\t', writeheader=false)
        modelSummary(parsed_args["regions"],parsed_args["modfile"],parsed_args["exclude"])
        run(`rm tmp.bed`)
    end
end

main()

# intro_individuals.jl
#
# @author Laura Colbran
#
# for each gene, identify individuals with an introgressed haplotype within 1Mb of the gene. works from gene_regions_neanSNPs.txt
# 
# requires system-installed bcftools
# julia1.1

using ArgParse
using CSV
using DataFrames
using GZip

# parses command-line arguments
function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--genotypes","-g"
            help = "path to vcf containing genotypes for individuals"
            arg_type = String
            required = true
        "--snps", "-s"
            help = "path to file with SNPs to filter. Assumes structured like gene_regions_neanSNPs.txt"
            arg_type = String
            required = true
        "--output_file","-o"
            help = "output file name"
            arg_type = String
            required = true
    end
    return parse_args(s)
end


function main()
    parsed_args = parseCommandLine()
    
    GZip.open(parsed_args["snps"]) do f
        prev = ""
        tmp = DataFrames.DataFrame(gene = String[],chr=Int64[],pos=Int64[],allele=String[])
        open(parsed_args["output_file"],"w") do g
            write(g,"#gene\tgtex_id\tnum_intro_snps\n")
            for line in eachline(f)
                if !startswith(line,"ENS") continue end #skip headers
                l = split(line, '\t')
                curr = l[1]
                if prev == "" prev = curr end
                if curr == prev 
                    push!(tmp, [curr,parse(Int64,split(l[2],'r')[2]),parse(Int64,l[3]),l[4]])
                    continue
                end
                if nrow(tmp) == 0 continue end
                CSV.write("tmp.txt", tmp[:,[:chr,:pos]];delim='\t',writeheader=false)
                cmd = `bcftools view $(parsed_args["genotypes"]) -R tmp.txt`
                id_dict = Dict{String,Int64}()
                ind_ids = Array{String,1}
                for snp in split(chomp(read(cmd,String)),'\n') # subset VCF to correct SNPs for this gene
                    if startswith(snp,"##") continue end #skip vcf header
                    if startswith(snp,"#")
                        ind_ids = split(snp,'\t')[10:end] #pull individual IDs
                        for i in 1:length(ind_ids)
                            id_dict[ind_ids[i]] = 0
                        end
                        continue
                    end
                    intro_allele = tmp[tmp[:pos] .== split(snp,'\t')[2],:allele] #id allele I care about
                    a = "alt"
                    if intro_allele == split(snp,'\t')[4] a = "ref" end #in case SNP is backwards
                    genotypes = split(snp,'\t')[10:end]
                    for i in 1:length(genotypes) # identify individuals with allele
                        dose = round(parse(Float64,split(genotypes[i],":")[3]))
                        if a == "ref" dose = 2-dose end #swap if backwards
                        if dose > 0 
                            id_dict[ind_ids[i]] = id_dict[ind_ids[i]]+ 1 # count SNPs
                        end
                    end
                end
                for ind in keys(id_dict)
                    individual = join(split(ind,"-")[1:2],"-")
                    write(g,"$(curr)\t$(individual)\t$(id_dict[ind])\n")
                end
                prev = curr
                tmp = DataFrames.DataFrame(gene = String[],chr=Int64[],pos=Int64[],allele=String[])
                push!(tmp, [curr,parse(Int64,split(l[2],'r')[2]),parse(Int64,l[3]),l[4]])
            end
        end
    end
end

main()
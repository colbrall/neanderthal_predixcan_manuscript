# gene_overlap.jl
# @author Laura Colbran
#
# assumes bed files are sorted
# usage: julia gene_overlap.jl ARGS
# look above each function for specific ARGS

using DataFrames
using GenomicFeatures

function readDF(f_path::String)
  println("Reading file $f_path as DataFrame.......")
  if ispath(f_path)
    df = readtable(f_path,separator='\t', normalizenames=false, allowcomments=true,header=false)
  else
    df = "NA"
  end
  return df
end

function geneDict(path::String)
    const dict = Dict{String,Array{String,1}}()
    open(path) do f
      for line in eachline(f)
          dict[split(chomp(line),"\t")[4]] = []
      end
    end
    return dict
end

# filters list of gene locations to just those we care about, calculates 2-Mb region for overlap, writes to tmp.bed
# ARGS: path/to/gene/locs/bed path/to/target/genes
function geneRegions(gene_bed::String,target_list::String,window::Int64)
    println("Finding target gene locations...")
    gene_df = readDF(realpath(gene_bed))
    for i in 1:nrow(gene_df)
        gene_df[i,1] = "$(gene_df[i,1][4:length(gene_df[i,1])])" #removes the chr from the name
        gene_df[i,4] = split(split(split(gene_df[i,4],";")[1]," ")[2],".")[1] #pull out gene_id
    end
    target_df = readDF(realpath(target_list))
    target_df[:x4] = target_df[:x1] #make gene_id column names match
    delete!(target_df,:x1)
    gene_df = join(gene_df, target_df, on=:x4, kind=:inner)
    genes = unique(target_df[:x4])
    println("\nFound location for $(length(unique(gene_df[:x4])))/$(length(genes)) target genes...")
    open("tmp.bed","w") do f
        for i in 1:length(genes)
            try #works for all genes whose locations were found
                tmp = gene_df[gene_df[:,:x4] .== genes[i],:]
                m = minimum(tmp[:,:x2])-window
                if (m < 0) m = 0 end
                write(f,"$(tmp[1,1])\t$(m)\t$(maximum(tmp[:,:x3])+window)\t$(genes[i])\n")
            end
        end
    end
    run(pipeline(`sortBed -i tmp.bed`,"t.bed"))
    run(`mv t.bed tmp.bed`)
end

# ARGS: path/to/gene/locs/bed path/to/target/genes interesting/snps/bed
function countSNPs(gene_regions::String,snps::String)
    println("\nIntersecting SNPs with genes...")
    gene_dict = geneDict(gene_regions)::Dict{String,Array{String,1}}
    for (gene,snp) in eachoverlap(open(BED.Reader,gene_regions),open(BED.Reader,snps))
        push!(gene_dict[BED.name(metadata(gene))], BED.name(metadata(snp)))
    end
    open("gene-snp_intersection_count.txt","w") do f
        write(f,"#gene\tsnp_count\n")
        for key in keys(gene_dict)
            write(f,"$(key)\t$(length(gene_dict[key]))\n")
        end
    end
    run(`rm tmp.bed`) #remove tmp file if using one
end

# pulls out genes whose regions overlap certain haplotypes
# assumes both are bed files
# ARGS: path/to/gene/locs/bed path/to/target/genes interesting/haplotypes/bed
function hapGenes(gene_regions::String,haps::String)
  println("\nIntersecting Genes with haplotypes...")
  f = open("tmpgenes.txt","w")
  for (gene,hap) in eachoverlap(open(BED.Reader,gene_regions),open(BED.Reader,haps))
      write(f,"$(BED.name(metadata(gene)))\n")
  end
  close(f)
  run(pipeline(`sort -u tmpgenes.txt`,"hap_genes.txt"))
  run(`rm tmpgenes.txt`)
end

# pulls out SNPs within certain haplotypes (based solely on location, not LD)
# assumes both are bed files
# ARGS: snps haps
# or: ARGS: path/to/gene/locs/bed path/to/target/genes interesting/haplotypes/bed interesting/snps/bed
function hapSNPs(snps::String,haps::String)
  println("\nIntersecting SNPs with haplotypes...")
  f = open("tmpsnps.bed","w")
  for (snp,hap) in eachoverlap(open(BED.Reader,snps),open(BED.Reader,haps))
      write(f,"$(BED.chrom(metadata(snp)))\t$(BED.chromstart(metadata(snp)))\t$(BED.chromend(metadata(snp)))\t$(BED.name(metadata(snp)))\n")
  end
  close(f)
  run(pipeline(`sort -u tmpsnps.bed`, "hap_snps.bed"))
  run(`rm tmpsnps.bed`)
end

# filters list of gene locations to just those we care about, gets location, writes to BED
# ARGS: path/to/gene/locs/bed path/to/target/genes
function geneBed(gene_bed::String,target_list::String)
    println("Finding target gene locations...")
    gene_df = readDF(realpath(gene_bed))
    for i in 1:nrow(gene_df)
        gene_df[i,1] = "$(gene_df[i,1][4:length(gene_df[i,1])])" #removes the chr from the name
        gene_df[i,4] = split(split(split(gene_df[i,4],";")[1]," ")[2],".")[1] #pull out gene_id
    end
    target_df = readDF(realpath(target_list))
    target_df[:x4] = target_df[:x1] #make gene_id column names match
    delete!(target_df,:x1)
    gene_df = join(gene_df, target_df, on=:x4, kind=:inner)
    genes = unique(target_df[:x4])
    println("\nFound location for $(length(unique(gene_df[:x4])))/$(length(genes)) target genes...")
    open("genes.bed","w") do f
        for i in 1:length(genes)
            try #works for all genes whose locations were found
                tmp = gene_df[gene_df[:,:x4] .== genes[i],:]
                write(f,"$(tmp[1,1])\t$(minimum(tmp[:,:x2]))\t$(maximum(tmp[:,:x3]))\t$(genes[i])\n")
            end
        end
    end
    run(pipeline(`sortBed -i genes.bed`,"t.bed"))
    run(`mv t.bed genes.bed`)
end


function main()
  args = ARGS[:,1]::Array{String,1}
  #geneRegions(args[1],args[2],1000000) #makes tmp.bed used by other functions
  #countSNPs("tmp.bed",args[3]) #saves file called gene-snp_intersection_count.txt
  #hapGenes("tmp.bed",args[3]) #saves file called hap_genes.txt
  #hapSNPs(args[4],args[3]) #saves file called hap_snps.bed
  geneBed(args[1],args[2]) #makes bed file of genes

  #run(`rm tmp.bed`)
end

main()

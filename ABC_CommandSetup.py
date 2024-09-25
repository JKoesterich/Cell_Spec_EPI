import os 
import os.path
import sys 
import commentjson

# read the input from the json file and set up the common paths 
jfil = open(sys.argv[1], 'r').read()
param = commentjson.loads(jfil)
param["OutDir"] = param["OutDir"].rstrip("/")
param["ABC_loc"] = param["ABC_loc"].rstrip("/")

crpar = " --regions_blocklist " + param["hg19_ref"]["BlockList"] + " --regions_includelist " + param["hg19_ref"]["IncludeList"] + " --peakExtendFromSummit " + param["ABC_Param"]["ExtendPeak"] + " --nStrongestPeaks " + param["ABC_Param"]["NStrongPeak"]

outfol = param["OutDir"] + "/" + param["JobName"] + "_ABCoutput/"
print("mkdir " + outfol)
print("source miniconda3/bin/deactivate")
print("source miniconda3/bin/activate miniconda3/envs/final-abc-env/")
print("module load samtools")
print("module load bedtools2")
# For each of the replicates we need to sort and call candidate regions from them
for npfi in range(0, len(param["NarrowPeak"])):
    outnp = outfol + os.path.basename(param["NarrowPeak"][npfi]) + ".sorted"
    npsor = "bedtools sort -faidx " + param["hg19_ref"]["GenomeLen"] + " -i " + param["NarrowPeak"][npfi] + " > " + outnp
    #npsor = "bedtools sort -faidx " + param["hg19_ref"]["GenomeLen"] + " -i " + param["NarrowPeak"][npfi] + " | cat - > " + outnp
    npcan = "python " + param["ABC_loc"] + "/src/makeCandidateRegions.py --narrowPeak " + outnp + " --bam " + param["ShiftBam"][npfi] + " --outDir " + outfol + " --chrom_sizes " + param["hg19_ref"]["GenomeLen"] + crpar
    print(npsor)
    print(npcan)
# Need to combine the different replicates, using if statements to determine how to combine them based on the number of replicates
if(len(param["NarrowPeak"]) == 2):
    np1 = outfol + os.path.basename(param["NarrowPeak"][0]) + ".sorted.candidateRegions.bed"
    np2 = outfol + os.path.basename(param["NarrowPeak"][1]) + ".sorted.candidateRegions.bed"
    inter = "bedtools intersect -a " + np1 + " -b " + np2 + " > " + outfol + param["JobName"] + "_CandRegionRepOverlaps.bed"
    insr = "bedtools sort -faidx " + param["hg19_ref"]["GenomeLen"] + " -i " + outfol + param["JobName"] + "_CandRegionRepOverlaps.bed > " + outfol + param["JobName"] + "_CandRegionRepOverlaps_sorted.bed"
    print(inter)
    print(insr)
elif(len(param["NarrowPeak"]) == 1):
    insr = "bedtools sort -faidx " + param["hg19_ref"]["GenomeLen"] + " -i " + os.path.basename(param["NarrowPeak"][0]) + ".sorted.candidateRegions.bed > " + outfol + param["JobName"] + "_CandRegionRepOverlaps_sorted.bed"
    print(insr)
else:
    in1 = outfol + os.path.basename(param["NarrowPeak"][0]) + ".sorted.candidateRegions.bed"
    for npcri in range(1, len(param["NarrowPeak"])):
        in2 = outfol + os.path.basename(param["NarrowPeak"][npcri]) + ".sorted.candidateRegions.bed"
        inter = "bedtools intersect -a " + in1 + " -b " + in2 + " > " + outfol + param["JobName"] + "_CandRegionInterReps.bed"
        interin = "cp " + outfol + param["JobName"] + "_CandRegionInterReps.bed " + outfol + param["JobName"] + "_CandRegionInterRepToGo.bed"
        in1 = outfol + param["JobName"] + "_CandRegionInterRepToGo.bed"
        print(inter)
        print(interin)
    inten = "mv " + outfol + param["JobName"] + "_CandRegionInterRepToGo.bed " + outfol + param["JobName"] + "_CandRegionRepOverlaps.bed"
    intrm = "rm " + outfol + param["JobName"] + "_CandRegionInterReps.bed"
    insr = "bedtools sort -faidx " + param["hg19_ref"]["GenomeLen"] + " -i " + outfol + param["JobName"] + "_CandRegionRepOverlaps.bed > " + outfol + param["JobName"] + "_CandRegionRepOverlaps_sorted.bed"
    print(inten)
    print(intrm)
    print(insr)

cr500 = "awk '{a=($2+$3)/2; printf($1"'"\\t%1.0f\\t%1.0f\\n"'",a-250,a+250)}' " + outfol + param["JobName"] + "_CandRegionRepOverlaps_sorted.bed" + " | bedtools sort -i stdin | bedtools merge -i stdin | " + "awk '{a=($2+$3)/2; printf($1"'"\\t%1.0f\\t%1.0f\\n"'",a-250,a+250)}' > " + outfol + param["JobName"] + "_CandRegionRepOverlaps_sorted_500bp.bed"
print(cr500)
# Running the neighborhood command 
if(len(param["ShiftBam"]) > 1):
    atdn = ' "' + '","'.join(param["ShiftBam"]) + '" '
else:
    atdn = ' "' + param["ShiftBam"][0] + '" '
if(len(param["H3K27ac"]) > 1):
    h3 = ' "' + '","'.join(param["H3K27ac"]) + '" '
else:
    h3 = ' "' + param["H3K27ac"][0] + '" '
if(param["TPM_Table"] == ""):
    gene_exp = ""
else:
    gene_exp = " --expression_table " + param["TPM_Table"] 

adh = " --" + param["ACtype"] + atdn + "--H3K27ac" + h3
rnp = "--chrom_sizes " + param["hg19_ref"]["GenomeLen"] + " --ubiquitously_expressed_genes " + param["hg19_ref"]["Ubiq_Genes"] + " --genes " + param["hg19_ref"]["GeneList"] + " --cellType " + param["JobName"] + " --outdir " + outfol
#rnnh = "python " + param["ABC_loc"] + "/src/run.neighborhoods.py --candidate_enhancer_regions " + outfol + param["JobName"] + "_CandRegionRepOverlaps_sorted.bed" + adh + rnp + gene_exp
rnnh = "python " + param["ABC_loc"] + "/src/run.neighborhoods.py --candidate_enhancer_regions " + outfol + param["JobName"] + "_CandRegionRepOverlaps_sorted_500bp.bed" + adh + rnp + gene_exp
print(rnnh)

# Since there can be many ways that the HiC data can be in I have to assume that it is correctly formatted having each chr separated into their own .gz file inside their own chr labelled folder

if(param["HiC"][0] == ""):
    predpr = " --score_column powerlaw.Score"
else:
    predpr = " --hic_resolution " + param["ABC_Param"]["hic_res"] + " --HiCdir " + param["HiC"][0] + " --hic_type " + param["HiC"][1]
    
pred = "python " + param["ABC_loc"] + "/src/predict.py --enhancers " + outfol + "EnhancerList.txt --genes " + outfol + "GeneList.txt --threshold " + param["ABC_Param"]["predict_thresh"] 
predo = pred + predpr + " --chrom_sizes " + param["hg19_ref"]["GenomeLen"] + " --make_all_putative --cellType " + param["JobName"] + " --outdir " + outfol
print(predo) 

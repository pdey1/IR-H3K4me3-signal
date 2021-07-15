###############################################################################################################
# Calculate mRNA exonic nucleosome densities
#
# Pinki Dey, July 2021
###############################################################################################################

# Prerequisite: Download reference genome data.
# - Get the gtf file of evidence-based annotation of the human genome (GRCh38), version 34 (Ensembl 100) - long non-coding RNAs
# - Get all exons from UCSC Table Browser UCSC refSeq table: hg38_refSeq_exons.bed
# - Get all introns+80 bp in each end from UCSC Table Browser UCSC refSeq table: hg38_refSeq_introns+80bp.bed

#########Extract the internal mRNA exons################

awk '{if($3=="exon"){print $0}}' gencode.v34.long_noncoding_RNAs.gtf > exons.gtf
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$4,$5,"NA","0",$7}' exons.gtf | sort -k1,1 -k2,2n > lnc_exons.bed

awk -F'\t' 'NR==FNR{c[$2]++;next};c[$2]' lnc_exons.bed hg38_refSeq_exons.bed > common_exons
sort -k1,1 -k2,2n common_exons > common_exons_sorted

#Extract all the exons excluding the common exons

awk 'NR==FNR{a[$0];next}!($0 in a)' common_exons_sorted hg38_refSeq_exons.bed > mRNA_exons

#Omit the fragmented parts of the genomic assembly

grep -v alt mRNA_exons | grep -v fix | grep -v random | grep -v chrUn > mRNA_clean.bed

grep -v alt hg38_refSeq_introns+80bp.bed | grep -v fix | grep -v random | grep -v chrUn > hg38_refSeq_introns+80bp_clean.bed

# Create a collection of first and last exons

grep -B1 exon_0_ mRNA_clean.bed | grep -v -- "^--$" > mRNA_first_and_last_exons.bed
tail -1 mRNA_clean.bed >> mRNA_first_and_last_exons.bed

# Collapse to unique introns and exons

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$6}' hg38_refSeq_introns+80bp_clean.bed | sort -k2,2n -k3,3n | uniq | sort > hg38_refSeq_introns+80bp_unique.bed

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$6}' mRNA_first_and_last_exons.bed | sort -k2,2n -k3,3n | uniq | sort > mRNA_first_and_last_exons_unique.bed

# Include exonic ends of 80 bp and intronic ends of 500 bp. Create different files for 3 prime and 5 prime splice sites.

awk 'BEGIN {FS="\t"; OFS="\t"} {if ($4=="+") {print $1,$2,$2+580,$4} else {print $1,$3-580,$3,$4} }' hg38_refSeq_introns+80bp_unique.bed | uniq > hg38_refSeq_3prSS_exonic80bp_intronic500bp_unique.bed

awk 'BEGIN {FS="\t"; OFS="\t"} {if ($4=="-") {print $1,$2,$2+580,$4} else {print $1,$3-580,$3,$4} }' hg38_refSeq_introns+80bp_unique.bed | uniq > hg38_refSeq_5prSS_exonic80bp_intronic500bp_unique.bed

# Convert the .bed files to six-field format 

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,"NA","0",$4}' hg38_refSeq_5prSS_exonic80bp_intronic500bp_unique.bed > hg38_refSeq_5prSS_exonic80bp_intronic500bp_unique_6col.bed

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,"NA","0",$4}' hg38_refSeq_3prSS_exonic80bp_intronic500bp_unique.bed > hg38_refSeq_3prSS_exonic80bp_intronic500bp_unique_6col.bed

# Prerequisite: Install bedtools 
# Remove all regions overlapping first or last exons with bedtools

bedtools intersect -v -a hg38_refSeq_3prSS_exonic80bp_intronic500bp_unique_6col.bed -b mRNA_first_and_last_exons_unique.bed > hg38_refSeq_3prSS_exonic80bp_intronic500bp_unique_internal.bed

bedtools intersect -v -a hg38_refSeq_5prSS_exonic80bp_intronic500bp_unique_6col.bed -b mRNA_first_and_last_exons_unique.bed > hg38_refSeq_5prSS_exonic80bp_intronic500bp_unique_internal.bed

# Prerequisite
# - download and Install bedtools 
# - download and install LiftOver.
# - download bedGraphToBigWig and bedItemOverlapCount from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/
# - download chromInfo from http://hgdownload.cse.ucsc.edu	/goldenpath/hg38/database/
# - download bwtool 
#
################################################################################################################
#
# Download the nucleosomal libraries 
# - Human total nucleosome library is obtained from NCBI Short Read Archive (SRA) accession number SRP000105.
# - Histone modified nucleosome libraries are obatined from https://dir.nhlbi.nih.gov/papers/lmi/epigenomes/hgtcell.aspx
#For each of the libraries, retrieve nucleosome density wrt exon position strand-specifically.
#
# Here, showing only one example library (H3K4me3 from Barski et al., 2007,
# available in https://dir.nhlbi.nih.gov/papers/lmi/epigenomes/hgtcell.aspx

# Notice that first it had to be converted with UCSC liftOver from hg18->hg19 and hg19->hg38
# The input to LiftOver Utility is obatined from http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz 
# and http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

sort -k1,1 -k2,2n H3K4me3.bed > Histone_sorted.bed

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,"NA","0",$6}' Histone_sorted.bed > Histone_sorted_6col.bed

liftOver  Histone_sorted_6col.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed
liftOver output.bed hg19ToHg38.over.chain.gz Histone_sorted_hg38.bed unlifted.bed

#Extend all the uniquely mapped short-read sequences of the nucleosome library to the expected 150 basepair (bp)

awk 'BEGIN{FS="\t";OFS="\t"} {if ($6=="+") {print $1,$2,$2+150,$4, $5, $6} else {print $1,$3-150,$3,$4,$5,$6}}' Histone_sorted_hg38.bed > Histone_sorted_ext.bed

grep -v alt Histone_sorted_ext.bed | grep -v fix | grep -v random | grep -v chrUn > Histone_sorted_ext_clean.bed

###Perform the nucleosome density analysis.


export f=Histone_sorted_ext_clean
	
	awk '{if($6=="+") print}' ${f}.bed | sort -k1,1 | bedItemOverlapCount hg38 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n > ${f}.plus.bedGraph
	awk '{if($6=="-") print}' ${f}.bed | sort -k1,1 | bedItemOverlapCount hg38 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > ${f}.minus.bedGraph
	
	
bedGraphToBigWig ${f}.plus.bedGraph chromInfo.txt ${f}.plus.bw
bedGraphToBigWig ${f}.minus.bedGraph chromInfo.txt ${f}.minus.bw

# Create strand- and splice site-specific matrices of nucleosome densities at exon/intron boundaries

export f=Histone_sorted_ext_clean

	export ref=hg38_refSeq_5prSS_exonic80bp_intronic500bp_unique_internal
	bwtool extract bed ${ref}.bed ${f}.plus.bw ${f}.5prSS.plus.matrix
	bwtool extract bed ${ref}.bed ${f}.minus.bw ${f}.5prSS.minus.matrix

	export ref=hg38_refSeq_3prSS_exonic80bp_intronic500bp_unique_internal
	bwtool extract bed ${ref}.bed ${f}.plus.bw ${f}.3prSS.plus.matrix
	bwtool extract bed ${ref}.bed ${f}.minus.bw ${f}.3prSS.minus.matrix

	cat ${f}.5prSS.plus.matrix ${f}.5prSS.minus.matrix > ${f}.5prSS.matrix
	cat ${f}.3prSS.plus.matrix ${f}.3prSS.minus.matrix > ${f}.3prSS.matrix

# Calculate the average nucleosomal density for upstram and downstream of the exons.

awk '{$1=$2=$3=$4=$5=$6=$7=""; print $0}' Histone_sorted_ext_clean.5prSS.matrix | tr -s 'NA' '0' | tr -d "-" | awk -F, '{for(i=1;i<=NF;i++)a[i]+=$i} END{for(i=1;i<=NF;i++)printf "%d%s", a[i], (i==NF?"\n":",")}'  > sum5s

awk '{$1=$2=$3=$4=$5=$6=$7=""; print $0}' Histone_sorted_ext_clean.3prSS.matrix | tr -s 'NA' '0' | tr -d "-" | awk -F, '{for(i=1;i<=NF;i++)a[i]+=$i} END{for(i=1;i<=NF;i++)printf "%d%s", a[i], (i==NF?"\n":",")}'  > sum3s

tr -s ',' '\n' < sum5s > sum5s_mod
tr -s ',' '\n' < sum3s > sum3s_mod

awk 'BEGIN {for (i=50; i<= 630; i++) print i}' > 3prime_bp_number
awk 'BEGIN {for (i=1; i<=580; i++) print -i}' | sort -n > 5prime_bp_number

paste 3prime_bp_number sum3s_mod > sum3s_mod1
paste 5prime_bp_number sum5s_mod > sum5s_mod1

wc -l Histone_sorted_ext_clean.5prSS.matrix > x1
wc -l Histone_sorted_ext_clean.3prSS.matrix > x2

#Get the avergae by dividing the sum by the total number of regions obtained from the above command.
awk '{print $1"\t"$2/x1}' sum5s_mod1 > sum5s_avg 
awk '{print $1"\t"$2/x2}' sum3s_mod1 > sum3s_avg

# - visualize 3pr SS on the left and 5pr SS on the right in the same chart






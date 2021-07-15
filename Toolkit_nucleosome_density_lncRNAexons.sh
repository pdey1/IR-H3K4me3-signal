###############################################################################################################
# Calculate lncRNA exonic nucleosome densities
#
# Pinki Dey, July 2021
###############################################################################################################

# Prerequisite: Download reference genome data from Gencode
# - Get the gtf file of evidence-based annotation of the human genome (GRCh38), version 34 (Ensembl 100) - long non-coding RNAs

##### Calculate the internal exons from the gtf file after removing the coordinates corresponding to the UTR regions.

awk '{if($3=="exon"){print $0}}' gencode.v34.long_noncoding_RNAs.gtf > exons.gtf
grep -v "alternative_5_UTR" exons.gtf | grep -v "alternative_3_UTR"> exons_no_UTR.gtf

#Print out the gene names, coordinates and exons number for calcultating the internal exons.

awk '{print $1"\t"$4"\t"$5"\t"$7"\t"$20"\t"$22}' exons_no_UTR.gtf > exons_no_UTR_exon_no.gtf # Remove the ; from the file for calculations.

#Count the exons
awk '{a[$6]++;} END{for(i in a) print a[i]" "i}' exons_no_UTR_exon_no.gtf | awk '($1!=2){print $1"\t"$2 }' | awk '($1!=1){print $1"\t"$2 }' > exons_no_UTR_exon_count.gtf 

#Get the internal exons
awk 'BEGIN { while(getline <"exons_no_UTR_exon_count.gtf") id[$2]=1; } id[$6]' exons_no_UTR_exon_no.gtf > exons_no_UTR_exon_positions.gtf

awk 'p{print $6-p}{p=$6}' exons_no_UTR_exon_positions.gtf > exons_no_UTR_exon_positions_diff.gtf

paste exons_no_UTR_exon_positions.gtf exons_no_UTR_exon_positions_diff.gtf > exons_positions.gtf ##Delete the last rwo rows as it had only two exons -first and last

awk '$7> 0 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' exons_positions.gtf | awk '($6 != 1){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > internal_exons.gtf

# Convert the internal exons from gtf to bed format
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,"NA","0",$4}' internal_exons.gtf > internal_exons.bed 

# Collapse to unique sorted exons
sort -k1,1 -k2,2n internal_exons.bed | uniq > internal_exons_sorted_uniq.bed

# Include exonic ends of 80 bp and intronic ends of 500 bp. Create different files for 3 prime and 5 prime splice sites.
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($6=="-") {print $1,$3-80,$3+500,$6} else {print $1,$2-500,$2+80,$6} }' internal_exons_sorted_uniq.bed | uniq > 5prSS_exonic80bp_intronic500bp_unique.bed 

awk 'BEGIN {FS="\t"; OFS="\t"} {if ($6=="+") {print $1,$3-80,$3+500,$6} else {print $1,$2-500,$2+80,$6} }' internal_exons_sorted_uniq.bed | uniq > 3prSS_exonic80bp_intronic500bp_unique.bed

# Convert the .bed files to six-field format 
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,"NA","0",$4}' 5prSS_exonic80bp_intronic500bp_unique.bed > 5prSS_exonic80bp_intronic500bp_unique_6col.bed

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,"NA","0",$4}' 3prSS_exonic80bp_intronic500bp_unique.bed > 3prSS_exonic80bp_intronic500bp_unique_6col.bed


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
	
	export ref=5prSS_exonic80bp_intronic500bp_unique_6col
	bwtool extract bed ${ref}.bed ${f}.plus.bw ${f}.5prSS.plus.matrix
	bwtool extract bed ${ref}.bed ${f}.minus.bw ${f}.5prSS.minus.matrix
	
	
	export ref=3prSS_exonic80bp_intronic500bp_unique_6col
	bwtool extract bed ${ref}.bed ${f}.plus.bw ${f}.3prSS.plus.matrix
	bwtool extract bed ${ref}.bed ${f}.minus.bw ${f}.3prSS.minus.matrix
	
	cat ${f}.5prSS.plus.matrix ${f}.5prSS.minus.matrix > ${f}.5prSS.matrix
	cat ${f}.3prSS.plus.matrix ${f}.3prSS.minus.matrix > ${f}.3prSS.matrix

# Calculate the average nucleosomal density for upstram and downstream of the exons.

awk '{$1=$2=$3=$4=$5=$6=$7=""; print $0}' Histone_sorted_ext_clean.5prSS.matrix | tr -s 'NA' '0' | tr -d "-" | awk -F, '{for(i=1;i<=NF;i++)a[i]+=$i} END{for(i=1;i<=NF;i++)printf "%d%s", a[i], (i==NF?"\n":",")}'  > sum5s

awk '{$1=$2=$3=$4=$5=$6=$7=""; print $0}' Histone_sorted_ext_clean.3prSS.matrix | tr -s 'NA' '0' | tr -d "-" | awk -F, '{for(i=1;i<=NF;i++)a[i]+=$i} END{for(i=1;i<=NF;i++)printf "%d%s", a[i], (i==NF?"\n":",")}'  > sum3s

tr -s ',' '\n' < sum5s > sum5s_mod
tr -s ',' '\n' < sum3s > sum3s_mod

awk 'BEGIN {for (i=50; i< 630; i++) print i}' > 3prime_bp_number
awk 'BEGIN {for (i=1; i<=580; i++) print -i}' | sort -n > 5prime_bp_number

paste 3prime_bp_number sum3s_mod > sum3s_mod1
paste 5prime_bp_number sum5s_mod > sum5s_mod1

wc -l Histone_sorted_ext_clean.5prSS.matrix > x1
wc -l Histone_sorted_ext_clean.3prSS.matrix > x2

#Get the avergae by dividing the sum by the total number of regions obtained from the above command.
awk '{print $1"\t"$2/x1}' sum5s_mod1 > sum5s_avg 
awk '{print $1"\t"$2/x2}' sum3s_mod1 > sum3s_avg


# - visualize 3pr SS on the left and 5pr SS on the right in the same chart






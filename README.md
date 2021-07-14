# IR-H3K4me3-signal

Toolkit for the calculation of nucleosomal occupancy on long non-coding RNA genes.

**#Dependencies**
The toolkit employs the following software dependencies

# - download and Install bedtools 
# - download and install LiftOver.
# - download bedGraphToBigWig and bedItemOverlapCount from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/
# - download chromInfo from http://hgdownload.cse.ucsc.edu	/goldenpath/hg38/database/
# - download bwtool 

**# Prerequisite** 
# - Download reference genome data from Gencode
# - Get the gtf file of evidence-based annotation of the human genome (GRCh38), version 34 (Ensembl 100) - long non-coding RNAs
# - Get all exons from UCSC Table Browser UCSC refSeq table: hg38_refSeq_exons.bed
# - Get all introns+80 bp in each end from UCSC Table Browser UCSC refSeq table: hg38_refSeq_introns+80bp.bed
# Notice that first it had to be converted with UCSC liftOver from hg18->hg19 and hg19->hg38
# The input to LiftOver Utility is obatined from http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz 
# and http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# - Human total nucleosome library is obtained from NCBI Short Read Archive (SRA) accession number SRP000105.
# - Histone modified nucleosome libraries are obatined from https://dir.nhlbi.nih.gov/papers/lmi/epigenomes/hgtcell.aspx
![image](https://user-images.githubusercontent.com/73449138/125560353-d66924d3-eb8f-468d-83c2-1a9462805482.png)






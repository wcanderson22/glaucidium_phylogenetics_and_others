#!/bin/bash

# Call snps in samtools

ref="/nas/will.sequences.dec.2024/ncbi_dataset/data/GCA_003259725.1/GCA_003259725.1_athCun1_genomic.fna"
bamdir="/data0/will.temp.2024.1/fastqs/sorted_bam_files/"

ID="Nopo_rerun" # This will be used as a prefix for the output file

echo "making a pileup file for" $ID >> log
#can also add the -R flag joined with a scaffold list to subset and parralel
bcftools mpileup -Ou -f $ref --ignore-RG -a AD,ADF,DP,SP,INFO/AD,INFO/ADF \
"$bamdir"*.bam | bcftools call -mv > "$ID"_raw_variants.vcf
#echo "removing all lines with two comment marks" >> log
#grep -v "##" "$ID"_snps_indels.vcf > "$ID"_snps_indels_short.vcf
echo "filtering low quality snps (<80)" >> log
awk '$1~/^#/ || $6 > 80 {print $0}' > \
"$ID"_snps_indels_filtered.vcf "$ID"_raw_variants.vcf
#echo "add the header and check length of column 4 and 5 to make sure
#they are snp type variants" >> log
#awk '$1~/^#/ || length($4)==1 && length($5)==1 {print $0}'> \
#"$ID"_snps_filtered.vcf "$ID"_snps_indels_filtered.vcf
samtools-snp-pipeline.sh (END)

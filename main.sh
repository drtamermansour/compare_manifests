mkdir compare_manifests

## Get manifest files
module load rclone ## Loading rclone/1.65.1
rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/Miscellaneous documents_standardbred/Equine80select_v2_1_HTS_20143333_B1_UCD.csv" --drive-shared-with-me .
mv Equine80select_v2_1_HTS_20143333_B1_UCD.csv Equine80select.v2.csv
cp $HOME/Horse_parentage_SNPs/backup_original/Equine80select_24_20067593_B1.csv Equine80select.v1.csv
cp $HOME/Equine80select_remapper/results/Equine80select_24_20067593_B1_remapped_equCab3.csv Equine80select.v1.remapped.csv

## run manifest_compare.py
conda activate remap
python manifest_compare.py
rclone -v copy output/ "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/compArr/"  --drive-shared-with-me
rclone -v copy logs/ "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/compArr/"  --drive-shared-with-me

## match_tier
tail -n+2 output/master_comparison_table.csv | cut -d"," -f8 | sort | uniq -c

## multi_match - Many-to-many positional matches
tail -n+2 output/master_comparison_table.csv | cut -d"," -f4,5,10 | sort | uniq -c
tail -n+2 output/master_comparison_table.csv | awk 'BEGIN{FS=OFS=","}{if($10=="True")print $3,$54,$15,$8}' 

## EC3 conflict — v1 native EC3 coords differ from v1.remapped coords
tail -n+2 output/master_comparison_table.csv | cut -d"," -f4,5,16 | sort | uniq -c ## 81
tail -n+2 output/master_comparison_table.csv | awk 'BEGIN{FS=OFS=","}{if($16=="True")print $3,$40":"$41,$54,$15,$8}' 

## flag_unreliable_position — MAPQ_TopGenomicSeq == 0
tail -n+2 output/master_comparison_table.csv | cut -d"," -f4,5,22 | sort | uniq -c ## 1361
tail -n+2 output/master_comparison_table.csv | awk 'BEGIN{FS=OFS=","}{if($22=="True")print $4,$5,$8}' | sort | uniq -c

#######
## Make list of perfect matching markers
tail -n+2 output/master_comparison_table.csv | awk 'BEGIN{FS=OFS=","}{if($8=="strict")print $3}' > strict_match.lst
cut -f3 $HOME/Equine80select_remapper/results/matchingSNPs_binary_consistantMapping.equCab3_map > remapped_QC.lst
cut -f2 $HOME/genDiv/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.bim > filtered.lst

wc -l *.lst
   # 58106 filtered.lst     === This is the old filtered list. We need to rerun genDiv after the final remapping
   # 79669 remapped_QC.lst
   # 79970 strict_match.lst
comm -12 <(sort remapped_QC.lst) <(sort strict_match.lst) | wc -l ## 79258
comm -12 <(sort filtered.lst) <(sort strict_match.lst) | wc -l ## 57739




github_page=https://raw.githubusercontent.com/vappiah/Transcriptomics/main/Tutorials/gene-expression-analysis/Episode-1
for script in assemble_transcripts.sh estimate_abundance.sh map_reads.sh prepare_data.sh sort_mapped_reads.sh visualize_deg.R 
do

wget "$github_page"/$script
done
wget https://raw.githubusercontent.com/vappiah/Transcriptomics/main/environments/environment_1.yaml

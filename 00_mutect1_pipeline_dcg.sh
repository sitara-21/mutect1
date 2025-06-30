#!/bin/bash
# This script processes tumor-normal BAM files to run MuTect1 for detecting somatic mutations.

# Ensure exactly three arguments are provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <step> <start_row> <stop_row>"
  exit 1
fi

###############################################################
# NOTE:
# MuTect1 only supports VCF version 4.1.
# Ensure input VCF files are converted to v4.1 format using:
# ./convert_vcf_to_v4.1.sh in.vcf out.vcf
###############################################################

# Parse input arguments
step=$1          # Pipeline step to execute (0 for setup, 1 for processing rows)
start_row=$2     # Start processing from this row in the CSV file
stop_row=$3      # Stop processing at this row

# File paths and configurations
csv="/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/mutation_calling/mutect1/tumor_normal_duplex.csv"               # Input CSV file
out_dir="/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/mutation_calling/mutect1_output"    # Output directory
email=""               # Email to receive job notifications

# Reference genome and gnomAD files (ensure that the gnomAD VCF is v4.1)
#ref_genome=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/b37/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta   # b37
#gnomad=/n/data1/hms/dbmi/gulhan/lab/software/helper_scripts/mutation_calling_scripts/af-only-gnomad.raw.sites.b37_vers4.1.vcf                    
ref_genome=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/hg38/cgap_matches/Homo_sapiens_assembly38.fa       # hg38
gnomad=/n/data1/hms/dbmi/gulhan/lab/software/helper_scripts/mutation_calling_scripts/af-only-gnomad.hg38_vers4.1.vcf.gz                  

# Path to interval BED files
## b37 ##
#interval_path=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/b37/Homo_sapiens_assembly19/bed               # chromosome interval
#interval_path=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/b37/Homo_sapiens_assembly19/bed/10mil_bp_bed  # 10 million bp interval

## hg38 ##
interval_path=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/hg38/high_confidence_regions/by_chr
# interval_path=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/hg38/high_confidence_regions/10million_bp

# Path to Mutect1 package
mutect1_path=/n/data1/hms/dbmi/gulhan/lab/software/helper_scripts/mutation_calling_scripts/mutect1

# Resource allocation
partitions="short"
time_limit="8:00:00"
mem_alloc="15G"


# Check input file and directory
if [[ ! -f ${csv} ]]; then
  echo "Error: Input CSV file not found at $csv"
  exit 1
fi

if [[ ! -d ${out_dir} ]]; then
  echo "Output directory ${out_dir} does not exist. Creating..."
  mkdir -p ${out_dir} ${out_dir}/tmp
fi


# Steps of the pipeline
case $step in

  ### Step 0: Environment Setup ###
  0)
    echo "module load gcc java/jdk-1.7u80"
  ;;

  ### Step 1: Run Mutect1 ###
  1)
    row=1
    while IFS=, read -r tumor_bam normal_bam sample_id; do
        if [[ $row -lt $start_row ]]; then
            row=$((row + 1))
            continue
        fi

        if [[ $row -gt $stop_row ]]; then
            break
        fi

        if [[ $row -eq 1 && "${sample_id}" == "SampleID" ]]; then
            echo "Skipping header row..."
            row=$((row + 1))
            continue
        fi

        echo "Processing row $row: $sample_id"
        row=$((row + 1))

        for bed_file in ${interval_path}/*bed; do
            segment_id=$(basename ${bed_file} | cut -d'.' -f1)
            output_file="${out_dir}/${sample_id}_${segment_id}.txt"

            if [[ -f "${output_file}" ]]; then
                echo "${output_file} exists. Skipping..."
                continue
            fi

            sbatch -J "mt1_${sample_id}_${segment_id}" \
                   -p ${partitions} \
                   -t ${time_limit} \
                   --mem=${mem_alloc} \
                   --mail-type=FAIL --mail-user=${email} \
                   -o "/n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/mutect1_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/mutect1_%j.err" \
                   --wrap="java -Xmx15g -Djava.io.tmpdir=${out_dir}/tmp \
                   -jar ${mutect1_path}/mutect-1.1.7.jar \
                   -T MuTect \
                   -R ${ref_genome} \
                   -I:tumor ${tumor_bam} \
                   -I:normal ${normal_bam} \
                   -L ${bed_file} \
                   -dcov 10000 \
                   --dbsnp ${gnomad} \
                   --max_alt_alleles_in_normal_count 10000 \
                   --max_alt_alleles_in_normal_qscore_sum 10000 \
                   --fraction_contamination 0 \
                   --tumor_f_pretest 0 \
                   --initial_tumor_lod 0 \
                   --tumor_lod 0 \
                   --out ${output_file}"
        done
    done < $csv
  ;;

  *)
    echo "Invalid step: $step"
    exit 1
  ;;
esac

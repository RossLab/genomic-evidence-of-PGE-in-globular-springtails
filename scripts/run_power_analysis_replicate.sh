#!/bin/bash

# qsub -o logs -e logs -cwd -N power_analysis -V -pe smp64 2 -b yes "scripts/run_power_analysis_replicate.sh 0.001 10 25 15"

# set -e

theta=$1
X_chromosomes=$2 # number of 1Mbp chromosomes that are X-linked.
autosomes=$((20 - $X_chromosomes)) # each genome is 20Mbp long, so 20 - X -> autosomes
maternal_coverage=$3
paternal_coverage=$4

source_dir=$(pwd)
sim=sim_"$JOB_ID"_h"$theta"_X"$X_chromosomes"_m"$maternal_coverage"_p"$paternal_coverage"

working_dir=/scratch/kjaron/"$sim"

reference_template=data/reference/Afus1/genome_short_headers.fa

mkdir -p $working_dir && cd $working_dir
mkdir -p scripts simulation tmp output
cp "$source_dir"/scripts/generate_reference_subset.py scripts

### generating reference for this simulation
# -l cen specify length of each chromosome, but it's 1Mbp by default, so no need to adjust that
python3 scripts/generate_reference_subset.py -g "$source_dir"/"$reference_template" -a "$source_dir"/tables/chr_assignments_Afus1.tsv -c "A" -n $autosomes -o simulation/sim_reference_A.fasta
python3 scripts/generate_reference_subset.py -g "$source_dir"/"$reference_template" -a "$source_dir"/tables/chr_assignments_Afus1.tsv -c "X" -n $X_chromosomes -o simulation/sim_reference_X.fasta

# the merged reference will used for mapping reads
cat simulation/sim_reference_A.fasta simulation/sim_reference_X.fasta > simulation/sim_reference.fasta
echo "Step 1 done: New reference has been generated (simulation/sim_reference.fasta)"

cp "$source_dir"/scripts/create_divergent_haplotypes.py scripts

# here I will need to sprinkle mutations on autosomes and X chromosomes separately, so I will use the corresponding files
python3 scripts/create_divergent_haplotypes.py -g simulation/sim_reference_A.fasta -het $theta -o simulation/sim_genome_A
python3 scripts/create_divergent_haplotypes.py -g simulation/sim_reference_X.fasta -het $theta -o simulation/sim_genome_X --monosomic

echo "Step 2 done: Maternal and Paternal genomes with mutations were simulated (simulation/sim_genome_?_[mp]at.fasta)"

cp "$source_dir"/scripts/simulate_reads.py scripts
# maternal A
python3 scripts/simulate_reads.py -g simulation/sim_genome_A_mat.fasta  -c $maternal_coverage -o simulation/reads_maternal_A
# maternal X
python3 scripts/simulate_reads.py -g simulation/sim_genome_X_mat.fasta -c $maternal_coverage -o simulation/reads_maternal_X
# paternal A
python3 scripts/simulate_reads.py -g simulation/sim_genome_A_pat.fasta  -c $paternal_coverage -o simulation/reads_paternal_A

cat simulation/reads_maternal_A_R1.fq simulation/reads_maternal_X_R1.fq simulation/reads_paternal_A_R1.fq > simulation/sim_reads_R1.fq
cat simulation/reads_maternal_A_R2.fq simulation/reads_maternal_X_R2.fq simulation/reads_paternal_A_R2.fq > simulation/sim_reads_R2.fq

echo "Step 3 done: Reads were generated (simulation/sim_reads_R?.fq)"

cp "$source_dir"/scripts/two_tissue_model.R "$source_dir"/scripts/two_tissue_model_functions.R scripts
ls simulation/sim_reads_R?.fq > simulation/FILES

kmc -k21 -t4 -m16 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram output/kmcdb_k21.hist -cx1000

Rscript scripts/two_tissue_model.R -i output/kmcdb_k21.hist -o output/two_tissue_model

echo "Step 4 done: k-mer histogram generated and the model fit (output/kmcdb_k21.hist and output/two_tissue_model*)"

cp "$source_dir"/scripts/depth2windows.pl "$source_dir"/scripts/mapping_coverage_tab2cov_estimates.R scripts

REFERENCE=simulation/sim_reference
bowtie2-build simulation/sim_reference.fasta $REFERENCE

R1=simulation/sim_reads_R1.fq
R2=simulation/sim_reads_R2.fq

bowtie2 --very-sensitive-local -p 4 -x $REFERENCE \
        -1 $R1 -2 $R2 \
        --rg-id 1 --rg SM:rand --rg PL:ILLUMINA --rg LB:LIB_rand \
        | samtools view -h -q -20 \
        | samtools sort -O bam - > simulation/rand.rg.sorted.bam

samtools depth simulation/rand.rg.sorted.bam | perl scripts/depth2windows.pl 10000 > output/per_window_coverage.tsv

Rscript scripts/mapping_coverage_tab2cov_estimates.R -i output/per_window_coverage.tsv -o output/mapping_coverage

echo "Step 5 done: mapping 1n and 2n coverages estimated"

freebayes -f simulation/sim_reference.fasta -b simulation/rand.rg.sorted.bam --standard-filters --min-coverage 5 -p 2 > output/snp_calls_raw.vcf

echo "Step 6 done: SNPs called"

rm -r simulation
rsync -av --remove-source-files $working_dir $source_dir/data/simulations/

echo "Done. With everything!"

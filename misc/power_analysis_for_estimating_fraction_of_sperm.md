### The actual power analysis

Here I am designing individual components of the pipeline. In the end I will wrap everything in a single script so each job will do all the tasks on local scatch disc,

#### Simulating the reference

Script `scripts/generate_reference_subset.py` can generate a subset of the genome reference with desired total length. They will be randomly sampled scaffolds that are already in the genome. And it was tested to have deterministic behavior when seed is specified. Example use would be

```
python3 scripts/generate_reference_subset.py -g data/reference/Afus1/genome_short_headers.fa -a tables/chr_assignments_Afus1.tsv -c "A" -l 1 -n 19 -o simulation/Afus_sim_reference_A.fasta
python3 scripts/generate_reference_subset.py -g data/reference/Afus1/genome_short_headers.fa -a tables/chr_assignments_Afus1.tsv -c "X" -l 1 -o simulation/Afus_sim_reference_X.fasta

cat simulation/sim_reference_A.fasta simulation/sim_reference_X.fasta > data/sim_reference_complete.fasta
```

which generates 1Mbp of X-linked reference sequences.

#### Simulating variants on top of the refenrece

```
python3 scripts/create_divergent_haplotypes.py -g simulation/sim_reference_A.fasta -het 0.003 -o simulation/sim_genome_A
python3 scripts/create_divergent_haplotypes.py -g simulation/sim_reference_X.fasta -het 0.003 -o simulation/sim_genome_X
```

#### Simulating the reads

```
# maternal A
python3 scripts/simulate_reads.py -g simulation/sim_genome_A_mat.fasta  -c 50 -r 150 -o simulation/maternal_A

# maternal X
python3 scripts/simulate_reads.py -g simulation/sim_genome_X_mat.fasta -c 50 -r 150 -o simulation/maternal_X

# paternal A
python3 scripts/simulate_reads.py -g simulation/sim_genome_A_pat.fasta  -c 50 -r 150 -o simulation/paternal_A

cat simulation/maternal_A_R1.fq simulation/maternal_X_R1.fq simulation/paternal_A_R1.fq > simulation/sim_reads_R1.fq
cat simulation/maternal_A_R2.fq simulation/maternal_X_R2.fq simulation/paternal_A_R2.fq > simulation/sim_reads_R2.fq
```

#### K-mer based estimate of 1n and 2n

```
mkdir -p tmp

# inidividual subpart libraries
ls simulation/maternal_A*.fq > simulation/FILES
kmc -k21 -t2 -m4 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram simulation/kmcdb_k21_maternal.hist -cx1000

ls simulation/maternal_X*.fq > simulation/FILES
kmc -k21 -t2 -m4 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram simulation/kmcdb_k21_maternal_X.hist -cx1000

ls simulation/paternal_A*.fq > simulation/FILES
kmc -k21 -t2 -m4 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram simulation/kmcdb_k21_paternal.hist -cx1000

# joint histogram
ls simulation/maternal_A*.fq simulation/maternal_X*.fq simulation/paternal_A*.fq > simulation/FILES
kmc -k21 -t2 -m4 -ci1 -cs1000 @simulation/FILES simulation/kmcdb tmp
kmc_tools transform simulation/kmcdb histogram simulation/kmcdb_k21.hist -cx1000
```

Alright, now we generated from the simualted data a k-mer histogram!

```{R}
tab_mat_A <- read.table('simulation/kmcdb_k21_maternal.hist')
tab_pat_A <- read.table('simulation/kmcdb_k21_paternal.hist')
tab_mat_X <- read.table('simulation/kmcdb_k21_maternal_X.hist')

tab <- read.table('simulation/kmcdb_k21.hist')

par(mfrow = c(2, 2))
plot(tab$V1[4:100], tab$V2[4:100])
plot(tab_mat_A$V1[4:100], tab_mat_A$V2[4:100])
plot(tab_pat_A$V1[4:100], tab_pat_A$V2[4:100])
plot(tab_mat_X$V1[4:100], tab_mat_X$V2[4:100])
# that all looks alright!

x <- tab_mat_A$V1[4:100]
y <- tab_mat_A$V2[4:100]
nlsLM(y ~ dnbinom(x, size = kmer_cov / bias, mu = kmer_cov) * length,  start = list(kmer_cov = 27, length = 19e6, bias = 1))


x <- tab_pat_A$V1[4:100]
y <- tab_pat_A$V2[4:100]
nlsLM(y ~ dnbinom(x, size = kmer_cov / bias, mu = kmer_cov) * length,  start = list(kmer_cov = 27, length = 19e6, bias = 1))


x <- tab_mat_X$V1[4:100]
y <- tab_mat_X$V2[4:100]
nlsLM(y ~ dnbinom(x, size = kmer_cov / bias, mu = kmer_cov) * length,  start = list(kmer_cov = 27, length = 1e6, bias = 1))
```

but let's see about two tissue model

```{R}
Rscript scripts/two_tissue_model.R -i simulation/kmcdb_k21.hist
```

#### Mapping reads and estimating the coverage

```
REFERENCE=simulation/generated_complete_sim_reference
bowtie2-build simulation_complete_sim_reference.fasta $REFERENCE

R1=simulation/sim_reads_R1.fq
R2=simulation/sim_reads_R2.fq

bowtie2 --very-sensitive-local -p 4 -x $REFERENCE \
        -1 $R1 -2 $R2 \
        --rg-id 1 --rg SM:rand --rg PL:ILLUMINA --rg LB:LIB_rand \
        | samtools view -h -q -20 \
        | samtools sort -O bam - > simulation/rand.rg.sorted.bam
```

Mapped reads to a table of coverages (coverage per 10,000 nt).

```
samtools depth simulation/rand.rg.sorted.bam | perl scripts/depth2windows.pl 10000 > simulation/per_window_coverage.tsv
```

```
Rscript scripts/mapping_coverage_tab2cov_estimates.R -i simulation/per_window_coverage.tsv -o simulation/test
```
### Explored parametric space

```
heterozygosity = 0.01, 0.1, 0.5 (autosomal average heterozygosity in %)
X_chromosomes = 1 2 5 10 (Mbp out of 20 Mbp reference)
fraction_of_sperm = 0, 1, 5, 10, 25, 50 (0 is practically non-PGE system)
sequencing_depth = 10x, 15x, 25x (1n coverage)

female = 25
male = 25, 24.75, 23.75, 22.5, 18.75, 12.5

female = 15
male = 15, 14.85, 14.25, 13.50, 11.25, 7.50

female = 10

```

### Running it on cluster

Let's create a conda enviorment for the analysis called `CSKS`.

```
conda create -n CSKS -c bioconda -c conda-forge msprime kmc python=3 numpy matplotlib wgsim pyfaidx samtools
# I had to add a few package later on:
# conda install -c bioconda wgsim
# conda install -c bioconda pyfaidx
# conda install -c bioconda samtools
# conda install -c bioconda freebayes
conda install -c r r
# and install 'argparse' and 'minpack.lm' packages
# install.packages('argparse')
# install.packages('minpack.lm')
```

now with this conda enviorment we should be able to test one replicate

```
qsub -o logs -e logs -cwd -N power_analysis -V -pe smp64 2 -b yes "scripts/run_power_analysis_replicate.sh 0.001 10 25 15"
```

### The actual power analysis

Here I am designing individual components of the pipeline. In the end I will wrap everything in a single script so each job will do all the tasks on local scatch disc,

#### Simulating the reference

Script `scripts/generate_reference_subset.py` can generate a subset of the genome reference with desired total length. They will be randomly sampled scaffolds that are already in the genome. And it was tested to have deterministic behavior when seed is specified. Example use would be

```
python3 scripts/generate_reference_subset.py -g data/reference/Afus1/genome_short_headers.fa -a tables/chr_assignments_Afus1.tsv -c "A" -l 1 -n 19 -o data/generated/Afus_sim_reference_A.fasta
python3 scripts/generate_reference_subset.py -g data/reference/Afus1/genome_short_headers.fa -a tables/chr_assignments_Afus1.tsv -c "X" -l 1 -o data/generated/Afus_sim_reference_X.fasta

cat data/generated/sim_reference_A.fasta data/generated/sim_reference_X.fasta > data/sim_reference_complete.fasta
```

which generates 1Mbp of X-linked reference sequences.

#### Simulating variants on top of the refenrece

```
python3 scripts/create_divergent_haplotypes.py -g data/generated/sim_reference_A.fasta -het 0.003 -o data/generated/sim_genome_A
python3 scripts/create_divergent_haplotypes.py -g data/generated/sim_reference_X.fasta -het 0.003 -o data/generated/sim_genome_X
```

#### Simulating the reads

```
# maternal A
python3 scripts/simulate_reads.py -g data/generated/sim_genome_A_mat.fasta  -c 50 -r 150 -o data/generated/maternal_A

# maternal X
python3 scripts/simulate_reads.py -g data/generated/Afus_reference_X_ch.fasta -c 50 -r 150 -o data/generated/maternal_X

# paternal A
python3 scripts/simulate_reads.py -g data/generated/Afus_reference_subset_alt.fasta  -c 50 -r 150 -o data/generated/paternal_A

cat data/generated/maternal_A_R1.fq data/generated/maternal_X_R1.fq data/generated/paternal_A_R1.fq > data/generated/sim_reads_R1.fq
cat data/generated/maternal_A_R2.fq data/generated/maternal_X_R2.fq data/generated/paternal_A_R2.fq > data/generated/sim_reads_R2.fq
```

#### K-mer based estimate of 1n and 2n

```
ls data/generated/maternal_A*.fq data/generated/maternal_X*.fq data/generated/paternal_A*.fq > data/generated/FILES
mkdir -p tmp

kmc -k21 -t2 -m4 -ci1 -cs1000 @data/generated/FILES data/generated/kmcdb tmp
kmc_tools transform data/generated/kmcdb histogram data/generated/kmcdb_k21.hist -cx1000
```

Alright, now we generated from the simualted data a k-mer histogram!

```{R}
tab <- read.table('data/generated/kmcdb_k21.hist')
plot(tab$V1[3:100], tab$V2[3:100])
```

but let's see about two tissue model

```{R}
Rscript scripts/two_tissue_model.R -i data/generated/kmcdb_k21.hist
```

#### Mapping reads and estimating the coverage

```
REFERENCE=data/generated/generated_complete_sim_reference
bowtie2-build data/generated_complete_sim_reference.fasta $REFERENCE

R1=data/generated/sim_reads_R1.fq
R2=data/generated/sim_reads_R2.fq

bowtie2 --very-sensitive-local -p 4 -x $REFERENCE \
        -1 $R1 -2 $R2 \
        --rg-id 1 --rg SM:rand --rg PL:ILLUMINA --rg LB:LIB_rand \
        | samtools view -h -q -20 \
        | samtools sort -O bam - > data/generated/rand.rg.sorted.bam
```

reads -> conclusions

### Explored parametric space

```
heterozygosity = 0.01, 0.1, 0.5 (autosomal average heterozygosity in %)
X_chromosomes = 1 2 5 10 (Mbp out of 20 Mbp reference)
fraction_of_sperm = 0, 1, 5, 10, 25, 50 (0 is practically non-PGE system)
sequencing_depth = 10x, 15x, 25x (1n coverage)
```

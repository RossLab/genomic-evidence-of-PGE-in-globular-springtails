### The actual power analysis

Here I am designing individual components of the pipeline. In the end I will wrap everything in a single script so each job will do all the tasks on local scatch disc,

#### Simulating the reference

Script `scripts/generate_reference_subset.py` can generate a subset of the genome reference with desired total length. They will be randomly sampled scaffolds that are already in the genome. And it was tested to have deterministic behavior when seed is specified. Example use would be

```
python3 scripts/generate_reference_subset.py -g data/reference/Afus1/genome_short_headers.fa -a tables/chr_assignments_Afus1.tsv -c "A" -l 1 -o data/generated/Afus_reference_subset.fasta
python3 scripts/generate_reference_subset.py -g data/reference/Afus1/genome_short_headers.fa -a tables/chr_assignments_Afus1.tsv -c "X" -l 1 -o data/generated/Afus_reference_X_ch.fasta
```

which generates 1Mbp of X-linked reference sequences.

#### Simulating variants on top of the refenrece

```
python3 scripts/create_divergent_haplotypes.py -g data/generated/Afus_reference_subset.fasta -het 0.003
```

#### Simulating the reads

```
# maternal A
python3 scripts/simulate_reads.py -g data/generated/Afus_reference_subset_ref.fasta  -c 50 -r 150 -o data/generated/maternal_A

# maternal X
python3 scripts/simulate_reads.py -g data/generated/Afus_reference_X_ch.fasta -c 50 -r 150 -o data/generated/maternal_X

# paternal A
python3 scripts/simulate_reads.py -g data/generated/Afus_reference_subset_alt.fasta  -c 50 -r 150 -o data/generated/paternal_A
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

#### Mapping reads and estimating the coverage


### Explored parametric space

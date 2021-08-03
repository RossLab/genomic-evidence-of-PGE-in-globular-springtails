### Stuart's alternative plan

```
qsub -o logs -e logs -cwd -N bam2sync -V -pe smp64 8 -b yes 'samtools mpileup -a --no-BAQ --fasta-ref data/reference/Afus1/genome.fa --output /scratch/$USER/Afus1.mpileup data/mapped_reads/Afus1.rg.sorted.rmdup.bam && java -jar --input ~/src/popoolation2/mpileup2sync.jar /scratch/$USER/Afus1.mpileup --threads 16 --output /scratch/$USER/Afus1.sync && rsync -av --remove-source-files /scratch/$USER/Afus1.sync data/sync_files/'
```

I gave the sync file to Stuart and he made some magic with it that resulted in a strange figure supporting 0.09 as the number. Not super sure what the number is (if difference in covrage supports, or ratio?). Anyway, to understand better what we are dealing with I decided to explore the sync file (`data/sync_files/Afus1.sync.gz`) or even the raw pileup files stored on `biggar`.

I checked the number of lines corresponsing with each nt in mpileup files and sync files and it matches. I wonder where the difference with fasta index comes from then, must be in the index or mpileup process.

```
zcat data/sync_files/Afus1.sync.gz | scripts/sync2variable_sites.py > data/sync_files/variant_sites.tsv
```


```R
variant_sites <- read.table("data/sync_files/variant_sites.tsv")
colnames(variant_sites) <- c('scf', 'pos', 'major_cov', 'minor_cov')

variant_sites$cov_ratio <- variant_sites$minor_cov / (variant_sites$minor_cov + variant_sites$major_cov)

hist(variant_sites$cov_ratio)
# only the error variants are properly visible, probably better somehow cleared. Ideas:
#  - min minor cov
#  - min total cov
#  - only assigned scf

asn_tab <- read.table('tables/chr_assignments_Afus1.tsv', header = T)
strlen <- max(nchar(c(asn_tab$scf, variant_sites$scf)))

names2tokens <- function(scf_name){
  sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
}

rownames(asn_tab) <- names2tokens(asn_tab$scf)
variant_sites$chr <- asn_tab[names2tokens(variant_sites$scf), 'chr']
variant_sites$len <- asn_tab[names2tokens(variant_sites$scf), 'len']

assigned_vars <- variant_sites[!is.na(variant_sites$chr), ]
assigned_vars$total_cov <- assigned_vars$major_cov + assigned_vars$minor_cov

hist(assigned_vars$total_cov[assigned_vars$total_cov < 140,])
hist(assigned_vars$total_cov[assigned_vars$total_cov < 140 & assigned_vars$chr == 'X'], add = T, col = 'red')
hist(assigned_vars$total_cov[assigned_vars$total_cov < 140 & assigned_vars$chr == 'A'], add = T, col = 'blue')

X_variants <- assigned_vars[assigned_vars$total_cov < 140 & assigned_vars$chr == 'X', ]
A_variants <- assigned_vars[assigned_vars$total_cov < 140 & assigned_vars$chr == 'A', ]

hist(X_variants$cov_ratio[X_variants$minor_cov > 5])
hist(A_variants$cov_ratio[A_variants$minor_cov > 5])

scfs <- unique(A_variants$scf)
scf <- scfs[4]
per_scf_cov_ratios <- sapply(scfs, function(scf){ mean(A_variants[A_variants$minor_cov > 5 & A_variants$scf == scf, 'cov_ratio'])} )

# plot(log10(asn_tab[names2tokens(scfs), 'len']) ~ per_scf_cov_ratios, pch = 1, cex = 2)
hist(per_scf_cov_ratios)
```
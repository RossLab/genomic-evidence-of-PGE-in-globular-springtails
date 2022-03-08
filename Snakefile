# Space for preprocessing

reseq_samples = []
raw_read_files = []
sample_accesions = dict()

with open('tables/Allacma_reseq_samples.tsv') as tab :
	tab.readline()
	for textline in tab :
		line = textline.split()
		ind = line[0]
		reseq_samples.append(ind)
		# print(line[1])
		# line[3] is SRA accession number
		for lib in line[1].split(',') :
			raw_lib_file = 'data/raw_reads/' + ind + '/' + lib + '_1.fastq.gz'
			trimmed_lib_file = 'data/trimmed_reads/' + ind + '/' + lib + '-trimmed-pair1.fastq.gz'
			raw_read_files.append(raw_lib_file)
			sample_accesions[ind] = sample_accesions.get(ind, []) + [trimmed_lib_file]

localrules: help, download_all_reads, map_all

##
## help : print this help
rule help :
	shell :
		"sed -n 's/^##//p' Snakefile"

##
## download_all_reads : download all raw reads from ENA
rule download_all_reads :
	input :
		raw_read_files

##
## map_all : map all reads to corresponding references
rule map_all :
	input :
		expand("data/mapped_reads/{ind}.rg.sorted.bam", ind=reseq_samples), "data/mapped_reads/Afus1.rg.sorted.bam"

rule download_reads :
	threads : 1
	resources : mem=2000000, tmp=30000
	output : "data/raw_reads/{ind}/{accesion}_1.fastq.gz", "data/raw_reads/{ind}/{accesion}_2.fastq.gz"
	shell : "scripts/download_reads.sh {wildcards.ind} {wildcards.accesion}"

rule trim_reads :
	threads : 16
	resources : mem=80000000, tmp=50000
	input : "data/raw_reads/{ind}/{accesion}_1.fastq.gz"
	output : "data/trimmed_reads/{ind}/{accesion}-trimmed-pair1.fastq.gz"
	shell : "scripts/trim_reads.sh data/raw_reads/{wildcards.ind}/{wildcards.accesion}_[1,2].fastq.gz data/trimmed_reads/{wildcards.ind}/{wildcards.accesion}"

rule index_reference :
	threads : 8
	resources : mem=2000000, tmp=30000
	input : "data/reference/{sp}/genome.fa.gz"
	output: "data/reference/{sp}/idx.1.bt2"
	shell : "bowtie2-build --threads 8 <(zcat {input}) data/reference/{wildcards.sp}/idx"

rule map_reads_to_Afus1 :
	threads : 16
	resources : mem=104857600, tmp=50000
	input : "data/reference/Afus1/idx.1.bt2", lambda wildcards: sample_accesions[wildcards.ind]
	output : "data/mapped_reads/{ind}.rg.sorted.bam"
	shell : "scripts/map_reads_to_Afus1.sh data/trimmed_reads/{wildcards.ind}/ {output}"

# rule get_data :
# 	threads : 1
# 	output : "data/{sp}.txt"
# 	shell : "scripts/get_data.sh {wildcards.sp} > {output}"
#
# rule trim_reads :
# 	threads : 1
# 	input : "data/hummingbird.txt", "data/mealybug.txt"
# 	output : "data/final_result.txt"
# 	shell : "wc -l {input} 1> {output}"

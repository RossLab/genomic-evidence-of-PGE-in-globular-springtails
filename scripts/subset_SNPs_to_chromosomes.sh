#!/bin/bash

# IN_VCF=data/SNP_calls/freebayes_all_samples_raw.vcf
# IN_VCF=data/SNP_calls/freebayes_Ocin2_raw.vcf
IN_VCF=$1
# CHR_ASN=tables/chrosmome_asn.tsv
# CHR_ASN=tables/Ocin_X-linage_assigment.tsv
CHR_ASN=$2

VCF_FILT=${IN_VCF%".vcf"}"_filt.vcf"
VCF_FILT_ANS=${IN_VCF%".vcf"}"_filt_asn_only.vcf"
VCF_STATS=${IN_VCF%".vcf"}"_filt.stats"

SCRATCH="/scratch/$USER/snp_filt"

mkdir -p $SCRATCH/$(dirname $IN_VCF)

vcffilter -f "QUAL > 20" -g "DP > 5" $IN_VCF > $SCRATCH/$VCF_FILT

bcftools stats $SCRATCH/$VCF_FILT > $SCRATCH/$VCF_STATS

bgzip $SCRATCH/$VCF_FILT

tail -n+2 $CHR_ASN | awk '{ print $1 "\t1\t" $2 }' > $SCRATCH/assigined_scaffolds.tsv

bcftools index "$SCRATCH/$VCF_FILT".gz
bcftools view -R $SCRATCH/assigined_scaffolds.tsv "$SCRATCH/$VCF_FILT".gz > $SCRATCH/$VCF_FILT_ANS

rsync -av --remove-source-files "$SCRATCH/$VCF_FILT".gz $SCRATCH/$VCF_FILT_ANS $SCRATCH/$VCF_STATS $(dirname $IN_VCF)

#!/usr/bin/env bash

mkdir -p logs

snakemake $@ -p --jobs 10 --latency-wait 120 --cluster "qsub \
    -N {rule} \
    -pe smp64 \
    {threads} \
    -cwd \
    -b y \
    -o \"logs/{rule}.{wildcards}.out\" \
    -e \"logs/{rule}.{wildcards}.err\""

# The previous memory managment was done via:
# -M {resources.mem} \
# On SGE, the memory should be managed via -l h_vmem=size
# However here: http://bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html
# they say it's per processor slot, while here: http://gridscheduler.sourceforge.net/htmlman/htmlman5/queue_conf.html?pathrev=V62u5_TAG
# They say "The per-job maximum memory limit in bytes."
# not sure what am I missing, but for now, I will just not specify memory at all

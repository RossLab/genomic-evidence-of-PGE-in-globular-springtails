#!/bin/bash

samtools depth $1 | scripts/depth2depth_per_contig_median.py > $2

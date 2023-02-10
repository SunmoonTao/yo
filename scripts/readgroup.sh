#!/bin/bash

id=$(basename -s _L001_R1_001.fastq.gz $1)
echo "Read Group @RG\tID:$id\tSM:$id"\tLB:$id"\tPL:ILLUMINA"

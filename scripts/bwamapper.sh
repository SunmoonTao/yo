#!/bin/bash

header=$(zcat $1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

id=$(basename -s _L001_R1_001.fastq.gz $1)
echo "Read Group @RG\tID:$id\tSM:$idLB:$id"\tPL:ILLUMINA"

ml BWA/0.7.16a-foss-2017a

bwa mem \
-M \
-t 8 \
-v 3 \
-R $(echo "Read Group @RG\tID:$id\tSM:$id\tLB:$id"\tPL:ILLUMINA") \
"$path_bwaindex_genome" \
$1 $2 | samblaster -M | samtools fixmate - - | samtools sort -O bam -o "$1.mapped-bwa.bam"
#!/bin/sh
PREFIX=../data
BINARY=jannovar/jannovar_0.26/jannovar-cli-0.26-SNAPSHOT.jar
REFSEQ=jannovar/jannovar_0.26/data/hg19_refseq.ser
REFFASTA=referenceGenome/data/human_g1k_v37.fasta

java -jar $PREFIX/$BINARY \
hgvs-to-vcf \
-d $PREFIX/$REFSEQ \
-r $PREFIX/$REFFASTA \
--server

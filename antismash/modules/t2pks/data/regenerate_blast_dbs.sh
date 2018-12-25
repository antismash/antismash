#!/bin/bash

for fname in *.fasta; do
    NAME=${fname/\.fasta/}
    rm $NAME.p??
    makeblastdb -dbtype prot -in $fname -title $NAME -out $NAME
done

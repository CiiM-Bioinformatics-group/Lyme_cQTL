#!/bin/bash

#Name of the QTL output from MatrixEQTL
input=

#1. Simple subsetting to be able to read the files at least
head -n 1 $input.txt > $input_0.05.txt
head -n 1 $input.txt > $input_5e-8.txt
echo 'headers added'
awk '{if($5 < 0.05)print}' $input.txt >> $input_0.05.txt
echo 'Nominal file created'
head $input_0.05.txt
awk '{if($5 < 5e-8)print}' $input_0.05.txt >> $input_5e-8.txt
echo 'GW file created'
head $input_5e-8.txt

#2. Split the files per cytokine
awk 'NR>1{print >> "$input_pergene/"$2}' $input.txt
echo 'Split per gene created'

#3. Input for FUMA
sort -k1,1n -k5,5g $input_5e-8.txt | awk '!a[$1] {a[$1] = $5} $5 == a[$1]' | sed 's/SNP/CHR:POS:REF:ALT%SNP/g' | tr ':' '\t' | tr '%' '\t' | sed 's/p-value/P/g' > $input_5e-8_fuma.txt
echo 'FUMA input created'


#4. File for mannhattan
sort -k1,1n -k5,5g $input.txt | awk '!a[$1] {a[$1] = $5} $5 == a[$1]' | sed 's/SNP/CHR:POS:REF:ALT%SNP/g' | tr ':' '\t' | tr '%' '\t' | sed 's/p-value/P/g' > $input_manhattan.txt
echo 'Manhattan input created'


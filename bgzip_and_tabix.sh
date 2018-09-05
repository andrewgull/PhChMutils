#!/bin/bash
# ungzipped copies of vcf files will remain unchanged 

for f in *.vcf; do
  #echo $f
  vcf-sort $f | bgzip -c > ${f}.gz ;
  #bgzip $f &&
  tabix -p vcf ${f}.gz
done

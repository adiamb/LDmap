# LDmap
A class method to map LD structure of top snps with 1000genomes and further indicate if top snps from GWAS meta were genotyped or not
# Example command
```
python MATCH_LD_TOP_VARIANTS_CLASS_METHOD.py -meta CHR2_LDSCRIPT_TEST_NOV2018.meta \
-ld CHR2_180342618_231763011.ld \
-refld CLEAN_1KG_CHR2_180342618_231763011.ld \
-chrnum 2 \
-outfile test_2LDscript 
-snplist KLS_GENOTYPES_FILE
```

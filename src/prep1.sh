prj=$(pwd)

cd raw/p1/gno

# standardize: only SNPs, only chromosome 1-23
plink --bfile 000 --snps-only just-acgt --make-bed --chr 1-23 --out S00

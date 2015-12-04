## prepare phase 2 twin genotype data
cd dat/p2.0

## rearrange Illumina GenomicStudio output
## 1) extract genotype variants, drop R and Theta values, and the the header
## 2) transform X, Y, XY and MT chromosome to numeric coding
## 3) sort the variants by chromosome position
cut g.txt -f1,2,3,$(seq -s, 4 3 $(head -n 1 g.txt | wc -w)) | tail -n +2| g.txt.1
sed g.txt.1 -e 's/\tX\t/\t23\t/' -e 's/\tY\t/\t24\t/' -e 's/\tXY\t/\t25\t/' -e 's/\tMT\t/\t26\t/' > g.txt.2
sort -k 2g -k 3g g.txt.2 > g.txt.3

## get variant map out of sorted genotype data, rewrite in plink
## map format
awk '{print($2,"\t",$1,"\t","0","\t",$3)}' g.txt.3 > g.plink.map

## get value matrix from sorted genotype data, rewrite in a format
## close to plink tped
cut g.txt.3 -f4- | sed -e 's/\([AB]\)\([AB]\)/\1 \2/g' -e 's/NC/0 0/g' > g.plink.tped.mtx

## get subject ID, in the original order
head -n 1 g.txt | cut --output-delimiter $'\n' -f$(seq -s, 4 3 $(head -n 1 g.txt | wc -w)) | tr -d '.GType' > g.txt.sbj

## patch up the plink tped file, then creat plink binary files
paste g.plink.map g.plink.tped.mtx >g.plink.tped
plink --tped g.plink.tped --tfam g.plink.fam --out g.plink.0

## haploid heterogeneity check
plink --bfile g.plink.0 --set-hh-missing --make-bed --out g.plink.1

## prepare phase 1 twin genotype data
cd dat/p1.0

## recode BED to dosage values, set set heterozygous haploids to NA
plink --bfile g --recode A-transpose --set-hh-missing --out g

## remove genetic distance, reference allele and alternative allel columns
## move the feature name to the first column
## use '3' to represent missing genotype (NA)
cut -f3,5,6 --complement g.traw | sed -e 's/^\([^\t]*\)\t\([^\t]*\)/\2\t\1/' -e '1,1 s/[0-9]*_//g' -e '2,$ s/\bNA\b/3/g' > g.dsg.0

## split genome feature map and dosage values matrix
cut -f1-3 g.dsg.0 > g.dsg.map && cut -f1-3 --complement g.dsg.0 > g.dsg.mtx


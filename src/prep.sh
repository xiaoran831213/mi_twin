## prepare phase 2 twin genotype data
cd dat/p2.0

## extract genotype columns, drop R and Theta value from GenomeStudio
cut g.txt -f1,2,3,$(seq -s, 4 3 $(head -n 1 g.txt | wc -w)) > g.nox

## transform AA, BB, Ab coding to dosage values;
## use '3' to represent NC (no calling)
sed g.nox -e '1,1 s/\.GType//g' -e '2,$ s/AA/0/g' -e 's/BB/2/g' -e 's/AB/1/g' -e 's/NC/3/g' > g.dsg.0

## transform X, Y, XY and MT chromosome to numeric coding
sed g.dsg.0 -e 's/\tX\t/\t23\t/' -e 's/\tY\t/\t24\t/' -e 's/\tXY\t/\t25\t/' -e 's/\tMT\t/\t26\t/' > g.dsg.1

## sort the genome features by chromosome and basepair loacation
sort -k 2g -k 3g g.dsg.1 > g.dsg.2

## split genome feature map and dosage matrix
cut -f1,2,3 g.dsg.2 > g.dsg.map && cut -f4- g.dsg.2 > g.dsg.3

## remove subject IDs in the head of dosage matrix
tail -n +2 g.dsg.3 > g.dsg.mtx

## transform texture dosage format to plink tped format
awk 'NR>1 {print($2,"\t",$1,"\t","0","\t",$3)}' g.dsg.map > g.plink.map
sed g.dsg.mtx -e 's/0/A A/g' -e 's/1/A B/g' -e 's/2/B B/g' -e 's/3/0 0/g' > g.plink.tped.mtx
paste g.plink.map g.plink.tped.mtx >g.plink.tped

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


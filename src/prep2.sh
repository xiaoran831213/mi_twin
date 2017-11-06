## prepare phase 2 twin genotype data
prj=$(pwd)
cd raw/p2/seq/Burt_2014_Analysis;
# extract basic data from all Ilumina batches into a working directory
mkdir -p wrk
for f in *Analysis/fdt.txt; do i=$(sed <<<${f} 's/^.*_\([0-9]\)_.*$/\1/'); dos2unix -n $f wrk/$i.txt; done
for f in *Analysis/snp.txt; do i=$(sed <<<${f} 's/^.*_\([0-9]\)_.*$/\1/'); dos2unix -n $f wrk/$i.snp; done

# separate subject id, genomic maps, and basic genotype
for f in wrk/*.txt; do
    echo $f
    # subject ID
    head -n 1 ${f%.*}.txt | cut -f2- --output-delimiter=$'\n' | sed 's/[.].*$//' > ${f%.*}.iid

    # SNP maps
    tail -n+2 ${f%.*}.snp > ${f%.*}.snp.headless; mv ${f%.*}.snp.headless ${f%.*}.snp

    # preliminary genotypes #1
    tail -n+2 ${f%.*}.txt | cut -f2- > ${f%.*}.gx1
done

# sanity check, the SNP files should be identical!
for f in wrk/*.snp; do md5sum $f; done

# combine all bathes
mkdir -p cmb
cat   wrk/*.iid > cmb/gmx.ii0
mv    wrk/1.snp   cmb/gmx.mp0
paste wrk/*.gx1 > cmb/gmx.tx0

cd cmb
# count A/B alleles and NC
tr -d "$'\t' TCGID-" < gmx.tx0 | while read l; do echo ${#l}; done > A.cnt
tr -d "$'\t'A CGID-" < gmx.tx0 | while read l; do echo ${#l}; done > T.cnt
tr -d "$'\t'AT GID-" < gmx.tx0 | while read l; do echo ${#l}; done > C.cnt
tr -d "$'\t'ATC ID-" < gmx.tx0 | while read l; do echo ${#l}; done > G.cnt
tr -d "$'\t'ATCG D-" < gmx.tx0 | while read l; do echo ${#l}; done > I.cnt
tr -d "$'\t'ATCGI -" < gmx.tx0 | while read l; do echo ${#l}; done > D.cnt
tr -d "$'\t'ATCGID " < gmx.tx0 | while read l; do echo ${#l}; done > N.cnt
paste A.cnt T.cnt C.cnt G.cnt I.cnt D.cnt N.cnt | awk -v OFS=$'\t' 'M=$1+$2+$3+$4+$5+$6+$7 {print $1,$2,$3,$4,$5,$6,$7,M}' > cnt.txt

# tranform combined data into plink format
# forge *.fam file
# IID, FID, SEX from phenotype file
tail -n+2 "$prj/raw/09_28_2017 MI_twin Data.csv" | cut -f1,2,5 -d, --output-delimiter=$'\t' | sort -k1 >/tmp/fm0
# IID from genotype
awk <gmx.ii0 -v OFS=$'\t' '{print $1,NR}' | sort -k1 >/tmp/ii0
# forge now:
join /tmp/fm0 /tmp/ii0 -j 1 -o 1.2,1.1,1.3,2.2 -t $'\t' | sort -k4n | awk '{print $1,$2,0,0,$3,1}' > gmx.fm0
awk <gmx.mp0 -v OFS=$'\t' '{if($2=="X") c=23; else if($2=="Y") c=24; else if($2=="XY") c=25; else if($2=="MT") c=26; else c=$2}; {print c,$1,0,$3}' >gmx.mp1
tr <gmx.tx0 '-' '0' | sed -r 's/([ATCGID0])([ATCGID0])/\1 \2/g' >gmx.tx1

paste gmx.mp1 gmx.tx1 > gmx.tped
cp gmx.fm0 gmx.tfam
plink --tfile gmx --maf 0.00001 --make-bed --out gmx


# move the plink binary data to the new directory, clean up the workplace
d=$prj/raw/p2/gno; mkdir -p $d
mv gmx.bed $d/000.bed
mv gmx.bim $d/000.bim
mv gmx.fam $d/000.fam
# cd ..; rm -rf wrk cmb

# ---------------------------------------------------------------------- #
# complete phase 2 monozygote twins' genotypes by duplicating the typed twin.                                                                                                             
# back to the top of project directory                                                                                                                                                    
d=raw/p2/gno
# 1) find out monozygote twins sorted by IID, it is ok to keep phase one twins in the list.
awk -v FS=, <raw/phe.csv '$7==1 {print $2,$1,$4}' | sort -k2b,2 >/tmp/phe.mzt

# 2) filter out typed twin in phase 2, sorted by IID                                                                                                                                      
plink --bfile $d/000 --keep /tmp/phe.mzt --make-bed --out gno.mzt
# letter "G" serves as a mark of having been "Genotyped"                                                                                                                                  
awk <gno.mzt.fam '{print $2,"G"}' | sort -k2,2 >/tmp/gno.mzt

# 3) compile statistics of twins and theire genotypes                                                                                                                                     
# join on IID -> 1.FID [2.IID] 3.TID (0 or 1) 4.GNO (G=genotyped, U=untyped)                                                                                                              
join /tmp/phe.mzt /tmp/gno.mzt -1 2 -2 1 -a 1 -o 1.1,1.2,1.3,2.2 -e U > /tmp/mzt.tx0

# join on FID -> [1.FID] 2.IID_0 3.TID_0 (always 0) 4.GNO_0 5.IID_1 6.TID_1 (always 1) 7.GNO_1                                                                                            
join <(sort /tmp/mzt.tx0 -k1b,1 | awk '$3==0') <(sort /tmp/mzt.tx0 -k1b,1 | awk '$3==1') > /tmp/mzt.tx1

# 4) find out twin pairs that has only one individual genotyped and another untyped,                                                                                                          
# then find out the genotyped one who can be duplicated,                                                                                                                       
# also find out the untyped one who should receive the duplication.
awk </tmp/mzt.tx1 '$4=="U" && $7=="G" {print $1,$5}'  >/tmp/mzg.tx0
awk </tmp/mzt.tx1 '$4=="G" && $7=="U" {print $1,$2}' >>/tmp/mzg.tx0
sort /tmp/mzg.tx0 -k1b,1 >/tmp/mzg.tx1 # sorted by FID                                                                                                                                    

awk </tmp/mzt.tx1 '$4=="U" && $7=="G" {print $1,$2}'  >/tmp/mzu.tx0
awk </tmp/mzt.tx1 '$4=="G" && $7=="U" {print $1,$5}' >>/tmp/mzu.tx0
sort /tmp/mzu.tx0 -k1b,1 >/tmp/mzu.tx1 # sorted by FID                                                                                                                                    

# 5) duplication                                                                                                                                                                          
# get the genotyped twin                                                                                                                                                                  
plink --bfile $d/000 --keep /tmp/mzg.tx1 --make-bed --out mzg
# change individual ID to another twin, while maintaining the order of genotype samples                                                                                                   
awk <mzg.fam '{print $1,$2,$3,$4,$5,$6,NR}' | sort -k1b,1 | join - /tmp/mzu.tx1 | sort -k7n,7 | awk '{print $1,$8,$3,$4,$5,$6}' > mzu.fam
# sanity check                                                                                                                                                                            
diff <(cut mzg.fam -f1 -d' ') <(cut mzu.fam -f1 -d' ') # FID, should print nothing                                                                                                        
diff <(cut mzg.fam -f1 -d' ') <(cut mzu.fam -f2 -d' ') # IID, should print a lot                                                                                                          
diff <(cut mzg.fam -f3 -d' ') <(cut mzu.fam -f3 -d' ') # PID, should print nothing                                                                                                        
diff <(cut mzg.fam -f4 -d' ') <(cut mzu.fam -f4 -d' ') # MID, should print nothing                                                                                                        
diff <(cut mzg.fam -f5 -d' ') <(cut mzu.fam -f5 -d' ') # SEX, should print nothing                                                                                                        
diff <(cut mzg.fam -f6 -d' ') <(cut mzu.fam -f6 -d' ') # PHE, should print nothing                                                                                                        
# rename and merge                                                                                                                                                                        
# mzu.fam already created                                                                                                                                                                 
cp mzg.bed mzu.bed
cp mzg.bim mzu.bim
plink --bfile raw/p2/gno/000 --bmerge mzu --make-bed --out raw/p2/gno/001


# check allele consistency between phase 1 and phase 2 genotypes.
join <(sort 000.bim -k2b,2) <(sort p01.bim -k2b,2) -j 2 -t $'\t' | cut -f5,6,10,11 | grep "^[ATCG]" | awk '$1==$3 {print}'


# standardize: only SNPs, only chromosome 1-23
plink --bfile $d/001 --snps-only just-acgt --make-bed --chr 1-23 --out $d/S00

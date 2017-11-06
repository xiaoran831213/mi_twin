## merge standardized phase 1 and phase 2 twin genotypes, only SNPs on chromosome 1~23 are included.
prj=$(pwd)
d=raw/p3/gno; mkdir -p $d

# 1) diagnosis
# identify mis-matched SNP location
join <(sort raw/p1/gno/S00.bim -k2b,2) <(sort raw/p2/gno/S00.bim -k2b,2) -t $'\t' -j 2 | awk '$2!=$7 || $4!=$9' # 6
join <(sort raw/p1/gno/S00.bim -k2b,2) <(sort raw/03_04_2015_ABC.gwa -k1b,1) -t $'\t' -1 2 -2 1 | awk '$4!=$8' # 0
join <(sort raw/p2/gno/S00.bim -k2b,2) <(sort raw/03_04_2015_ABC.gwa -k1b,1) -t $'\t' -1 2 -2 1 | awk '$4!=$8' # 1
# phase 2 has one more poistion mismatch than phase 1: rs208961

# identify strain mis-matched allele
join <(sort raw/p1/gno/S00.bim -k2b,2) <(sort raw/p2/gno/S00.bim -k2b,2) -t $'\t' -j 2 | awk '!($5==$10 && $6==$11 || $5==$11 && $6==$10)' # 1
join <(sort raw/p1/gno/S00.bim -k2b,2) <(sort raw/03_04_2015_ABC.gwa -k1b,1) -t $'\t' -1 2 -2 1 | awk '!($5==$9 && $6==$10 || $5==$10 && $6==$9)' # 12
join <(sort raw/p2/gno/S00.bim -k2b,2) <(sort raw/03_04_2015_ABC.gwa -k1b,1) -t $'\t' -1 2 -2 1 | awk '!($5==$9 && $6==$10 || $5==$10 && $6==$9)' # 13
# phase two have one more strand mismatch than phase 1: rs208961
# ** decision **
# for positions, align phase 2 to phase 1 which is slightly closer to the panel of Broad ABC GWAS.
# for strands, drop the only 1 mismatched SNP


# 2) merge phase 1 and 2
# 2.1) select shared SNPs, excluding SNPs with mismatched strand
join <(sort raw/p1/gno/S00.bim -k2b,2) <(sort raw/p2/gno/S00.bim -k2b,2) -t $'\t' -j 2 | awk '$5==$10 && $6==$11 || $5==$11 && $6==$10 {print $1}' >/tmp/p30.snp
plink --bfile raw/p1/gno/S00 --extract /tmp/p30.snp --make-bed --out p1s
plink --bfile raw/p2/gno/S00 --extract /tmp/p30.snp --make-bed --out p2s
# sanity check for strand mismatch, should be empty
join <(sort p1s.bim -k2b,2) <(sort p2s.bim -k2b,2) -t $'\t' -j 2 | awk '!($5==$10 && $6==$11 || $5==$11 && $6==$10)' # 0 (good!)
join <(sort p1s.bim -k2b,2) <(sort p2s.bim -k2b,2) -t $'\t' -j 2 | awk '!($5==$10 && $6==$11)' # 6948 flipped minor/major allele
# 2.2) expand binary (*.{bed,bim,fam}) format to texture format (*.{ped, map}), then concatinate
plink --bfile p1s --recode --out p1s
plink --bfile p2s --recode --out p2s
cat p1s.ped p2s.ped >p30.ped	# concatinate
cp p1s.fam >p30.fam		# use phase 1 SNP positions.
plink --file p30 --make-bed --out $d/S00
# rm p[12]s.{ped,map,log,hh,bed,bim,fam} p30.*

# 3) some quality control for the merged data
cd $d
# maximum subject missing rate : 0.1
# maximum SNP missing rate:      0.1
# switching the above QC makes huge difference that only a handful of SNPs would remain due to a handful (8) of poorly typed individuals,
# removing these individual before checking the SNP-wise missing rate is more preferred.
# setting minimum MAF to 1e-5 removes degenerated SNPs.
plink --bfile $d/S00 --geno .1 --mind .1 --maf .00001 --make-bed --out $d/S01

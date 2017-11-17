# phase 3.0 analysis
# merged MI twin study from phase 1 and phase 2 genotype
mkdir -p dat/p3_0

# Calculation of Polygenic Risk Scores
# the basic GWAS output remain the same: Broad_ABC
bas=raw/03_04_2015_ABC.gwa
# target genotye that merged from phase 1 & 2
tgt=raw/p3/gno/S01
# use dummy phenotype to trick Prisce1.2 into calculating a complete set of risk scores
phe=dat/p3_0/dmy.txt; awk -v OFS=$'\t' < $tgt.fam '{print $2,1}' > $phe

# user.covariate.file $cvr
R --file=src/PRSice_v1.21.R -q --args \
    plink ~/bin/plink \
    base "$bas" \
    target "$tgt" \
    pheno.file "$phe" \
    binary.target F \
    covary F \
    slower 0 \
    sinc 0.01 \
    supper 0.5 \
    clump.r2 0.25 \
    clump.kb 500 \
    report.individual.scores T \
    ggfig F \
    no.regression T
mv PRSice_SCORES_AT_ALL_THRESHOLDS.txt dat/p3_0/prs.txt


# PCA
module load Eigensoft/6.0.1
d=dat/p3_0/PCA; mkdir -p $d
awk -v OFS=$'\t' <$tgt.fam '{print $1,$2,$3,$4,$5,1}' > $d/pop.pedind
cp $tgt.bed $d/pop.bed
cp $tgt.bim $d/pop.pedsnp

cd $d
smartpca.perl -i pop.bed -a pop.pedsnp -b pop.pedind -o pop.pcs -p pop_pc2 -e pop.egv -l pop.log -k 20 -t 0
cd ../../..
tail -n+2 $d/pop.pcs.evec | sed 's/^  *//; s/:/\t/; s/  */\t/g' | sed 's/\tControl$//' > dat/p3_0/pcs_t00.txt

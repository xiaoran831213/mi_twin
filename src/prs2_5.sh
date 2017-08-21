# merge MI twin study phase 1 and 2 genotype
mkdir -p dat/p2_5

# unify major/minor allele into 1/2 first, because phase 2 coding were still in A/B format
plink --bfile dat/p2_1/g0 --maf 0.00001 --recode 12 --out dat/p2_5/p1
plink --bfile dat/p2_0/g2 --maf 0.00001 --recode 12 --out dat/p2_5/p2

# merge
# phase 1 major/minor coding overwrite phase 2
plink --file dat/p2_5/p2 --merge dat/p2_5/p1 --merge-mode 5 --make-bed --out dat/p2_5/p21
# phase 2 major/minor coding overwrite phase 1
plink --file dat/p2_5/p1 --merge dat/p2_5/p2 --merge-mode 5 --make-bed --out dat/p2_5/p12

# phase 2.5 analysis, merged phase 1 and 2
# the basic GWAS output remain the same
bas=dat/p1_1/Base_BroadABC_final_forPRS_03042015.txt

# target genotye, merged 1 & 2
tgt=dat/p2_5/p2

# extract CBCL_Mom_cond_r for r2 calculation, whites only, ignore family structure issue for now
dir=dat/p2_5/CBCL_Mom_cond_r; mkdir -p $dir; phe=$dir/y.txt
# awk -v OFS=$'\t' <dat/p1_0/phe.txt '$6==6  {print $2,$8}'  > $phe # phase 1
awk -v OFS=$'\t' <dat/p2_0/phe.txt '$8==6 {print $1,$10}' > $phe # phase 2

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
    no.regression F
# mv PRSice_SCORES_AT_ALL_THRESHOLDS.txt dat/p1_2/psc_gender1_allsnp.txt
# mv PRSice_RAW_RESULTS_DATA.txt       out/p1_2/psc_gender1_raw_result.txt
# mv PRSice_SCORES_AT_BEST-FIT-PRS.txt out/p1_2/psc_gender1_best_score.txt

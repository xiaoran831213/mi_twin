# replicate meta-analysis result in the MI Twin cohort for Alex Burt Lab, v2.1
# v2.x deal with the 2nd phase of MI twin study

# polygenic risk score based on previous GWAS
# phase 2 study takes the same base GWAS score used by phase 1 analysis,
# allele recoded to 1/2 format, where A1=1 for minor allele, A2=2 for reference allele.
bas=dat/p1_1/Base_BroadABC_final_forPRS_03042015.txt
awk <$bas 'NR==1 {print}; NR>1 {print $1,$2,$3,1,2,$6,$7}' > dat/p2_1/bas.txt
bas=dat/p2_1/bas.txt

# target genotye, which is designed for phase 2 MI twin study, allele re-coded as 1/2
tgt=dat/p2_0/g3

# extract CBCL_Mom_cond_r for r2 calculation, whites only, ignore family structure issue for now
dir=dat/p2_1/CBCL_Mom_cond_r; mkdir -p $dir; phe=$dir/y.txt
awk -v OFS=$'\t' <dat/p2_0/phe.txt '$8==6 {print $1,$10}' >$phe

# covariat are useless in this analysis, since Prsice1.2 only use them to infer ancestory info!

# analyze now
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

R --file=src/PRSice_v1.21.R -q --args \
    plink ~/bin/plink \
    base "$bas2" \
    target "dat/p1_0/g" \
    pheno.file "$phe2" \
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
# mv PRSice_SCORES_AT_ALL_THRESHOLDS.txt dat/p1_2/psc_gender2_allsnp.txt


# polygenic risk scores based on top SNPs, gender specific
# the reports were manually prepared.
R --file=src/PRSice_v1.21.R -q --args \
    plink ~/bin/plink \
    base "dat/p1_2/TopSnp_gender1" \
    target "dat/p1_0/g" \
    pheno.file "$ped1" \
    binary.target F \
    covary F \
    slower 0 \
    sinc 0.01 \
    supper 0.5 \
    clump.r2 0.25 clump.kb 500 report.individual.scores T no.regression T
mv PRSice_SCORES_AT_ALL_THRESHOLDS.txt psc_topSnp_gender1.txt

R --file=src/PRSice_v1.21.R -q --args \
    plink ~/bin/plink \
    base "dat/p1_2/TopSnp_gender2" \
    target "dat/p1_0/g" \
    pheno.file "$ped2" \
    binary.target F \
    covary F \
    slower 0 \
    sinc 0.01 \
    supper 0.5 \
    clump.r2 0.25 clump.kb 500 report.individual.scores T no.regression T
mv PRSice_SCORES_AT_ALL_THRESHOLDS.txt psc_topSnp_gender2.txt

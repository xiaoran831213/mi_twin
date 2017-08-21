# replicate meta-analysis result in the MI Twin cohort for Alex Burt Lab, v2.1
# v2.x deal with the 2nd phase of MI twin study

# extract dummy phenotype file, for score calculation only.
ped=/tmp/ped.txt; awk < dat/p2_0/phe.txt -v 'OFS=\t' 'NR<2 {print "IID","PHE"}; NR>1 && $6==6 {print $2,"0"}' > $ped

# polygenic risk score based on previous GWAS
# the raw report must be sorted, and given the header required by PRSics, here we take the same base score used by
# phase 1 analysis
# echo "SNP CHR BP A1 A2 BETA P" > $bas1; tail -n+2 dat/p1_2/BroadABC_Males2017 | sort -k2,2n -k3,3n | tr 'atcg' 'ATCG' >> $bas1
bas=dat/p1_1/Base_BroadABC_final_forPRS_03042015.txt
# awk <$bas 'NR==1 {print}; NR>1 {print $1,$2,$3,2,1,$6,$7}' > dat/p1_1/bas.txt
# bas=dat/p1_1/bas.txt

# target genotye for phase 1 MI twin, allele MAY BE recoded into 1/2;
tgt=dat/p1_0/g0
#tgt=dat/p1_0/g1

# extract CBCL_Mom_cond_r for r2 calculation, whites only, ignore family structure issue for now
dir=dat/p1_1/CBCL_Mom_cond_r; mkdir -p $dir; phe=$dir/y.txt
awk -v OFS=$'\t' <dat/p1_0/phe.txt '$6==6 {print $2,$8}' >$phe # must not contain header liner.

# extract genter, sex of European Americans, keep file header
# !!! Turns out covariat is not usefult in GLM analysis, Prcise only use them for ancestry info !!
# cvr=dat/p2_1/cvr.txt
# awk -v OFS=$'\t' <dat/p2_0/phe.txt '$8==6 || NR==1 {print $1,$6,$7}' > $cvr

# extract TRF_rule_r for r2 calculation, whites only, ignore family structure issue for now
# dir=dat/p2_1/TRF_rule_r; mkdir -p $dir; phe=$dir/y.txt
# awk <dat/p2_0/phe.txt '$8==6 {OFS="\t"; print $1,$15}' >$phe

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

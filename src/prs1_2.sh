# replicate meta-analysis result in the MI Twin cohort

# extract dummy phenotype file, for score calculation only.
ped1=/tmp/ped1.txt; awk < dat/p1_0/phe.txt -v 'OFS=\t' 'NR<2 {print "IID","PHE"}; NR>1 && $4==1 && $6==6 {print $2,"0"}' > $ped1
ped2=/tmp/ped2.txt; awk < dat/p1_0/phe.txt -v 'OFS=\t' 'NR<2 {print "IID","PHE"}; NR>1 && $4==2 && $6==6 {print $2,"0"}' > $ped2

# polygenic risk score based on all SNPs, gender specific
# the raw report must be sorted, and given the header required by PRSics
bas1=dat/p1_2/BroadABC_gender1
bas2=dat/p1_2/BroadABC_gender2
# echo "SNP CHR BP A1 A2 BETA P" > $bas1; tail -n+2 dat/p1_2/BroadABC_Males2017 | sort -k2,2n -k3,3n | tr 'atcg' 'ATCG' >> $bas1
# echo "SNP CHR BP A1 A2 BETA P" > $bas2; tail -n+2 dat/p1_2/BroadABC_Females   | sort -k2,2n -k3,3n | tr 'atcg' 'ATCG' >> $bas2

# extract CBCL_Mom_cond_r for r2 calculation, ignore family structure issue.
phe1=dat/p1_2/CBCL_Mom_cond_r_gender1.txt
phe2=dat/p1_2/CBCL_Mom_cond_r_gender2.txt
awk <dat/p1_0/phe.txt '$4==1 && $6==6 {OFS="\t"; print $2,$8}' >$phe1
awk <dat/p1_0/phe.txt '$4==2 && $6==6 {OFS="\t"; print $2,$8}' >$phe2

phe1=dat/p1_2/TRF_rule_r_gender1.txt
phe2=dat/p1_2/TRF_rule_r_gender2.txt
awk <dat/p1_0/phe.txt '$4==1 && $6==6 {OFS="\t"; print $2,$15}' >$phe1
awk <dat/p1_0/phe.txt '$4==2 && $6==6 {OFS="\t"; print $2,$15}' >$phe2

R --file=src/PRSice_v1.21.R -q --args \
    plink ~/bin/plink \
    base "$bas1" \
    target "dat/p1_0/g" \
    pheno.file "$phe1" \
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

## extract genter, sex of European Americans, keep file header
awk '$6==6 || NR==1 {print $2"\t"$4"\t"$5"\t"$8}' <dat/sav/rv1.phe > cov

## extract CBCL_Mom_cond_r phenotype of European Amerians, remove header
awk '$6==6 {print $2"\t"$6"\t"$8}' <dat/sav/rv1.phe > phe

## run PRSice
R --file=src/PRSice_v1.21.R -q --args \
    plink ~/bin/plink \
    base bck/Base_BroadABC_final_forPRS_03042015.txt \
    target dat/sav/rv1 \
    pheno.file phe.eur \
    binary.target F \
    covary F \
    slower 0 \
    sinc 0.01 \
    supper 0.5 \
    clump.r2 0.25 \
    clump.kb 500 \
    report.individual.scores T \
    no.regression T

## run PRSice Summary-Summary statistics
# sz=$(wc -l < phe)
# sz=$(($sz-1))
# R --file=src/PRSice_v1.21.R -q --args \
#     plink ~/bin/plink \
#     sumsum T \
#     base bck/Base_BroadABC_final_forPRS_03042015.txt \
#     target t5 \
#     pheno.file phe \
#     size.targ $sz \
#     binary.target F \
#     covary T \
#     user.covariate.file cov \
#     covariates gender,age \
#     slower 0 \
#     sinc 0.01 \
#     supper 0.5 \
#     clump.ref dat/rv2 \
#     clump.r2 0.25 \
#     clump.kb 500 \
#     best.thresh.on.bar T \
#     report.individual.scores T \

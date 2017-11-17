## responses
rsp={SCICA_cond_r,SCICA_aggr_r,SCICA_aggs_r,SCICA_rules_r,CBCL_Mom_cond_r,CBCL_Mom_rule_r,CBCL_Mom_aggr_r,ASR_mom_aggr_r,ASR_mom_rule_r,ASR_mom_antis_r,CBCL_Dad_cond_r,CBCL_Dad_rule_r,CBCL_Dad_aggr_r,ASR_dad_aggr_r,ASR_dad_rule_r,ASR_dad_antis_r,TRF_cond_r,TRF_rule_r,TRF_aggr_r,CP_all,CP_family,RB_adults,AGG_adults}
print "%s\n" {TRF_rule_r,CBCL_Mom_cond_r}


## self replication trial
y=TRF_rule_r
hpcwp "Rscript -e 'source(\"prs3_0.R\"); main(50, \"$y\", once.psi, npc=10, sav=\"{n:04d}.rds\")'" -d tmp/bs1 -i 20 -t2 -m4 -b4 -q1 --cp src/prs3_0.R --ln dat,raw --tag bs1

y=TRF_rule_r
hpcwp "Rscript -e 'source(\"prs3_0.R\"); main(50, \"$y\", once.bin, npc=10, sav=\"{n:04d}.rds\")'" -d tmp/bs2 -i 20 -t2 -m4 -b4 -q1 --cp src/prs3_0.R --ln dat,raw --tag bs2

y=CBCL_Mom_rule_r
hpcwp "Rscript -e 'source(\"prs3_0.R\"); main(50, \"$y\", once.psi, npc=10, sav=\"{n:04d}.rds\")'" -d tmp/bs3 -i 20 -t2 -m4 -b4 -q1 --cp src/prs3_0.R --ln dat,raw --tag bs3

y=CBCL_Mom_cond_r
hpcwp "Rscript -e 'source(\"prs3_0.R\"); main(50, \"$y\", once.psi, npc=10, sav=\"{n:04d}.rds\")'" -d tmp/bs4 -i 20 -t2 -m4 -b4 -q1 --cp src/prs3_0.R --ln dat,raw --tag bs4

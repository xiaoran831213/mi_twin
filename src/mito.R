source("src/helper.R")
source("src/analysis.R")

## read subject phenotypes
if(!exists('phe', inherits = F))
{
    library(foreign)
    phe = read.spss('dat/rv2.phe.sav', to.data.frame = T)
    phe = within(phe,{
                     FID <- as.integer(famidun);
                     IID <- as.integer(id);
                     id <- NULL;
                     famidun <- NULL;})
    hdr = names(phe)
    hdr = c(tail(hdr, 2), head(hdr, length(hdr) - 2))
    phe = phe[, hdr]
}
                                        # read genotype and meta data
if(!exists('gmx', inherits = F))
{
    bin<-paste('dat/gno.bin')
    if(file.exists(bin))
    {
        load(file = bin);
    }
    else
    {
        ## read genotype individual list
        idv<-paste('dat/rv1.fam')
        idv<-read.table(
            file = idv,
            col.names = c("FID", "IID", "PID", "MID", "SEX", "PHE"),
            as.is = T)
        
        ## read genotype map and allel frequency information
        frq<-read.table(
            file = 'dat/rv1.frq', header = T, as.is = T)
        
        ## read physical map
        map<-paste('dat/rv1.bim')
        map<-read.table(
            file = 'dat/rv1.bim',
            header = F,
            col.names = c("CHR", "SNP", "GDS", "POS", "A1", "A2"))
        map$MAF<-frq$MAF;
        
        ## read genotype matrix
        gmx<-sprintf("tail -n +2 %s | cut -f7- -d ' '", 'dat/rv1.raw')
        tmp<-pipe(gmx, "r");
        gmx<-scan(file = tmp, what = integer(0));
        close(tmp);
        rm(tmp);
        dim(gmx)<-c(nrow(map), nrow(idv));
        
        ## save R binary
        save(gmx, map, idv, file = bin);
        rm(bin);
    }
    gc();
}

## take out only mitocondria genome
mito.idx = which(map$CHR == 26)
map <- map[mito.idx, ]
gmx <- gmx[mito.idx, ]

## define genotype checking function
## gvr ---- genomic variant
## map ---- map information
## nes ---- number of effective sample
gck<-function(gvr, map, nes)
{
    err = T;
    msg = NA;
    if(map$CHR < 26L) # only mitochondria
    {
        msg = 'CHR<26';
    }
    else if(map$MAF < 0.05) # fail MAF test
    {
        msg = 'MAF<.05';
    }
    else if(nes<100L) # fail ESF test (number of effective sample)
    {
        msg = 'NES<100';
    }
    else
    {
        err = F;
    }
    list(err=err, msg=msg);
}

## prepare file surfix and covariant
cov<-c('gender', 'age', 'mage')

if(!exists('tst', inherits = F))
{
    tst <- list(
        "SCICA_cond_r",
        "CBCL_Mom_cond_r",
        "CBCL_Mom_rule_r",
        "CBCL_Mom_aggr_r",
        "CBCL_Dad_cond_r",
        "CBCL_Dad_rule_r",
        "CBCL_Dad_aggr_r",
        "TRF_cond_r",
        "TRF_rule_r",
        "TRF_aggr_r",
        "CP_all",
        "CP_family",
        "RB_adults",
        "AGG_adults")
}

## start analysis
while(length(tst)>0)
{
    rsp<-tst[[1]];
    print(paste('tst for phe: ', rsp));
    out<-run_std(gmx, map, idv, phe, rsp, cov, gck);
    whr<-sprintf('dat/MT.mage.%s.out', rsp);
    write.table(out, whr, row.names = F, quote = F, sep = '\t');
    
    rpt<-data.frame(
        SNP=out$SNP, CHR=out$CHR, POS=out$POS, EFF_ALL=out$A1, NONEFF_ALL=out$A2,
        BETA   = round(out$EST, 4L),
        SE     = round(out$SE,  4L),
        P      = format(out$P, digits = 3L, scientific = T, trim=T),
        AF_COD = round(out$MAF, 4L), 
        HWE=NA, IMP=0L, INFO=1L, INFO_TYPE=NA, N_EFF=out$NES,
        row.names = NULL, stringsAsFactors = F);
    whr<-sprintf('dat/MT.mage.%s.rpt', rsp);
    write.table(rpt, whr, row.names = F, quote = F, sep = '\t');
    
    rm(whr, rsp, out, rpt);
    tst[[1]]<-NULL;
}
gc();

if(exists('rut', inherits = F)) rm(rut);
if(exists('sex', inherits = F)) rm(sex);
if(exists('cov', inherits = F)) rm(cov);
if(exists('sfx', inherits = F)) rm(sfx);

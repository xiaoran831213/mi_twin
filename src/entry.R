source("src/helper.R")
source("src/analysis.R")

## data file root name.
rut<-"dat/rv1";
##sex<-2L;

## read subject phenotypes
if(!exists('phe', inherits = F))
{
    phe<-paste(rut, 'phe', sep='.');
    phe<-read.table(file = phe, header = T, as.is = T)
    #phe<-phe[phe$ethnicity==6L, ]; # only take sample of european origin
    if(!is.null(sex))
    {
        phe<-phe[phe$gender==sex, ]; # only task sample of one gender
    }
}

## read genotype and meta data
if(!exists('gmx', inherits = F))
{
    bin<-paste(rut, 'gno', 'bin', sep='.');
    if(file.exists(bin))
    {
        load(file = bin);
    }
    else
    {
        ## read genotype individual list
        idv<-paste(rut, 'fam', sep='.');
        idv<-read.table(file = idv, col.names = c("FID", "IID", "PID", "MID", "SEX", "PHE") , as.is = T)
        
        ## read genotype map
        ## read allel frequency information
        frq<-paste(rut, 'frq', sep='.');
        frq<-read.table(file = frq, header = T, as.is = T)
        
        ## read physical map
        cls<-c("integer", "character", "NULL", "integer", "character", "character");
        tag<-c("CHR", "SNP", "GDS", "POS", "A1", "A2");
        map<-paste(rut, 'bim', sep='.');
        map<-read.table(file = map, header = F, colClasses =  cls, col.names = tag, )
        map$MAF<-frq$MAF;
        rm(cls, tag, frq);
        
        ## read genotype matrix
        gmx<-sprintf("tail -n +2 %s.%s | cut -f7- -d ' '", rut, 'raw')
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

## define genotype checking function
## gvr ---- genomic variant
## map ---- map information
## nes ---- number of effective sample
gck<-function(gvr, map, nes)
{
    err = T;
    msg = NA;
    if(sex == 2L & map$CHR == 24L)
    {
        msg = 'NO_XY';
    }
    else if(map$CHR > 25L) ## ignore mitochondria
    {
        msg = 'CHR>25';
    }
    else if(map$MAF < 0.05) ## fail MAF test
    {
        msg = 'MAF<.05';
    }
    else if(nes<100L) ## fail ESF test (number of effective sample)
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
if(is.null(sex))
{
    cov<-c('gender', 'age');
    sfx<-'sxc';	
}else
{
    cov<-c('age');
    sfx<-paste('sx', sex, sep='');
}

if(!exists('tst', inherits = F))
{
    tst<-as.list(colnames(phe)[0:-6]);
}

## start analysis
while(length(tst)>0)
{
    rsp<-tst[[1]];
    print(paste('tst for phe: ', rsp));
    out<-run_std(gmx, map, idv, phe, rsp, cov, gck);
    whr<-sprintf('%s.%s.%s.%s', rut, rsp, sfx, 'out');
    write.table(out, whr, row.names = F, quote = F, sep = '\t');
    
    rpt<-data.frame(
        SNP=out$SNP, CHR=out$CHR, POS=out$POS, EFF_ALL=out$A1, NONEFF_ALL=out$A2,
        BETA   = round(out$EST, 4L),
        SE     = round(out$SE,  4L),
        P      = format(out$P, digits = 3L, scientific = T, trim=T),
        AF_COD = round(out$MAF, 4L), 
        HWE=NA, IMP=0L, INFO=1L, INFO_TYPE=NA, N_EFF=out$NES,
        row.names = NULL, stringsAsFactors = F);
    whr<-sprintf('%s.%s.%s.%s', rut, rsp, sfx, 'rpt');
    write.table(rpt, whr, row.names = F, quote = F, sep = '\t');
    
    rm(whr, rsp, out, rpt);
    tst[[1]]<-NULL;
}
gc();

if(exists('rut', inherits = F)) rm(rut);
if(exists('sex', inherits = F)) rm(sex);
if(exists('cov', inherits = F)) rm(cov);
if(exists('sfx', inherits = F)) rm(sfx);

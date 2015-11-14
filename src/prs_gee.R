library(gee)

## read phenotype, and polygenic risk score
phe <- read.table(file='dat/rv1.phe', header=T, as.is=T);
prs <- read.table(file='PRSice_SCORES_AT_ALL_THRESHOLDS.txt', header=T, as.is=T);

## align data
iid <- sort(Reduce(f=intersect, x=list(phe$IID, prs$IID)))
phe=phe[match(table=phe$IID, x=iid, nomatch=FALSE), ]
prs=prs[match(table=prs$IID, x=iid, nomatch=FALSE), ]

##
phe[phe==-9]=NA

## pick data for GEE analysis
dat <- data.frame(
    FID=phe$FID, IID=phe$IID,
    SEX=phe$gender, AGE=phe$age,
    PHE=phe$CBCL_Mom_cond_r)

out <- matrix(double(), nrow=ncol(prs)-1L, ncol=4L);
for(i in 2L:ncol(prs))
{
    dat$PRS <- prs[,i];
    h <- colnames(prs)[i];
    h <- sub('^[^0-9]*', '', h);
    h <- as.double(h);
    
    ## GEE model fit
    m <- gee(PHE~PRS+SEX+AGE, id=FID, data=dat, family='poisson', silent = T)
    s <- summary(m);
    z <- s$coefficients['PRS', 'Robust z'];
    s <- s$coefficients['PRS', 'Robust S.E.'];
    p <- 2*(1-pnorm(abs(z)));
    out[i-1L,] <- c(h, z, s, p);
}

out <- data.frame(PT=out[,1L], Z=out[,2L], SE=out[,3L], P=out[,4L]);
write.table(out, "PRS_GEE", quote = F, sep='\t', row.names=F);
write.table(out[order(out$P),][1L,], "PRS_GEE_TOP", quote = F, sep='\t', row.names=F);

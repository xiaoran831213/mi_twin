source("src/helper.R")
require(gee);

run_std<-function(gmx, map, idv, phe, rsp, cov, chk=NULL)
{
    if(is.null(chk))
    {
        chk<-function(g, m, i, ...) # dummy checking
        {
            list(err=F, msg=NA);
        }
    }
    phe<-phe[, c("FID","IID", rsp, cov)];
    phe<-deNA(phe);

    ## prevent gnome matrix from degenerating to a row vector.
    if(!is.matrix(gmx))
        gmx <- matrix(gmx);
    
    tmp <- align_idv(gmx, idv, phe, odr = T);
    gdx <- tmp$gdx;
    pdx <- tmp$pdx;
    rm(tmp);
    
    ## integrity check
    stopifnot(identical(idv[gdx, 'FID'], phe[pdx, 'FID']));
    stopifnot(identical(idv[gdx, 'IID'], phe[pdx, 'IID']));
    
    ## number of effective sample
    nes<-sapply(1L : nrow(gmx), FUN = function(i) sum(!is.na(gmx[i, gdx])));
    
    ## phenotype
    rsp<-phe[pdx, rsp]   # response variable
    cov<-phe[pdx, cov]   # covariants
    rm(phe);
    
    ## genotype
    fml<-idv[gdx, 'FID']; # family cluster
    ngv<-nrow(gmx);
    
    ## construct statistical modle: P = Beta*G + Gama*C + miu + ...
    mdl<-paste('cov[,',1L:ncol(cov),']', collapse = '+');
    mdl<-paste('rsp~g', mdl, sep = '+');
    mdl<-as.formula(mdl);

    ## prepare output
    out<-matrix(data = NA, nrow = ngv, ncol = 4L, dimnames = list(NULL, c('EST', 'SE', 'Z', 'P')));
    err<-rep.int(NA, ngv);
    
    ## iterate genome variants
    for(i in 1L:ngv)
    {
        if(bitwAnd(i, 0x7FFF) == 1L)
        {
            cat(i, ngv, '\n');
        }
        g<-gmx[i,gdx];
        check<-chk(g, map[i,,drop=T], nes[i]);
        if(check$err)
        {
            err[i]<-check$msg;
            next;
        }

        ## call GEE.
        sink('/dev/null');
        suppressWarnings(suppressMessages(
            r <-try(
                gee(mdl, fml, family=gaussian, corstr="exchangeable",
                        maxiter=100, na.action=na.omit), silent = T)));
        sink(NULL);
        if (inherits(r, "try-error"))
        {
            err[i]<-err_msg(r);
            next;
        }
        
        ## summary may also generate error.
        suppressWarnings(suppressMessages(r <- try(summary(r))));
        if (inherits(r, "try-error"))
        {
            err[i] <- err_msg(r);
            next;
        }
        
        ## get p value.
        r <- unname(r$coefficients[2L, ]);
        p <- try(2*pnorm(-abs(r[5L]))); #5.Robust z
        if (inherits(r, "try-error"))
        {
            err[i] <- err_msg(p);
            next;
        }
        
                                        # 1.Estimate, 4.Robust S.E., 5.Robust z
        out[i, ] <- c(r[1L] , r[4L], r[5L], p);
    }
    out<-data.frame(map, out, NES=nes, ERR=err, stringsAsFactors = F, row.names = NULL);
    out;
}

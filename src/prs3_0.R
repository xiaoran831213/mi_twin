library(gee)
library(ggplot2)
library(magrittr)
library(dplyr)

options(stringsAsFactors=FALSE)
get.data <- function(npc=10)
{
    phe <- read.csv('raw/09_28_2017 MI_twin Data.csv')
    ## identify factors
    phe <- lapply(phe, function(x)
    {
        tb <- table(x)
        nm <- names(tb)
        nm <- as.numeric(nm)
        if(length(nm) < 2)
            x <- NULL
        else if(min(nm) == 0 && max(nm) == 1)
        {
            x <- as.factor(x)
            levels(x) <- (1:length(nm)) - 1
        }
        else if(is.integer(nm))
            x <- as.integer(x)
        x
    })
    factors <- c('gender', 'zygosity', 'ethnicity', 'phase')
    phe[factors] <- lapply(phe[factors], as.factor)
    phe <- phe[!sapply(phe, is.null)]
    phe <- do.call(data.frame, c(phe, check.names=FALSE))
    rownames(phe) <- phe$id

    ## polygenic risk score based on ABC, 1991 subjects
    rsc <- read.table('dat/p3_0/prs.txt', T)
    rownames(rsc) <- rsc$IID
    rsc$IID <- NULL
    ths <- as.numeric(sub('^pT_', '', names(rsc)))
    names(rsc) <- sub('0[.]', '', sprintf("T%.2f", ths))

    ## principle components of genotype
    pcs <- read.table('dat/p3_0/pcs_t00.txt', F)
    rownames(pcs) <- pcs[, 2]
    pcs <- pcs[, -1:-2][, 1:20]
    names(pcs) <- sprintf('P%02d', 1:ncol(pcs))

    ## merge
    iid <- Reduce(intersect, list(rownames(phe), rownames(rsc), rownames(pcs)))
    phe <- cbind(phe[iid, ], rsc[iid, ], pcs[iid, ])
    ## phe <- subset(phe, ethnicity == 6)  # EA only

    ## return
    env <- c("neigh.pov_2", "fam.pov_2", "ecv_2", "disadvan_2")

    list(
        phe=phe,
        env=env,
        cvr=c("gender", "age", head(names(pcs), npc)),
        rsc=names(rsc))
}

nmq <- function(y)
{
    qnorm((rank(y)-0.5)/length(y))
}
nmq.log <- function(y)
{
    y <- log(y)
    qnorm((rank(y)-0.5)/length(y))
}
hist.y <- function(out.dir='~', fun=nmq.log)
{
    dat <- get.data()
    rsp <- with(dat, phe[, rsp])
    graphics.off()
    for(nm in names(rsp))
    {
        y <- rsp[, nm]
        if(!is.null(out.dir))
        {
            fo <- file.path(out.dir, paste0(nm, '.png'))
            png(fo, 7, 7, units='in', res=300)
        }
        else
            dev.new()
        tag <- deparse(substitute(fun))
        hist(fun(y), main=tag)
        if(!is.null(out.dir))
            dev.off()
    }
}

## a single analysis
once.gau <- function(y, r, e, c, dt1, dt2=dt1)
{
    ## clean up
    dt1 <- dt1[c(y, r, e, c, 'famidun')] %>% na.omit
    dt2 <- dt2[c(y, r, e, c, 'famidun')] %>% na.omit
    dt1[, y] <- log(1 + dt1[, y])       # log scale
    dt2[, y] <- log(1 + dt2[, y])       # log scale

    ## model formula
    rxe <- paste(paste(r, '*', e), collapse=' + ')
    mdl.full <- paste(y, '~', rxe, '+', paste(c, collapse=' + '))
    mdl.full <- as.formula(mdl.full)
    mdl.base <- paste(y, '~', r, '+', paste(c, collapse=' + '))
    mdl.base <- as.formula(mdl.base)

    ## do GEE:
    gee.fam <- 'gaussian'               # poisson(log)
    mft.full <- gee(mdl.full, family=gee.fam, data=dt1, id=dt1$famidun)
    mft.base <- gee(mdl.base, family=gee.fam, data=dt1, id=dt1$famidun)
    coef.full <- coef(mft.full)
    coef.base <- coef(mft.base)
    dof <- length(coef.full) - length(coef.base)
    
    ## deviance on training data and testing data
    mm1.full <- model.matrix(mdl.full, dt1)
    mm1.base <- model.matrix(mdl.base, dt1)
    ## dv1.full <- (mm1.full %*% coef.full) %>% exp %>% {-2 * dpois(dt1[, y], ., log=TRUE)}
    ## dv1.base <- (mm1.base %*% coef.base) %>% exp %>% {-2 * dpois(dt1[, y], ., log=TRUE)}
    dv1.full <- (mm1.full %*% coef.full - dt1[, y])^2
    dv1.base <- (mm1.base %*% coef.base - dt1[, y])^2
    dv1 <- dv1.base - dv1.full
    dv1 <- sum(dv1[!is.na(dv1) & !is.infinite(dv1)])
    
    mm2.full <- model.matrix(mdl.full, dt2)
    mm2.base <- model.matrix(mdl.base, dt2)
    ## dv2.full <- (mm2.full %*% coef.full) %>% exp %>% {-2 * dpois(dt2[, y], ., log=TRUE)}
    ## dv2.base <- (mm2.base %*% coef.base) %>% exp %>% {-2 * dpois(dt2[, y], ., log=TRUE)}
    dv2.full <- (mm2.full %*% coef.full - dt2[, y])^2
    dv2.base <- (mm2.base %*% coef.base - dt2[, y])^2
    dv2 <- dv2.base - dv2.full
    dv2 <- sum(dv2[!is.na(dv2) & !is.infinite(dv2)])
    
    r <- data.frame(
        Y=y, T=as.numeric(sub('^T', '.', r)), E=e,
        N=nobs(mft.base), dv1=dv1, dv2=dv2, dof=dof)
    r
}

## a single analysis for binomials
once.bin <- function(y, r, e, c, dt1, dt2=dt1)
{
    ## clean up
    dt1 <- dt1[c(y, r, e, c, 'famidun')] %>% na.omit
    dt2 <- dt2[c(y, r, e, c, 'famidun')] %>% na.omit
    dt1[, y] <- dt1[, y] > 0            # dicotomize
    dt2[, y] <- dt2[, y] > 0            # dicotomize
    
    ## model formula
    rxe <- paste(paste(r, '*', e), collapse=' + ')
    mdl.full <- paste(y, '~', rxe, '+', paste(c, collapse=' + '))
    mdl.full <- as.formula(mdl.full)
    mdl.base <- paste(y, '~', e, '+', paste(c, collapse=' + '))
    mdl.base <- as.formula(mdl.base)
    
    ## do GEE:
    gee.fam <- binomial()
    lnk.inv <- gee.fam$linkinv
    mft.full <- gee(mdl.full, family=gee.fam, data=dt1, id=dt1$famidun) %>% suppressMessages
    mft.base <- gee(mdl.base, family=gee.fam, data=dt1, id=dt1$famidun) %>% suppressMessages
    coef.full <- coef(mft.full)
    coef.base <- coef(mft.base)
    dof <- length(coef.full) - length(coef.base)
    
    ## deviance on training data and testing data
    mm1.full <- model.matrix(mdl.full, dt1)
    mm1.base <- model.matrix(mdl.base, dt1)
    xb1.full <- mm1.full %*% coef.full
    xb1.base <- mm1.base %*% coef.base
    dv1.full <- xb1.full %>% lnk.inv %>% {-2 * dbinom(dt1[, y], 1, ., log=TRUE)}
    dv1.base <- xb1.base %>% lnk.inv %>% {-2 * dbinom(dt1[, y], 1, ., log=TRUE)}
    dv1 <- dv1.base - dv1.full
    dv1 <- sum(dv1)        # sum(dv1[!is.na(dv1) & !is.infinite(dv1)])
    
    mm2.full <- model.matrix(mdl.full, dt2)
    mm2.base <- model.matrix(mdl.base, dt2)
    xb2.full <- mm2.full %*% coef.full
    xb2.base <- mm2.base %*% coef.base
    dv2.full <- xb2.full %>% lnk.inv %>% {-2 * dbinom(dt2[, y], 1, ., log=TRUE)}
    dv2.base <- xb2.base %>% lnk.inv %>% {-2 * dbinom(dt2[, y], 1, ., log=TRUE)}
    dv2 <- dv2.base - dv2.full
    dv2 <- sum(dv2)        # sum(dv2[!is.na(dv2) & !is.infinite(dv2)])
    
    r <- data.frame(
        Y=y, T=as.numeric(sub('^T', '.', r)), E=e,
        N=nobs(mft.base), dv1=dv1, dv2=dv2, dof=dof,
        stringsAsFactors=FALSE)
    r
}

## a single analysis for binomials
once.psi <- function(y, r, e, c, dt1, dt2=dt1)
{
    ## clean up
    dt1 <- dt1[c(y, r, e, c, 'famidun')] %>% na.omit
    dt2 <- dt2[c(y, r, e, c, 'famidun')] %>% na.omit
    
    ## model formula
    rxe <- paste(paste(r, '*', e), collapse=' + ')
    mdl.full <- paste(y, '~', rxe, '+', paste(c, collapse=' + '))
    mdl.full <- as.formula(mdl.full)
    mdl.base <- paste(y, '~', e, '+', paste(c, collapse=' + '))
    mdl.base <- as.formula(mdl.base)
    
    ## do GEE:
    gee.fam <- poisson()
    lnk.inv <- gee.fam$linkinv
    mft.full <- gee(mdl.full, family=gee.fam, data=dt1, id=dt1$famidun)
    mft.base <- gee(mdl.base, family=gee.fam, data=dt1, id=dt1$famidun)
    coef.full <- coef(mft.full)
    coef.base <- coef(mft.base)
    dof <- length(coef.full) - length(coef.base)
    
    ## deviance on training data and testing data
    mm1.full <- model.matrix(mdl.full, dt1)
    mm1.base <- model.matrix(mdl.base, dt1)
    xb1.full <- mm1.full %*% coef.full
    xb1.base <- mm1.base %*% coef.base
    dv1.full <- xb1.full %>% lnk.inv %>% {-2 * dpois(dt1[, y], ., log=TRUE)}
    dv1.base <- xb1.base %>% lnk.inv %>% {-2 * dpois(dt1[, y], ., log=TRUE)}
    dv1 <- dv1.base - dv1.full
    dv1 <- sum(dv1)        # sum(dv1[!is.na(dv1) & !is.infinite(dv1)])
    
    mm2.full <- model.matrix(mdl.full, dt2)
    mm2.base <- model.matrix(mdl.base, dt2)
    xb2.full <- mm2.full %*% coef.full
    xb2.base <- mm2.base %*% coef.base
    dv2.full <- xb2.full %>% lnk.inv %>% {-2 * dpois(dt2[, y], ., log=TRUE)}
    dv2.base <- xb2.base %>% lnk.inv %>% {-2 * dpois(dt2[, y], ., log=TRUE)}
    dv2 <- dv2.base - dv2.full
    dv2 <- sum(dv2)        # sum(dv2[!is.na(dv2) & !is.infinite(dv2)])
    
    r <- data.frame(
        Y=y, T=as.numeric(sub('^T', '.', r)), E=e,
        N=nobs(mft.base), dv1=dv1, dv2=dv2, dof=dof,
        stringsAsFactors=FALSE)
    r
}

main <- function(rep=1, rsp="TRF_rule_r", fun=once.psi, npc=5, sav=NULL)
{
    ## read data
    . <- get.data(npc=npc)
    phe <- .$phe
    ## phe <- subset(.$phe, phase == 2)
    cvr <- .$cvr
    cfg <- with(., expand.grid(r=rsc, e=env, stringsAsFactors=FALSE))
    
    ret <- replicate(rep,
    {
        ## divide families
        fid <- unique(phe$famidun)
        fdx <- sample(fid, length(fid)/2)
        dt1 <- subset(phe,  famidun %in% fdx)
        dt2 <- subset(phe, !famidun %in% fdx)
        
        rpt <- mapply(function(r, e)
        {
            r <- try(fun(rsp, r, e, cvr, dt1, dt2))
            if(inherits(r, 'try-error'))
                NULL
            else r
        }, cfg$r, cfg$e, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        do.call(rbind, rpt)
    }, simplify=FALSE)

    ## combine, save, and return
    ret <- do.call(rbind, ret) %>%
        mutate(pvl.dvp=pchisq(dv1, dof, low=FALSE), pvl.evl=pchisq(dv2, dof, low=FALSE)) %>%
        rename(chq.dvp=dv1, chq.evl=dv2)
    if(!is.null(sav))
        saveRDS(ret, sav)
    invisible(ret)
}

plot.agg <- function(agg, out=NULL, ...)
{
    dot <- list(...)
    ## 1 SD upper and lower bound
    agg <- agg %>% mutate(lw.pvl.dvp = pmax(0, mu.pvl.dvp - sd.pvl.dvp),
                          up.pvl.dvp = pmin(1, mu.pvl.dvp + sd.pvl.dvp),
                          lw.pvl.evl = pmax(0, mu.pvl.evl - sd.pvl.evl),
                          up.pvl.evl = pmin(1, mu.pvl.evl + sd.pvl.evl))
    itr <- round(mean(agg$it))

    plt <- ggplot(agg, aes(x=T))
    plt <- plt + geom_line(aes(y=mu.pvl.dvp, color='dvp'))
    plt <- plt + geom_line(aes(y=mu.pvl.evl, color='evl'))
    if('1sd' %in% names(dot) && dot[['1sd']])
    {
        plt <- plt + geom_ribbon(aes(ymin=lw.pvl.dvp, ymax=up.pvl.dvp, fill='dvp'), alpha=0.2)
        plt <- plt + geom_ribbon(aes(ymin=lw.pvl.evl, ymax=up.pvl.evl, fill='evl'), alpha=0.2)
    }
    plt <- plt + facet_grid(Y~E)
    plt <- plt + ggtitle('Interaction of Enviroment and Polygenic Risk', paste('repetition:', itr))
    plt <- plt + xlab('PRS-thresholds') + ylab('p-values')

    ## no title for legends
    plt <- plt + theme(legend.title=element_blank())
    plt <- plt + scale_fill_discrete(
                     name="Experimental\nCondition",
                     labels=c("developing", "evaluation"))

    if(!is.null(out))
        ggsave(out, plt, ...)
    plt
}

aggr.rpt <- function(rpt)
{
    if(is.character(rpt))
        rpt <- readRDS(rpt)

    ## mean and SD of p-values
    agg <- rpt %>% group_by_at(vars(Y, T, E)) %>%
        summarise(mu.pvl.dvp = mean(pvl.dvp), sd.pvl.dvp = sd(pvl.dvp),
                  mu.pvl.evl = mean(pvl.evl), sd.pvl.evl = sd(pvl.evl),
                  it=n()) %>% ungroup
    agg
}

read.rpt <- function(fs)
{
    rpt <- lapply(dir(fs, 'rds$', full.names=TRUE), readRDS)
    rpt <- do.call(rbind, rpt) %>% as.tbl
    rpt
}

get.all <- function(fs)
{
    r <- read.rpt(fs)
    a <- aggr.rpt(r)
    p <- plot.agg(a)
    list(rpt=r, agg=a, plt=p)
}

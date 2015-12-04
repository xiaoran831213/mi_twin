## calculate missing rate
.msr <- function(x) sum(is.na(x))/length(x)

prep1.0 <- function()
{
    ## surfix appender
    sfx <- function(x) sprintf('dat/p1.0/g.dsg.%s', x)

    ## read map and genomic matrix
    map <- read.table(
        sfx('map'), header = T, sep = '\t',
        col.names = c('vid', 'chr', 'bp1'))
    vid <- map$vid
    gmx <- scan(sfx('mtx'), integer(), skip = 1, na.strings = '3')

    ## read subject id: the first line of geno matrix file
    sid <- scan(sfx('mtx'), integer(), nlines=1)
    dim(gmx) <- c(length(sid), nrow(map))
    gmx <- t(gmx)
    dimnames(gmx) <- list(vid = vid, sid = sid)

    ## calculate MAF, flip genotype if necessary
    maf <- rowMeans(gmx, na.rm = T)/2L
    m50 <- which(maf > .5)
    gmx[m50,] <- 2L - gmx[m50,]
    maf[m50] <- 1 - maf[m50]

    ## calculate missing rates
    ms.gvr <- apply(gmx, 1, .msr)
    ms.sbj <- apply(gmx, 2, .msr)

    ## compile and return
    obj <- list(map=map, gmx=gmx, maf=maf, ms.gvr=ms.gvr, ms.sbj=ms.sbj)
    saveRDS(obj, sfx('rds'))

    obj <- within(obj, {gmx <- NULL; sid <- sid; vid <- vid})
    saveRDS(obj, sfx('stt.rds'))
}


prep2.0 <- function()
{
    ## surfix appender
    sfx <- function(x) sprintf('dat/p2.0/g.dsg.%s', x)

    ## read map and genomic matrix
    map <- read.table(
        sfx('map'), header = T, sep = '\t',
        col.names = c('vid', 'chr', 'bp1'))
    vid <- map$vid
    gmx <- scan(sfx('mtx'), integer(), skip = 1, na.strings = '3')

    ## read subject id: the first line of geno matrix file
    sid <- scan(sfx('mtx'), integer(), nlines=1)
    dim(gmx) <- c(length(sid), nrow(map))
    gmx <- t(gmx)
    dimnames(gmx) <- list(vid = vid, sid = sid)

    ## calculate MAF, flip genotype if necessary
    maf <- rowMeans(gmx, na.rm = T)/2L
    m50 <- which(maf > .5)
    gmx[m50,] <- 2L - gmx[m50,]
    maf[m50] <- 1 - maf[m50]

    ## calculate missing rates
    ms.gvr <- apply(gmx, 1, .msr)
    ms.sbj <- apply(gmx, 2, .msr)

    ## compile and return
    obj <- list(map=map, gmx=gmx, maf=maf, ms.gvr=ms.gvr, ms.sbj=ms.sbj)
    saveRDS(obj, sfx('rds'))

    obj <- within(obj, {gmx <- NULL; sid <- sid; vid <- vid})
    saveRDS(obj, sfx('stt.rds'))
}

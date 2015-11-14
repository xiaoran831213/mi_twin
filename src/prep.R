prepT2 <- function()
{
    ## surfix appender
    sfx <- function(x) sprintf('dat/p2.0/g.dsg.%s', x)

    ## read map and genomic matrix
    map <- read.table(sfx('map'), header = T, sep = '\t')
    gmx <- scan(sfx(mtx), integer(), skip = 1, na.strings = '3')

    ## read subject id: the first line of geno matrix file
    sid <- scan(sfx(mtx), integer(), nlines=1)
    dim(gmx) <- c(length(sid), nrow(map))
    gmx <- t(gmx)
    dimnames(gmx) <- list(vid = map$Name, sid = sid)

    ## caculate MAF, flip genotype if necessary
    maf <- rowMeans(gmx, na.rm = T)/2L

    ## compile and return
    obj <- list(map = map, gmx = gmx)
    saveRDS(obj, 'dat/p2.0/g.dsg.rds')
}

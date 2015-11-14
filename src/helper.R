#remove individuals with any missing phenotype
deNA<-function(phe, naf = function(x) x==-9 | is.na(x))
{
	msk<-rep(x=F, times=nrow(phe));
	for(i in 3 : ncol(phe))
		msk = msk | naf(phe[,i]);
	phe<-phe[!msk,];
}

align_idv<-function(gmx, idv, phe, odr=T)
{
	#initialize genotype and phenotype indices.
	gdx <- 1:nrow(idv);
	pdx <- 1:nrow(phe);
	
	# align phenotype to genotype by individual ID
	idx <- match(idv$IID[gdx], phe$IID[pdx], 0L);
	idx <- idx[idx > 0L];
	pdx <- idx;
	
	# align genotype to phenotype by individual ID
	idx <- match(phe$IID[pdx], idv$IID[gdx], 0L);
	idx <- idx[idx > 0L];
	gdx <- idx;
	
	# sort by family ID
	if(odr)
	{		
		idx <- sort.int(x = idv$FID[gdx], index.return = T);
	 	idx <- idx$ix;
	 	pdx <- pdx[idx];
	 	gdx <- gdx[idx];
	}
	# return
	list(gdx=gdx, pdx=pdx);
}

# get error message from try-error object
err_msg<-function(try_error)
{
	sub('\\(converted from warning\\) ', '', attr(try_error, 'condition')[['message']])
}
# data file root name.
rut<-"out/rv1";
sex<-NULL;

if(!exists('fix', inherits = F))
{
	fix<-paste(rut, 'phe', sep='.');
	fix<-read.table(fix, nrows = 1, header = F, sep ='\t', stringsAsFactors = F)
	fix<-unlist(fix);
	fix<-as.list(fix[0L:-6L]);
}

# prepare file surfix
if(is.null(sex))
{
	sfx<-'sxc';	
}else
{
	sfx<-paste('sx', sex, sep='');
}

# start analysis
while(length(fix)>0)
{
	rsp<-fix[[1]];
	whr<-sprintf('%s.%s.%s.%s', rut, rsp, sfx, 'out');
	out<-read.table(whr, header = T, sep = '\t');
	rpt<-data.frame(
		SNP=out$SNP, CHR=out$CHR, POS=out$POS, EFF_ALL=out$A1, NONEFF_ALL=out$A2,
		BETA   = round(out$EST, 4L),
		SE     = round(out$SE,  4L),
		P      = format(out$P, digits = 3L, scientific = T, trim=T),
		AF_COD = round(out$MAF, 4L), 
		HWE=NA, IMP=0L, INFO=1L, INFO_TYPE=NA, N_EFF=out$NEF, # later, this would be NES
		row.names = NULL, stringsAsFactors = F);
	whr<-sprintf('%s.%s.%s.%s', rut, rsp, sfx, 'rpt');
	write.table(rpt, whr, row.names = F, quote = F, sep = '\t');
	
	rm(whr, rsp, out, rpt);
	fix[[1]]<-NULL;
}

if(exists('rut', inherits = F)) rm(rut);
if(exists('sex', inherits = F)) rm(sex);
if(exists('sfx', inherits = F)) rm(sfx);

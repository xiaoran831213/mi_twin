## prepare phase 2 twin genotype data
cd dat/p2.0

## extract genotype columns, drop R and Theta value from GenomeStudio
cut g.txt -f1,2,3,$(seq -s, 4 3 $(head -n 1 g.txt | wc -w)) > g.nox

## transform AA, BB, Ab and NC coding to dosage values
sed g.nox -e '1,1 s/\.GType//g' -e '2,$ s/AA/0/g' -e 's/BB/2/g' -e 's/AB/1/g' -e 's/NC/3/g' > g.dsg.0

## transform X, Y, XY and MT chromosome to numeric coding
sed g.dsg.0 -e 's/\tX\t/\t23\t/' -e 's/\tY\t/\t24\t/' -e 's/\tXY\t/\t25\t/' -e 's/\tMT\t/\t26\t/' > g.dsg.1

## sort the genome features by chromosome and basepair loacation
sort -k 2g -k 3g g.dsg.1 > g.dsg.2

## split off genome feature map dosage values matrix
cut -f1,2,3 g.dsg.2 > g.dsg.map && cut -f4- g.dsg.2 > g.dsg.mtx



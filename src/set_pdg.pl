use Getopt::Long;
my $pedPath = '';    # I, Plink PED path
my $outPath = '';    # O, Plink PED path.
my $setExpr = '';

GetOptions(
	"ped=s" => \$pedPath,
	"out:s" => \$outPath,
	"exp=s" => \$setExpr
  )
  or die("Error in command line arguments\n");

if ( !$outPath )
{
	( $outPath = $pedPath ) =~ s/\.[^.]+$//;
}
if ( $outPath !~ /\.[^.]+$/ )
{
	$outPath = $outPath . '.pdg';
}

# twin id copy map, k: sib_1 (existing), v: sib_2 (copy to)
open my $fPed, "<", $pedPath or die "Could not open $pedPath: $!";
open my $fOut, ">", $outPath or die "Could not open $outPath: $!";
my $phe = 2;
while (<$fPed>)
{
	my @idv = split( /\s/, $_, 7 );    # individual description
	my $fid = $idv[0];
	my $iid = $idv[1];
	my $pid = $idv[2];
	my $mid = $idv[3];
	my $sex = $idv[4];
	my $phe = $idv[5];
	my $gno = $idv[6];
	if ($setExpr)
	{
		eval $setExpr;
	}
	my $pdg = join( " ", $fid, $iid, $pid, $mid, $sex, $phe);
	$pdg =~ s/\s+/ /g;
	$pdg =~ s/^\s|\s$//;
	print $fOut $pdg;
	if($gno)
	{
		print $fOut ' '.$gno;
	}
	else
	{
		print $fOut "\n";
	}
}

close $fOut;
close $fPed;


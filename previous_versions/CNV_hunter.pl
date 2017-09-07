#! usr/share/bin/perl -w
# 
#
# This script analyses CNV files for identity and rare copy number identification using Fishers exact test
# 
# Written by Roel Sterken June18 2009 for Sim1 ;-)
#

use strict;
use lib "/Users/simone-sanna-cherchi/lib/";
#use Text::NSP::Measures::2D::Fisher2::left;
#use lib "/Library/Perl/5.8.8/";
use Text::NSP::Measures::2D::Fisher2::right;

=head1

CNV_hunter.pl -case <Case rawcnv file> -control <Control rawcnv file> -pc <% identity for call e.g. 90> -out <outfile> -dump <pairwise matches output>

The Fishers exact test per CNV is a left-sided test for the contingency table

	  Case	Control
   
Present    a	   b  	| a+b

Absent     x	   y	| x+y
	________________

	  a+x 	  b+y

	

=cut


#Check if there are arguments given, otherwise print usage-message via the subroutine "usage"
&usage("\n\tError -> please submit all parameters \n") if(scalar(@ARGV)<10);
my %params;
my ($Casefile,$Controlfile,$percent,$Outfile,$Dumpfile);

#check if output files exist in the current directory						


%params=@ARGV;
$params{"-case"} || die "Parameter \"-case\" is missing.Please specify a case file.\n";
$params{"-control"} || die "Parameter \"-control\" is missing.Please specify a control file.\n";
$params{"-pc"} || die "Parameter \"-pc\" is missing.Please specify percentage of identity to call.\n";
$params{"-out"} || die "Parameter \"-out\" is missing.Please specify an output file.\n";
$params{"-dump"} || die "Parameter \"-dump\" is missing. Please specify a pairwise matching output file.\n";

$Casefile = $params{"-case"};
$Controlfile = $params{"-control"};
$percent = $params{"-pc"};
$percent = ($percent)/100; 
$Outfile = $params{"-out"};
$Dumpfile = $params{"-dump"};

if(-e $Outfile or -e $Dumpfile){
	print 'ERROR: output files already exist',"\n"; exit; 
}

open (CASE, $Casefile) || die "Couldn't open $Casefile\n";
my @Cases = <CASE>;
close CASE;
open (CONTROL, $Controlfile)|| die "Couldn't open $Controlfile\n";
my @Controls = <CONTROL>;
close CONTROL;

# Step 1 -> I make a hash for all the controls by sample / cn / Chr / Pos
# Like this I can count the presence/absence of cnvs for the fisher test 
# and I can limit the number of tests to make by identical cn

my %CASES = &Hashit(\@Cases);
my %CONTROLS = &Hashit(\@Controls);

# Now I start looking for identical CNVs within the cases
# I test sample per sample. I first check if the cn# is identical, 
# then if the chrom is identical, then the positions are being matched
# To make sure I test every CNV, and output every CNV, I use a for loop on the array

my ($ChrCAS,$Adjust,$Start,$Uplimit,$InlimitUP,$Stop,$Downlimit,$InlimitDown,$CNV,$Type);

open (TESTOUT,">> $Outfile");
print TESTOUT "Case_CNV\tPos_Case\tNeg_Case\tPos_Control\tNeg_Control\tFisher_p\tNI_Bigger_Cases\tNI_Smaller_Cases\tNI_Big&SmallCases\tNI_Bigger_Controls\tNI_Smaller_Control\tNI_Big&SmallControl\n";

for (my $x = 0; $x < scalar @Cases;$x++){
	my $Caseline = $Cases[$x];
	chomp $Caseline;
	if ($Caseline =~ m/^chr/i) {
		my @CaseData = split(/\s+/,$Caseline);
		my ($State,$cn) = split(/,/,$CaseData[3]);
		# The subscript &Hunter will count how many samples have identical hits, non identical overlaps
		# and cases without identical overlap. It needs to get the "type" of comparison for the dump
		my ($PresCase,$AbsCase,$NI_BIGCases,$NI_SMALLCases,$NI_BSCases) = &Hunter(\%CASES,$Caseline,"Case-Case",$Dumpfile,$percent);
		my ($PresCont,$AbsCont,$NI_BIGControl,$NI_SMALLControl,$NI_BSControl) = &Hunter(\%CONTROLS,$Caseline,"Case-Control",$Dumpfile,$percent);
		
		# ... now up to the fisher exact test ...
		my $npp = $PresCase+$AbsCase+$PresCont+$AbsCont;
		my $n1p = $PresCase+$PresCont;
		my $np1 = $PresCase+$AbsCase;
		my $n11 = $PresCase;
		#my $leftFisher = Text::NSP::Measures::2D::Fisher2::left->new();
		my $rightFisher_value = calculateStatistic( n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
		#if((my $errorCode = $leftFisher->getErrorCode())){
    		#	print STDERR $errorCode." - ".$leftFisher->getErrorMessage();
  		#}
 		#else{
    		#	print $leftFisher_value;
  		#}
		print TESTOUT "$Caseline\t$PresCase\t$AbsCase\t$PresCont\t$AbsCont\t$rightFisher_value\t$NI_BIGCases\t$NI_SMALLCases\t$NI_BSCases\t$NI_BIGControl\t$NI_SMALLControl\t$NI_BSControl\n";
		
	}
}

exit;


#------------------------------------------------------------
sub usage( $ )
{
	print "$_[0]\n";
	system("pod2text $0"); #the command pod2text calls the commando given between the =head and =cut 
	exit(1);
}
#------------------------------------------------------------
sub roundup ( $ )
{
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}
#------------------------------------------------------------
sub Hashit ( $ )
{
	my @Input = @{$_[0]};
	chomp @Input;
	my %HASH;
	foreach my $CNVline(@Input){
		#All the CNV files have to start with a CNV identifier, which has the structure ChrX:Start-Stop
		chomp $CNVline;
		if ($CNVline =~ m/^chr/i) {
			my @CNVData = split(/\s+/,$CNVline);
			my @State = split(/,/,$CNVData[3]);
			my $CNV = $State[1];
			my ($Chr,$Pos) = split(/[:]/,$CNVData[0]);
			my $CNVID = $CNVData[0];
			my $Sample = $CNVData[4];
			$HASH{$Sample}{$CNV}{$Chr}{$Pos}=$CNVline;
		}
	}
return (%HASH);
}
#-------------------------------------------------------------
sub Hunter ( $ $ $ $ )
{

	my %HASH = %{$_[0]};
	my $line = $_[1];
	my $Type = $_[2];
	my $Dump = $_[3];
	my $percent = $_[4];
	open (OUT, ">> $Dump");
	my @lineData = split(/\s+/,$line);
	my ($Chr,$Start,$Stop) = split(/[:-]/,$lineData[0]);
	my @State = split(/,/,$lineData[3]);
	$CNV = $State[1];
	my $Sample = $lineData[4];
	#The Positive and negative counters will keep track of the samples found with an identical match (positive) or not
	my $Positive = 0;
	my $Negative = 0;
	my $NI_BIG = 0;
	my $NI_SMALL = 0;
	my $NI_BS = 0;
	my $sampleNo = 0; 
	my ($TestSample,$Identical,$NIoverlap_SM,$NIoverlap_BIG);
	foreach $TestSample (sort keys %HASH){
		$Identical = 0;
		$NIoverlap_SM = 0;
		$NIoverlap_BIG = 0;
		$sampleNo++;
		if ($TestSample eq $Sample){
			$Positive++;
			next;
		}
		#I have to add an if-not to avoid counting the line itself
		if (!($TestSample eq $Sample)){
			# I can limit the testing to all positions with identical Chr & cnv
			foreach my $TestPos (sort keys %{$HASH{$TestSample}{$CNV}{$Chr}}){
				my $match = $HASH{$TestSample}{$CNV}{$Chr}{$TestPos};
				# The &Typing routine will test how the TestPos is located to the CNV
				my $kind = &Typing($lineData[0],$percent,$TestPos);
				#Now that the kind is determined, I can start doing the counting
				if (!defined $kind){
					print "Something isn't working on $lineData[0] ,$percent ,$TestPos\n";
				}
				if ($kind eq "Identical"){
					$Identical++;
					print OUT "$line\tIdentical\t$Type\t$match\n";
				}
				elsif ($kind eq "NIoverlap_bigger"){
					$NIoverlap_BIG++;
					print OUT "$line\tNIoverlap_bigger\t$Type\t$match\n";
				}
				elsif ($kind eq "NIoverlap_smaller"){
					$NIoverlap_SM++;
					print OUT "$line\tNIoverlap_smaller\t$Type\t$match\n";
				}
			}
			#Now if I have found at least one Identical or NIoverlap, I have to count that
			if ($Identical > 0){
			$Positive++;
			}
			if (($NIoverlap_BIG > 0) && ($NIoverlap_SM == 0)){
			$NI_BIG++;
			}
			if (($NIoverlap_BIG == 0) && ($NIoverlap_SM > 0)){
			$NI_SMALL++;
			}
			if (($NIoverlap_BIG > 0) && ($NIoverlap_SM > 0)){
			$NI_BS++;
			}
			if ($Identical == 0){
			$Negative++;
			}
		}
	}
	# A short internal control to see if nothing is counted twice	
	if ($sampleNo != ($Positive + $Negative)){
		print STDOUT "Testing $line for $Type in $Sample The scripts counts $Positive positive and $Negative negative samples in $sampleNo samples ... this doesn't compute! Check the script\n";
		die;
	}
	close OUT;

return ($Positive,$Negative,$NI_BIG,$NI_SMALL,$NI_BS);

}

#----------------------------------------------------------------
sub Typing ( $ $ $)
{
	my ($Chr,$Start,$Stop) = split(/[:-]/,$_[0]);
	my $percent = $_[1];
	my $Size = $Stop-$Start+1;
	my ($TestStart,$TestStop) = split (/-/,$_[2]);
	my $SizeTest = $TestStop-$TestStart+1;
	my ($Uplimit,$InlimitUP,$Downlimit,$InlimitDOWN,$Adjust,$ratio,$NewPercent,$kind,$pMargins);
	if (($TestStart == $Start) && ($TestStop == $Stop)){
		$kind = "Identical";
	}
	if (($TestStop <= $Start) || ($TestStart >= $Stop)){
		 $kind = "Null";
	}
	elsif (($SizeTest < ($percent*$Size)) && ((($TestStop <= $Stop) && ($TestStop > $Start)) || (($TestStart >= $Start) && ($TestStart < $Stop)))){
		# If the testCNV is smaller than x% the CNV, it can only be a NIoverlap, but I have to test for overlap
		$kind = "NIoverlap_smaller";
	}
	elsif ($SizeTest > (1/$percent)*$Size) {
		#Now I have to test if the overlap between CNV and TestCNV is smaller than x% or not 
		$InlimitUP = $Start + ((1-$percent)*$Size);
		$InlimitDOWN = $Stop - ((1-$percent)*$Size);
		if (($TestStart < $Start) && ($TestStop > $Stop)){
			$kind = "NIoverlap_bigger";
		}
		elsif (($TestStart >= $Start) && ($TestStart <= $InlimitUP) && ($TestStop > $Stop)){
			$kind = "NIoverlap_bigger";
		}
		elsif (($TestStart < $Start) && ($TestStop >= $InlimitDOWN) && ($TestStop <= $Stop)){
			$kind = "NIoverlap_bigger";
		}
		elsif (($TestStart < $Start) && ($TestStop < $InlimitDOWN) && ($TestStop >= $Start)){
			$kind = "NIoverlap_smaller";
		}
		elsif (($InlimitUP < $TestStart) && ($TestStart <= $Stop) && ($TestStop > $Stop)){
			$kind = "NIoverlap_smaller";
		}
		else { print STDOUT "Here - $SizeTest = Sizetest, $Size = Size, $percent = percent \n There is something wrong in the script, as it doesn't know what to do with $Chr\:$Start\-$Stop and $TestStart\-$TestStop\n";die}
	}
	elsif (($SizeTest <= ((1/$percent)*$Size)) && ($SizeTest >= ($percent*$Size))){
		# In this case there is a possibility to be identical. The CNV/TestCNV size ratio will determine the positions of the overlap limits
		# The Newpercentage is needed to set the boundaries of overlap, but is adjusted from the ratio
		# I need to make a distiction if the TestCNV is bigger or smaller than the CNV to see which one has to be 90% overlapping with the other
		if ($SizeTest < $Size){
			$ratio = (1-$SizeTest/$Size);
			$pMargins = (1-$percent);
			$NewPercent = ($pMargins-$ratio);
			$Adjust = ($NewPercent*$Size);
			$Adjust = &roundup($Adjust);
			$Uplimit = $Start-$Adjust;
			$Downlimit = $Stop+$Adjust;
			if (($TestStart >= $Uplimit) && ($TestStop <= $Downlimit)){
			$kind = "Identical";
			}
			elsif (($TestStart >= $Uplimit) && ($TestStart < $Stop) && ($TestStop > $Downlimit)){
			$kind = "NIoverlap_smaller";
			}
			elsif (($TestStart < $Uplimit) && ($TestStop <= $Downlimit) && ($TestStop > $Start)){
			$kind = "NIoverlap_smaller";
			}
			else {print STDOUT "1 - $SizeTest = Sizetest, $Size = Size, $ratio = ratio, $NewPercent = NP, $percent = percent \n There is something wrong in the script, as it doesn't know what to do with $Chr\:$Start\-$Stop and $TestStart\-$TestStop\n";
				die;}
		}
		if ($SizeTest > $Size){
			# Now I'm actually gonna test if the CNV is 90% overlapping the TestCNV
			# But I have to see what percentage of the CNV is covered to see if we're talking about a smaller overlap or a bigger overlap
			# When < x% of the CNV is covered, it's smaller and this would mean that in the cases it's possible the extra chunk of CNV is actually causing the phenotype
			$ratio = (1-$Size/$SizeTest);
			$pMargins = (1-$percent);
			$NewPercent = ($pMargins-$ratio);
			$Adjust = ($NewPercent*$SizeTest);
			$Adjust = &roundup($Adjust);
			$Uplimit = $TestStart-$Adjust;
			$Downlimit = $TestStop+$Adjust;
			if (($Start >= $Uplimit) && ($Stop <= $Downlimit)){
			$kind = "Identical";
			}
			elsif (($Start >= $Uplimit) && ($Start < $TestStop) && ($Stop > $Downlimit)){
			$kind = "NIoverlap_smaller";
			}
			elsif (($Start < $Uplimit) && ($Stop <= $Downlimit) && ($Stop > $TestStart)){
			$kind = "NIoverlap_smaller";
			}
			else {print STDOUT "2 - $SizeTest = Sizetest, $Size = Size, $ratio = ratio, $NewPercent = NP, $percent = percent \n There is something wrong in the script, as it doesn't know what to do with $Chr\:$Start\-$Stop and $TestStart\-$TestStop\n";die;}
		}
	}
	else {print STDOUT "3 - $SizeTest $TestStart $TestStop = Sizetest, $Size $Start $Stop = Size, $ratio = ratio, $NewPercent = NP, $percent = percent \n There is something wrong in the script, as it doesn't know what to do with $Chr\:$Start\-$Stop and $TestStart\-$TestStop\n";die;}
	return ($kind);	
}



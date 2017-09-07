#!/usr/bin/perl
#
# This script analyses CNV files to identify rare CNV's.
# By David A Fasel daf2139<at>columbia.edu 11/2012
# Adapted from code by Roel Sterken 6/2009
# 
# Basic usage:
# ============
# perl CNV_hunter2.pl -case CaseFile.txt
#
# Usage with optional parameters:
# ===============================
# perl CNV_hunter2.pl -case CaseFile.txt -control ControlFile.txt -out OutputPath -dump DumpPath -pc 99 -freq 10
#
# -case (required)
#  the path to the whitespace delimited Cases file (see "Format of Input Files" below)
#
# -control 
#  Optionally specifies a control file to check each case against.
#  If omitted, the cases will only be checked against themselves
#
# -out 
#  Specifies the path for the output files. Defaults to the case file path
#  plus ".matches.txt" i.e. "cases.txt.matches.txt" 
#
# -pc 
#  Specifies the percent for the match to be considered "identical";  
#  defaults to 100 (see notes below for types of matches)
#
# -freq
#  The dump file will show all types of matches for each CNV in the cases file;
#  This will potentially create a large dump file.  Setting this value will output 
#  only matches when the match frequency is less than that percent. Default is 100.
#
#
# Format of Input File(s):
# ========================
# The input file(s) must be whitespace delimited and in this format:
# <chromosome:basepair_range> <number of snps> <length of range> <copy number type> <sample id> <starting SNP> <ending SNP> <confidence score>
# Example:
# chr16:32438272-32505195 91 66924 state2,cn=1 7 cnvi0020587 cnvi0020628 86.717
# 
# Notes: 
# -- "number of snps"; "length of range"; "starting SNP"; "ending SNP" and "confidence score" 
#     are not used but must contain a value
# 
# -- "copy number type" must contain <cn=#>, any other information in this field is ignored
#
# -- "sample id" identifies the individual and can be any type of string, i.e. "7" or "person_1".  It cannot be blank or 0.
# 
# Format of Output File:
# ======================
# The format of the output file is the **entire line from the case file** plus:
# * Cases Identical: Number of identical matches in the Case file.  
#                    There could be more than one match per sample id if <pc> is set to less 
#                    than 100 (see match descriptions below).
# * Cases Non-match: Total samples that do not have at least one identical match.
# * Cases Frequency: Frequency of cases (samples with at least one identical match / total sample id's)
# * Controls Identical: see above
# * Controls Non-match: see above
# * Controls Frequency: see above
# * Cases Non-identical Smaller: the case CNV is contained within the comparison CNV
# * Cases Non-identical Larger: the comparison CNV is contained within the case CNV
# * Cases Non-identical: none of the matches above apply, however, there is at least some overlap between the CNV's 
# * Controls Non-identical Smaller: see above
# * Controls Non-identical Larger: see above
# * Controls Non-identical: see above
# * Case identical matches: comma separated list of sample ID's that contain at least 1 identical match
# * Control identical matches: see above
# 
# Format of Dump File:
# ====================
# TODO
#
# Match Types:
# ============
# Note: Matches are not performed within the same sample ID
#
# "Identical" - matches within the <pc> percent specified.  This checks for a two-way match
#  where both the CNV and the comparison CNV are both at least percent <pc> matches of each other
# 
#     |----------------|          [test start/end]
#     TS               TE
#      |---------------|         [case start/end]
#      S                E
#
# "Non-identical larger" - the case CNV is contained within the comparison CNV
# 
#     |----------------|          [test start/end]
#     TS               TE
#         |-------|               [case start/end]
#         S       E
#
# "Non-identical smaller" - the comparison CNV is contained within the case CNV
# 
#     |----------|          [test start/end]
#     TS         TE
#   |-----------------|     [case start/end]
#   S                 E
#
# "Non-identical" - none of the matches above apply, 
#  however, there is at least some overlap between the CNV's 
#  
#     |----------|          [test start/end]
#     TS         TE
#           |---------|     [case start/end]
# 



use strict;
use warnings;
use Time::HiRes;
#use Text::NSP::Measures::2D::Fisher2::right;

=head1 USAGE
=====
CNV_counter.pl 
 -case <case cnv file, white space delimited>
 -control <control cnv file, white space delimited (defaults to case file> 
 -pc <percent for match to be considered "identical" (defaults to 100)> 
 -out <Output File Path (defaults to Casefile.matches.txt)>
 -dump <Output Dump File Path (defaults to Output.matches_dump.txt)>  
 -freq <If frequency of matches is greater than this number, suppress listing 
        them in output & dump file (defaults to 100 - no filtering)>

=cut

#Check if there are arguments given, otherwise print usage-message via the subroutine "usage"
&usage("\n\tError -> please submit all parameters \n") if(scalar(@ARGV) < 1);

# check that all arguments are present
my %params = @ARGV;
$params{"-case"} || die "Parameter \"-case\" is missing.Please specify a case.\n";
# $params{"-control"} || die "Parameter \"-control\" is missing.Please specify a control file.\n";
# $params{"-out"} || die "Parameter \"-out\" is missing.Please specify an output file.\n";

my $casePath = $params{"-case"};
my $controlPath = $params{"-control"} || $params{"-case"};;
my $outPath = $params{"-out"} || "$casePath.matches";
my $dumpPath = $params{"-dump"} || $outPath . "_dump";
$outPath = $outPath . ".txt";
$dumpPath = $dumpPath . ".txt";
my $PERCENT = $params{"-pc"} || 90;
$PERCENT = ($PERCENT)/100; 
my $FREQ_FILTER = $params{"-freq"} || 100;
$FREQ_FILTER = ($FREQ_FILTER)/100;
my  ($Case_counter, $Control_counter, $Counter); 

# check if output files exist in the current directory
if(-e $outPath){
	print 'ERROR: output files already exist',"\n"; exit; 
}

#my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
my $start = time();
my ($sec, $min, $hour) = localtime();
printf "Started: %d:%02d:%02d.\n", $hour, $min, $sec;



# Make a hash for all the cases by SampleID / cn / Chromosome / Position
open (CASE, $casePath) || die "Couldn't open $casePath\n";
my @casesFile = <CASE>;
close CASE;
chomp (@casesFile);
my %cases = &HashCNV("Case", \@casesFile);

# do the same for controls if a control file was specified
open (CONTROL, $controlPath)|| die "Couldn't open $controlPath\n";
my @controlsFile = <CONTROL>;
close CONTROL;
chomp (@controlsFile);
my %controls = &HashCNV("Control", \@controlsFile);
print "Found $Case_counter valid case lines and $Control_counter valid control lines.\n";

# use Data::Dumper;
# my $temp = $cases{"7"};
# print Dumper $temp;

open (OUTPUTFILE,"> $outPath");
# add Fisher   \tFisher_p
print OUTPUTFILE "Position\tNum_SNPs\tBP_range\tCNV_State\tSample_ID\tStart_SNP\tEnd_SNP\tConf_Score";
print OUTPUTFILE "\tPos_Case\tNeg_Case\tFreq_Case\tPos_Control\tNeg_Control\tFreq_Control";
print OUTPUTFILE "\tNI_Bigger_Cases\tNI_Smaller_Cases\tNI_Cases";
print OUTPUTFILE "\tNI_Bigger_Controls\tNI_Smaller_Control\tNI_Control";
print OUTPUTFILE "\tID_Matches_Cases\tID_Matches_Controls\n";

open (DUMPFILE,"> $dumpPath");


for my $caseline (@casesFile) 
{
	chomp($caseline);
	my $totalCases = keys %cases;
	my $totalCont = keys %controls;
	
	my ($identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical, $identicalMatchIDs, $dumplines) = 
	    &findMatches(\%cases, $caseline, "Case-Case");
	# if the line was bad, &findMatches returns "badline", so go to next line    
 	if ($identical eq "badline"){
 		next;
 	}
	    
	# TODO: currently if control path is not specified, the case file is also used as the control.
	# May want to change the following code to omit the control if it is not specified
	# (also add the if statement for this when creating %controls above)
	my ($identicalCont, $nonIdentSmallerCont, $nonIdentLargerCont, $nonIdenticalCont, $identicalMatchIDsCont, $dumplinesCont) = 
	    &findMatches(\%controls, $caseline, "Case-Control");

	#frequency is determined by number of unique samples that have at least one match divided by total samples
	my ($caseFreq, $controlFreq) = ((scalar @$identicalMatchIDs / $totalCases), (scalar @$identicalMatchIDsCont / $totalCont));
	# nonMatches are samples that do not contain at least one match
	my ($nonMatch, $nonMatchCont) = (($totalCases - scalar @$identicalMatchIDs), ($totalCont - scalar @$identicalMatchIDsCont));
	
	# prepare the list of identical matches
	if (!$caseFreq) {@$identicalMatchIDs = ("none") }
	elsif ($caseFreq > $FREQ_FILTER){@$identicalMatchIDs = ("exceeds_freq_filter")}
	
	if (!$controlFreq) {@$identicalMatchIDsCont = ("none") }
	elsif ($controlFreq > $FREQ_FILTER) {@$identicalMatchIDsCont = ("exceeds_freq_filter")}
	
	#print $nonMatchCont;
	print OUTPUTFILE "$caseline\t$identical\t$nonMatch\t$caseFreq\t$identicalCont\t$nonMatchCont\t$controlFreq\t";
	print OUTPUTFILE "$nonIdentSmaller\t$nonIdentLarger\t$nonIdentical\t$nonIdentSmallerCont\t$nonIdentLargerCont\t$nonIdenticalCont\t";
	print OUTPUTFILE join(",", sort @$identicalMatchIDs) . "\t" . join(",", sort @$identicalMatchIDsCont) . "\n";
	
	
	# TODO dumpfile header
	if ($caseFreq <= $FREQ_FILTER)
	{
		print DUMPFILE @$dumplines;
	}
	if ($controlFreq <= $FREQ_FILTER) {
		print DUMPFILE @$dumplinesCont;
	}
	
	$Counter++;
	print "\rProcessed $Counter of $Case_counter case lines.";
	$| = 1;
	
}
($sec, $min, $hour) = localtime();
#printf "\nFinished: %d:%02d:%02d.\n", $hour, $min, $sec;
$sec = time() - $start;
$hour = int($sec / 3600);
$min = int(($sec % 3600) / 60);
$sec = $sec % 60;
printf "\nTime ellapsed: %dh %02dm %02ds.\n", $hour, $min, $sec;

close OUTPUTFILE;
close DUMPFILE;
exit;

sub HashCNV ($)
{
	my ($file_type, @lines) = ($_[0], @{$_[1]}) ;
	chomp @lines;
	my (%HASH, $lcounter);	
	
	foreach my $CNVline (@lines) 
	{
		++$lcounter;
		
		#skip lines that only consist of whitespace characters;
		if ($CNVline =~ /^\s*$/){
			next;
		}
		
		my ($sample, $state, $chr, $pos) = ParseLine($CNVline);
		if ($sample && $chr && $pos && $state){ 
		    $HASH{$sample}{$state}{$chr}{$pos}=$CNVline;
		    if ($file_type eq "Case") {
		    	$Case_counter++; 
		    }
		    else { 
		    	$Control_counter++;
		    }
	    }
	    # print a warning if the line is bad (unless it's the header)
	    elsif ($lcounter != 1) {
	    	print "Warning:\n  $file_type file line $lcounter is missing data and was skipped.\n";
	    	print "  None of these fields can be 0 or empty.\n";
	    	print "  Sample:$sample, State:$state, Chr:$chr, Pos:$pos\n";
	    }
	    
	}
    return (%HASH);
}

sub findMatches ($ $ $ $)
{
	my %hash = %{$_[0]};
	my ($sample, $state, $chr, $position) = ParseLine($_[1]);
	if (!$sample || !$state || !$chr || !$position){ 
		return "badline";
	}
	
	my ($start, $end) = split(/-/, $position);
	my ($testStart, $testEnd);
	my ($identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical, $nonMatch) = (0, 0, 0, 0, 0);
	my (@matchesForDumpfile, $matchtype);
	my %identicalMatchIDs;
	# for each sample id (which represents a unique individual), we check to see
	# if there are any matches. 
	for my $key_sample (keys %hash) 
	{
		#print "$sample; $key_sample\n";
	    # don't check the same person for matches
		if ($key_sample eq $sample)
		{
			#uncomment this line if you don't want to perform matches in the same sample id
			#next;  
		}
		# check each position within each individual for a match
		for my $key_position (keys %{$hash{$key_sample}{$state}{$chr}})
		{
			($testStart, $testEnd) = split(/-/, $key_position);
			#print "$_[3]: $sample: $start-$end; $key_sample:$testStart-$testEnd\n";
			#print "$key_position\n";
			
			# if there's no overlap at all, go to the next position
			if ($start > $testEnd || $end < $testStart) 
			{
				$nonMatch++;
				#print "nm\n";
				next;	
			}
			#print "match\n";		
			# check for an 'exact' match (a match within the percent specified)
			# the overlap between the CNV regions is found by sorting the base pair
			# positions and subtracting the 2nd largest value from the 3rd largest value.
			my @bp = ($start, $end, $testStart, $testEnd);
			@bp = sort {$a <=> $b} @bp;
			my $overlap = $bp[2] - $bp[1];
			# divide the overlap between the CNV's by the length of each CNV to
			# determine the match percentage.
			my $match_case = $overlap / ($end - $start);
			my $match_test = $overlap / ($testEnd - $testStart);
			#print "@bp; $PERCENT; $match_case; $match_test\n";
			if ($match_case >= $PERCENT && $match_test >= $PERCENT)
			{
				$identical++;
				$matchtype = "Identical";
				$identicalMatchIDs{$key_sample} = 1;
			}
			
			# If not an 'exact' match, check for the three types of matches below:
			
			# CNV position is contained within the control CNV
			#     |----------------|          [test start/end]
			#     TS               TE
			#         |-------|               [case start/end]
			#         S       E
			elsif ($start >= $testStart && $end <= $testEnd)
			{
				$nonIdentLarger++;
				$matchtype = "Non_Identical_Larger";
			}
			# CNV position is greater than the control CNV 
			#     |----------|          [test start/end]
			#     TS         TE
			#   |-----------------|     [case start/end]
			#   S                 E
			elsif ($start <= $testStart && $end >= $testEnd)
			{
				$nonIdentSmaller++;
				$matchtype = "Non_Identical_Smaller";
			}
			# Some overlap, but doesn't meet the specified match percentage 
			#     |----------|          [test start/end]
			#     TS         TE
			#           |---------|     [case start/end]
			#           S         E
			else 
			{
				$nonIdentical++;
				$matchtype = "Non_Identical";
			}
			
			# store information about the match to output to the dumpfile
			my @caseline = split(/\s+/, $_[1]);
			my @testline = split(/\s+/, $hash{$key_sample}{$state}{$chr}{$key_position});
			my $dumpline = "$caseline[0]\t" . ($end - $start) . 
				"\t$caseline[3]\t$caseline[4]\t" .
				"$matchtype\t$_[2]\t$testline[0]\t" . ($testEnd - $testStart) . 
				"\t$testline[3]\t$testline[4]\t$testline[7]\n";
			push (@matchesForDumpfile, $dumpline);
		}

	}
	my @identicalMatchIDs = keys %identicalMatchIDs;
	return ($identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical, \@identicalMatchIDs, \@matchesForDumpfile);
}

sub ParseLine ($) 
{
	my @fields;
	my ($sample, $state, $chr, $pos) = ('', '', '', '');
	my $line = $_[0];
	chomp $line;
	# All the CNV files have to start with a CNV identifier, which has the 
	# structure ChrX:Start-Stop
	# split the line into fields by whitespace
	@fields = split(/\s+/, $line);
	if ($fields[0]) { 
		($chr,$pos) = split(/[:]/, $fields[0]); 
		if (!$pos) {$pos = '';}
	}
	if ($fields[3]) { 
		$state = ($fields[3] =~ m/(cn=\d)/i) ? $1 : "";
	}
	if ($fields[4]) { 
		$sample = $fields[4];
	}
		
	return ($sample, $state, $chr, $pos);
}

#------------------------------------------------------------
sub usage( $ )
{
	#the command pod2text calls the commando given between the =head and =cut 
	print "$_[0]\n";
	system("pod2text $0"); 
	exit(1);
}

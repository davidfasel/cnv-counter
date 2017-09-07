#!/usr/bin/perl
#
# This script analyses copy number variation (CNV) files to identify rare CNV's.
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
# -case <file path>
#  Required. The path to the Cases file (see "Format of Input Files" below)
#
# -control <file path>
#  Optionally specifies a control file to check each case against.
#  If omitted, the cases will only be checked against themselves
#
# -out <file path>
#  Specifies the path for the output files. Defaults to the case file path
#  plus ".matches.txt" i.e. "cases.txt.matches.txt" 
#
# -pc <float>
#  Specifies the percent for the match to be considered "identical";  
#  defaults to 100 (see notes below for types of matches).
#  Since match percentages are not rounded, if you set <pc> to 70, matches 
#  of 0.69999999 will not be considered exact matches.  If you want this behavior, set <pc> to 69.5 (for example).  
#
# -freq <float>
#  Default is 100 (no filtering).
#  Use this filter to prevent bloated output and dump files. 
#  If the match frequency is greater than this value for a particular case:
#    -Suppress listing Id's of identical matches in the output file (see Format of Output File).  
#    -Suppress listing a line for each match in the dump file. 
#
#
# Format of Input File(s):
# ========================
# The input file(s) must be whitespace delimited and in this format:
# <chromosome:start-end> <# of snps> <# of bp's> <copy number type> <sample id> <first SNP> <last SNP> <confidence score>
# Example:
# chr16:32438272-32505195 91 66924 state2,cn=1 7 cnvi0020587 cnvi0020628 86.717
# 
# Notes: 
# -- "number of snps"; "length of range"; "starting SNP"; "ending SNP" and "confidence score" 
#     are not used but must contain a value
# 
# -- "copy number type" must contain <cn=#>, any other information in this field is ignored
#
# -- "sample id" identifies the individual and can be any type of string, i.e. "7" or "person_1".  
#     It cannot be blank or 0.
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
# <entire line from case file> 
#
# Match Types:
# ============
# Note: Matches are not performed within the same sample ID
#
# "Identical" - matches within the <pc> percent specified.  Checks for a two-way match
#  where both the CNV and the tested CNV are both at least <pc> matches of each other
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
# Fisher Test
# ===========
# Note: you must have the Fisher2 perl module installed in an included Perl directory.
# Module:  Text::NSP::Measures::2D::Fisher2::right
# The Fishers exact test per CNV is a right-sided test for the contingency table
# 
#          Case  Control    | Totals
#                           |
# Present   a     b         | a+b
#                           |
# Absent    x     y         | x+y
# __________________________|________
#                           |
# Totals    a+x   b+y       | a+b+x+y

use strict;
use warnings;
use Text::NSP::Measures::2D::Fisher2::right;

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

#Check if there are arguments given, otherwise print usage-message
&usage("\n\tError -> please submit all parameters \n") if(scalar(@ARGV) < 1);

# check that the case file path has been specified
my %params = @ARGV;
$params{"-case"} || die "Parameter \"-case\" is missing.Please specify a case.\n";

# get or create defaults for the rest of the params
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
my  ($Case_counter, $Control_counter, $Counter) = (0,0,0); 

# check if output files exist in the current directory
if(-e $outPath){
	print "Warning: output files already exist \n"; 
	print "Overwrite? y/n \n";
	my $ans = <STDIN>;
	exit if ($ans ne "y\n");
}

my $start = time();
my ($sec, $min, $hour) = localtime();
printf "Started: %d:%02d:%02d.\n", $hour, $min, $sec;

# Make a hash for all the cases by SampleID / cn / Chromosome / Position
open (CASE, $casePath) || die "Couldn't open $casePath\n";
my @casesFile = <CASE>;
close CASE;
chomp (@casesFile);
my %cases = &HashCNV("Case", \@casesFile);
print "Found $Case_counter valid case lines. \n";

# do the same for controls if a control file was specified
open (CONTROL, $controlPath)|| die "Couldn't open $controlPath\n";
my @controlsFile = <CONTROL>;
close CONTROL;
chomp (@controlsFile);
my %controls = &HashCNV("Control", \@controlsFile);
print "Found $Control_counter valid control lines.\n";

open (OUTPUTFILE,"> $outPath");

print OUTPUTFILE 
	"Position\tNum_SNPs\tBP_range\tCNV_State\tSample_ID\tStart_SNP\tEnd_SNP\tConf_Score".
	"\tPos_Case\tNeg_Case\tFreq_Case\tPos_Control\tNeg_Control\tFreq_Control\tFisher_p".
	"\tNI_Bigger_Cases\tNI_Smaller_Cases\tNI_Cases".
	"\tNI_Bigger_Controls\tNI_Smaller_Control\tNI_Control".
	"\tID_Matches_Cases\tID_Matches_Controls\n";

open (DUMPFILE,"> $dumpPath");

for my $caseline (@casesFile) 
{
	my $totalCases = keys %cases;
	my $totalCont = keys %controls;
	my ($identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical, $identicalMatchIDs, $dumplines) = 
	    &findMatches(\%cases, $caseline, "Case-Case");
	# if the line was bad, &findMatches returns "badline", so go to next line    
 	next if ($identical eq "badline");
	    
	# if control path is not specified, the case file is also used as the control.
	my ($identicalCont, $nonIdentSmallerCont, $nonIdentLargerCont, $nonIdenticalCont, $identicalMatchIDsCont, $dumplinesCont) = 
	    &findMatches(\%controls, $caseline, "Case-Control");

	# frequency is determined by number of unique samples that have 
	# at least one match, divided by total samples
	my $caseFreq = scalar @$identicalMatchIDs / $totalCases;
	my $controlFreq = scalar @$identicalMatchIDsCont / $totalCont;
	
	# nonMatches are samples that do not contain at least one match
	my $nonMatch = $totalCases - scalar @$identicalMatchIDs;
	my $nonMatchCont = $totalCont - scalar @$identicalMatchIDsCont;
	
	# prepare the list of identical matches
	if (!$caseFreq) {@$identicalMatchIDs = ("none") }
	elsif ($caseFreq > $FREQ_FILTER){@$identicalMatchIDs = ("exceeds_freq_filter")}
	if (!$controlFreq) {@$identicalMatchIDsCont = ("none") }
	elsif ($controlFreq > $FREQ_FILTER) {@$identicalMatchIDsCont = ("exceeds_freq_filter")}
	
	# Calculate the fisher exact test
	my $npp = $totalCases + $totalCont;        #$PresCase+$AbsCase+$PresCont+$AbsCont;
	my $n1p = $identical + $identicalCont;     #$PresCase+$PresCont;
	my $np1 = $identical + $nonMatch;          #$PresCase+$AbsCase;
	my $n11 = $identical;                      #$PresCase;
	my $rightFisher = calculateStatistic( n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);

	print OUTPUTFILE "$caseline\t$identical\t$nonMatch\t$caseFreq\t".
		"$identicalCont\t$nonMatchCont\t$controlFreq\t$rightFisher\t".
		"$nonIdentLarger\t$nonIdentSmaller\t".
		"$nonIdentical\t$nonIdentLargerCont\t$nonIdentSmallerCont\t$nonIdenticalCont\t".
		join(",", sort @$identicalMatchIDs)."\t". 
		join(",", sort @$identicalMatchIDsCont) . "\n";
	
	# output the dump file if frequency is less than the filter set by user
	print DUMPFILE @$dumplines if $caseFreq <= $FREQ_FILTER;
	print DUMPFILE @$dumplinesCont if ($controlFreq <= $FREQ_FILTER);
	
	
	($sec, $min, $hour) = localtime();
	my $sec = time() - $start;
	$hour = int($sec / 3600);
	$sec = $sec % 3600;
	$min = int($sec / 60);
	$sec = $sec % 60;


	# update the status in the terminal and flush the cache to force terminal to display text
	$Counter++;
	printf ("\rProcessed $Counter of $Case_counter case lines. Time ellapsed: %d:%02d:%02d. ", $hour, $min, $sec);
	$| = 1;
}


#printf "\nTime ellapsed: %dh %02dm %02ds.\n", $hour, $min, $sec;
print "\nDone.\n";

close OUTPUTFILE;
close DUMPFILE;
exit;

# ----------------------------------------------------------------------------
# %hash = &HashCNV( $file_type, \@lines )
# ----------------------------------------------------------------------------
sub HashCNV
{
	my ($file_type, @lines) = ($_[0], @{$_[1]}) ;
	chomp @lines;
	my (%hash, $lcounter);	
	
	foreach my $CNVline (@lines) 
	{
		++$lcounter;
		
		#skip lines that only consist of whitespace characters;
		next if ($CNVline =~ /^\s*$/);
		
		# Parse the line and if it's valid, put it in the hash and count it
		my ($sample, $state, $chr, $pos) = ParseLine($CNVline);
		if ($sample && $chr && $pos && $state){ 
		    $hash{$sample}{$state}{$chr}{$pos}=$CNVline;
		    $file_type eq "Case" ? $Case_counter++ : $Control_counter++;
	    }
	    # print a warning if the line is bad (unless it's the header)
	    elsif ($lcounter != 1) {
	    	print "Warning:\n $file_type file line $lcounter is missing data and was skipped.\n";
	    	print "  None of these fields can be 0 or empty.\n";
	    	print "  Sample:$sample, State:$state, Chr:$chr, Pos:$pos\n";
	    }
	}
    return (%hash);
}

# ----------------------------------------------------------------------------
# @values = &findMatches( \%hash, $line, $file_type )
# ----------------------------------------------------------------------------
sub findMatches
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
		# check each position within each individual for a match
		for my $key_position (keys %{$hash{$key_sample}{$state}{$chr}})
		{
			($testStart, $testEnd) = split(/-/, $key_position);
			
			# if there's no overlap at all, go to the next position
			if ($start > $testEnd || $end < $testStart) {
				$nonMatch++;
				next;	
			}
					
			# check for an 'exact' match (a match within the percent specified)
			# the overlap between the CNV regions is found by sorting the base pair
			# positions and subtracting the 2nd largest value from the 3rd largest value.
			my @bp = ($start, $end, $testStart, $testEnd);
			@bp = sort {$a <=> $b} @bp;
			my $overlap = $bp[2] - $bp[1];
			# divide the overlap between the CNV's by the length of each CNV to
			# determine the match percentage.
			#my $match_case = sprintf("%.3f", $overlap / ($end - $start));
			#my $match_test = sprintf("%.3f", $overlap / ($testEnd - $testStart));
			
			my $match_case = $overlap / ($end - $start);
			my $match_test = $overlap / ($testEnd - $testStart);
			
			if ($match_case >= $PERCENT && $match_test >= $PERCENT){
				$identical++;
				$matchtype = "Identical";
				$identicalMatchIDs{$key_sample} = 1;
			}
			# If not an 'exact' match, check for the three types of matches below:
			# For definitions of these matches, see comments at the top of this file.
			elsif ($start >= $testStart && $end <= $testEnd) {
				$nonIdentLarger++;
				$matchtype = "Non_Identical_Larger";
			}
			elsif ($start <= $testStart && $end >= $testEnd) {
				$nonIdentSmaller++;
				$matchtype = "Non_Identical_Smaller";
			}
			else {
				$nonIdentical++;
				$matchtype = "Non_Identical";
			}
			
			# store information about the match for the dumpfile
			my @caseline = split(/\s+/, $_[1]);
			my @testline = split(/\s+/, $hash{$key_sample}{$state}{$chr}{$key_position});
			my $dumpline = "$caseline[0]\t" . ($end - $start) . 
				"\t$caseline[3]\t$caseline[4]\t" .
				"$matchtype\t$_[2]\t$testline[0]\t" . ($testEnd - $testStart) . 
				"\t$testline[3]\t$testline[4]\t$testline[7]\n";
			push (@matchesForDumpfile, $dumpline);
			
		} # [end] for key_position
	} # [end] for key_sample
	
	my @identicalMatchIDs = keys %identicalMatchIDs;
	return ($identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical, 
				\@identicalMatchIDs, \@matchesForDumpfile);
}

# ----------------------------------------------------------------------------
# @parsed_values = &ParseLine( $line )
# ----------------------------------------------------------------------------
sub ParseLine 
{
	# split the line into fields by whitespace
	my @fields = split (/\s+/, $_[0]);
	my ($sample, $state, $chr, $pos) = ('', '', '', '');
	
	# the first field should have the structure ChrX:Start-Stop
	if ($fields[0]) { 
		($chr,$pos) = split(/[:]/, $fields[0]); 
		$pos = '' if !$pos;
	}
	# fourth field should have cn=#
	if ($fields[3]) { 
		$state = ($fields[3] =~ m/(cn=\d)/i) ? $1 : "";
	}
	# fifth field should be the sample ID
	if ($fields[4]) { 
		$sample = $fields[4];
	}
		
	return ($sample, $state, $chr, $pos);
}

#------------------------------------------------------------
sub usage( $ ){
	#the command pod2text calls the commando given between the =head and =cut 
	print "$_[0]\n";
	system("pod2text $0"); 
	exit(1);
}

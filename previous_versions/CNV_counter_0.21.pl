#!/usr/bin/perl

=head1 USAGE CNV_counter.pl
================================================================================

 Required:
 -case <path>       Case cnv file, white space delimited

 Optional:
 -control <path>    Control cnv file, white space delimited
 -pc <1-100>        Percent for match to be considered "identical" (defaults to 90)>
 -out <path>        Output file path (defaults to Casefile.matches.txt)>
 -dump <path>       Dump file path (defaults to Output.matches_dump.txt)
 -freq <0-1>        Match frequency filter. If frequency of matches is greater than this
                    number, suppress listing the matches in output & dump file
                    (defaults to 1 - no filtering)
 -freqcase <0-1>    Only apply filter to matches from the Case file
 -freqcontrol <0-1> Only apply filter to matches from the Control file


 See the comments at the top of this file for the format of
 input and output files and descriptions of match types.
=cut

#
# This script analyses copy number variation (CNV) files to identify rare CNV's.
# By David A Fasel daf2139<at>columbia.edu 11/2012
# Adapted from code by Roel Sterken 6/2009
#
# Basic usage:
# ============
# perl CNV_counter.pl -case CaseFile.txt
#
# Usage with optional parameters:
# ===============================
# perl CNV_counter.pl -case CaseFile.txt -control ControlFile.txt -out OutputPath -dump DumpPath -pc 85 -freq .1
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
#  defaults to 90 (see notes below for types of matches).
#  Since match percentages are not rounded, if you set <pc> to 70, matches
#  of 0.69999999 will not be considered exact matches.  If you want this behavior, set <pc> to 69.5 (for example).
#
# -freq <float>
#  Default is 100 (no filtering).
#  Use this filter to prevent bloated output and dump files.
#  If the match frequency is greater than this value for a particular case:
#    -Suppress listing Id's of identical matches in the output file (see "Format of Output File" below).
#    -Suppress listing each match for the case in the dump file.
#
# -freqcase
#  Same as above, but only apply to matches from the cases file
#
# -freqcontrol
#  Same as above, but only apply to matches from the controls file
#
# Format of Input File(s):
# ========================
# The input file(s) must be whitespace delimited and in this format:
# <chromosome:start-end> <# of snps> <# of bp's> <copy number type> <sample id> <first SNP> <last SNP> <confidence score>
# Example:
# chr16:32438272-32505195 91 66924 state2,cn=1 7 cnvi0020587 cnvi0020628 86.717
#
# Notes:
# -- "number of snps"; "length of range"; "starting SNP"; "ending SNP"
#     are not used but must contain a value
#
# -- "copy number type" must contain <cn=#>, any other information in this field is ignored
#
# -- "sample id" identifies the individual and can be any string, i.e. "7" or "person_1".
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
# * Controls Identical
# * Controls Non-match
# * Controls Frequency
# * Cases Non-identical Smaller: the case CNV is contained within the comparison CNV
# * Cases Non-identical Larger: the comparison CNV is contained within the case CNV
# * Cases Non-identical: none of the matches above apply, however, there is at least some overlap between the CNV's
# * Controls Non-identical Smaller
# * Controls Non-identical Larger
# * Controls Non-identical
# * Case identical matches: comma separated list of sample ID's that contain at least 1 identical match
# * Control identical matches: see above
#
# Format of Dump File:
# ====================
# Case CNV Position
# Case Number of Base Pairs in CNV
# Case CNV State
# Case Sample ID
# Match Type (defined below)
# Matched CNV: File (Case or Control input file)
# Matched CNV: Position
# Matched CNV: Number of Base Pairs in CNV
# Matched CNV: CNV State
# Matched CNV: Sample ID
# Matched CNV: Confidence Score
#
# Match Types:
# ============
#
# "Identical" - matches within the <pc> percent specified.  Checks for a two-way match
#  where both the CNV and the tested CNV are both at least <pc> matches of each other
#
#     |----------------|         [test start/end]
#     TS               TE
#      |---------------|         [case start/end]
#      S                E
#
# "Non-identical larger" - the case CNV is contained within the comparison CNV
#
#     |----------------|         [test start/end]
#     TS               TE
#         |-------|              [case start/end]
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
#  or they are directly adjacent
#
#     |----------|          [test start/end]
#     TS         TE
#           |---------|     [case start/end]
#
# Fisher Test
# ===========
# Note: you must have the Fisher perl module installed in an included Perl directory.
# Module:  Text::NSP::Measures::2D::Fisher::right
# The Fishers exact test per CNV is a right-sided test for the contingency table
#
#          Case  Control    | Totals
# Present   a     b         | a+b
# Absent    x     y         | x+y
# __________________________|________
# Totals    a+x   b+y       | a+b+x+y

use strict;
use warnings;
use IO::Compress::Zip qw(zip $ZipError) ;
use Text::NSP::Measures::2D::Fisher::right;
use Set::IntervalTree;
use feature 'say';

#Check if there are arguments given, otherwise print usage-message
&usage("\n\tError -> please submit all parameters \n") if(scalar(@ARGV) < 1);

# check that the case file path has been specified
my %params = @ARGV;
$params{"-case"} or die "Parameter '-case' is missing. Please specify a case.\n";
my $Counter = 0;
my $start_time;
# get or create defaults for the rest of the params
my $casePath = $params{"-case"};
my $controlPath = $params{"-control"};
my $outPath = $params{"-out"} || "$casePath.matches";
my $dumpPath = $params{"-dump"} || $outPath . "_dump";
my $PERCENT = $params{"-pc"} || 90;
my $FREQ_FILTER = defined $params{"-freq"} ? $params{"-freq"} : 1;
my $FREQ_FILTER_CASE = defined $params{"-freqcase"} ? $params{"-freqcase"} : $FREQ_FILTER;
my $FREQ_FILTER_CTRL = defined $params{"-freqcontrol"} ? $params{"-freqcontrol"} : $FREQ_FILTER;

&validateNumber($PERCENT, 1, 100) or die "Match percent must be between 1 and 100";
&validateNumber($FREQ_FILTER, 0, 1) or die "Frequency filter must be between 0 and 1";
&validateNumber($FREQ_FILTER_CASE, 0, 1) or die "Case frequency filter must be between 0 and 1";
&validateNumber($FREQ_FILTER_CTRL, 0, 1) or die "Control frequency filter must be between 0 and 1";

$outPath .= ".txt";
# check if output files exist in the current directory
if(-e $outPath){
    say "\nWarning: output files already exist";
    print "Overwrite (y/n)?:";
    exit if (<STDIN> ne "y\n");
}

say "\nMatch percent: $PERCENT%  ";
say "(Matches greater than $PERCENT% will be considered 'Identical')\n";
printf "Case frequency filter: %g", $FREQ_FILTER_CASE;
print " (no filter)" if ($FREQ_FILTER_CASE == 100);
say "";
printf "Control max frequency filter: %g", $FREQ_FILTER_CTRL;
print " (no filter)" if $FREQ_FILTER_CTRL == 1;
say "\n(To suppress listing matches in ID list and dump file for \n".
    "common CNV's, set frequency filters to less than 1.)  \n";

$PERCENT = ($PERCENT)/100;

$start_time = time();
my ($sec, $min, $hour) = localtime();
printf "Started: %d:%02d:%02d.\n", $hour, $min, $sec;

open (CASE, $casePath) || die "Couldn't open $casePath\n";
open (CONTROL, $controlPath)|| die "Couldn't open $controlPath\n";
my @casesFile = <CASE>;
my @controlsFile = <CONTROL>;
close CASE;
close CONTROL;
chomp (@casesFile);
chomp (@controlsFile);

# Make a hash for all the cases by CNV State => Chromosome => Interval Tree
my ($cases, $caselines, $totalCases) = &HashCNV("Case", \@casesFile);
my ($controls, $ctrllines, $totalCont) = &HashCNV("Control", \@controlsFile);

say "Found $caselines valid case lines. ";
say "Found $ctrllines valid control lines.";

open (my $OUTPUTFILE, "> $outPath") or
    die "Couldn't open output file $outPath. Check folder permissions.\n";
#open (DUMPFILE,"> $dumpPath") or
#   die "Couldn't open dump file $dumpPath. Check folder permissions.\n";
my $DUMPFILE = new IO::Compress::Zip "$dumpPath.zip", Name => "$dumpPath.txt"
    or die "zip failed: $ZipError\n";

# Headers for the output files
my @headers = qw(
    Position  Num_SNPs  BP_range  CNV_State  Sample_ID
    Start_SNP  End_SNP  Conf_Score
    Pos_Case  Neg_Case  Freq_Case
    Pos_Ctrl Neg_Ctrl  Freq_Ctrl  Fisher_p
    NI_Bigger_Cases  NI_Smaller_Cases NI_Cases CoverageMap_Cases Coverage_Cases
    NI_Bigger_Ctrls  NI_Smaller_Ctrls  NI_Ctrls CoverageMap_Ctrls Coverage_Ctrls
    ID_Matches_Cases  ID_Matches_Controls
);
my @dumpheaders = qw(
    Position  Size CNV_State  Sample_ID  Match_Type  Match_File
    Match_Pos  Match_Size  Match_State  Match_SampleID  Conf_Score
);
# create headers for any additional columns from the input file
my $columns = scalar(split(/\s+/, $casesFile[0]));
for(my $i = 9; $i <= $columns; $i++) {
    splice(@headers, $i-1, 0, "Col_$i");
}
say $OUTPUTFILE join("\t", @headers);
say $DUMPFILE join("\t", @dumpheaders);


for my $caseline (@casesFile)
{
    my ($identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical, $coverageMap, $coverage, $identicalMatchIDs, $dumplines) =
        &findMatches($cases, $caseline, "Case");

    #my $CaseMatches = &findMatches($cases, $caseline, "Case")

    # if the line was bad, &findMatches returns "badline", so go to next line
    next if ($identical eq "badline");

    # if control path is not specified, the case file is also used as the control.
    my ($identicalCont, $nonIdentSmallerCont, $nonIdentLargerCont, $nonIdenticalCont, $coverageMapCont, $coverageCont, $identicalMatchIDsCont, $dumplinesCont) =
        &findMatches($controls, $caseline, "Control");

    # frequency is determined by number of unique samples that have
    # at least one match, divided by total samples
    my $caseFreq = scalar @$identicalMatchIDs / $totalCases;
    my $controlFreq = scalar @$identicalMatchIDsCont / $totalCont;

    # nonMatches are samples that do not contain at least one match
    my $nonMatch = $totalCases - scalar @$identicalMatchIDs;
    my $nonMatchCont = $totalCont - scalar @$identicalMatchIDsCont;

    # Calculate the fisher exact test
    #            case     ctrl
    #  present    n11      n12 | n1p
    #  absent     n21      n22 | n2p
    #            --------------
    #             np1      np2   npp
    my $n11 = (scalar @$identicalMatchIDs);                                    #PresCase;
    my $n1p = (scalar @$identicalMatchIDs) + (scalar @$identicalMatchIDsCont); #PresCase+PresCont;
    my $np1 = $totalCases;                                                     #PresCase+AbsCase;
    my $npp = $totalCases + $totalCont;                                        #PresCase+AbsCase+PresCont+AbsCont;
    #Text::NSP::Measures::2D::Fisher::right::calculateStatistic
    my $rightFisher = calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
    if( (my $errorCode = getErrorCode())){
        say STDERR "\n".$errorCode." - ".Text::NSP::Measures::getErrorMessage();
    }

    # prepare the list of identical matches
    if (!$caseFreq) {@$identicalMatchIDs = ("none") }
    elsif ($caseFreq >= $FREQ_FILTER_CASE){@$identicalMatchIDs = ("exceeds_freq_filter")}
    if (!$controlFreq) {@$identicalMatchIDsCont = ("none") }
    elsif ($controlFreq >= $FREQ_FILTER_CTRL) {@$identicalMatchIDsCont = ("exceeds_freq_filter")}

    # subtract one from each interval in the Coverage Map, since self is counted in cases
    my @arr = map {$_ - 1} split('\.', $coverageMap);
    $coverageMap = join('.', @arr);

    # clean up the original case line by removing any spaces around tabs
    $caseline =~ s:\s*\t\s*:\t:g;
    # build the output line and print it to the output file
    my @output = (
        $caseline, $identical, $nonMatch, $caseFreq,
        $identicalCont, $nonMatchCont, $controlFreq, $rightFisher,
        $nonIdentLarger, $nonIdentSmaller, $nonIdentical, $coverageMap, $coverage,
        $nonIdentLargerCont, $nonIdentSmallerCont, $nonIdenticalCont, $coverageMapCont, $coverageCont,
        join(",", sort @$identicalMatchIDs),
        join(",", sort @$identicalMatchIDsCont)
    );
    say $OUTPUTFILE join("\t", @output);

    # output the dump file if frequency is less than the filter set by user
    print $DUMPFILE @$dumplines if $caseFreq <= $FREQ_FILTER;
    print $DUMPFILE @$dumplinesCont if ($controlFreq <= $FREQ_FILTER);

    # update the status in the terminal and flush the cache to force terminal to display text
    $Counter++;
    ($hour, $min, $sec) = &EllapsedTime($start_time);
    printf "\rProcessed $Counter of $caselines case lines. Time ellapsed: %d:%02d:%02d. ",
        $hour, $min, $sec;
    $| = 1;
}

say "\nDone.\n";
close $OUTPUTFILE;
close $DUMPFILE;
exit;

# ----------------------------------------------------------------------------
# ($hash_ref, line_count) = &HashCNV( $file_type, \@lines )
# ----------------------------------------------------------------------------
sub HashCNV
{
    my ($file_type, @lines) = ($_[0], @{$_[1]}) ;
    my ($line_number, $valid_lines) = (-1, 0);
    my (%hash, %sampleIDs);
    chomp @lines;

    foreach my $CNVline (@lines)
    {
        $line_number++;

        #skip lines that only consist of whitespace characters;
        next if ($CNVline =~ /^\s*$/);

        # Parse the line
        my ($sample, $state, $chr, $pos, $conf) = ParseLine($CNVline);
        my ($start, $end) = split(/-/, $pos);

        # if it's valid, put it in the hash of trees and count it
        if ($sample && $chr && $pos && $state){
            #$hash{$sample}{$state}{$chr}{$pos}=$CNVline;

            if (!$hash{$state}{$chr}) {
              my $tree = Set::IntervalTree->new;
              $hash{$state}{$chr} = $tree;
            }
            my $tree = $hash{$state}{$chr};
            my $line_ref = {
                'line' => $line_number,
                'sample' => $sample,
                'start' => $start,
                'end' => $end,
                'conf' => $conf,
            };
            # insert the node into the tree.  we add 1 to the end so that we will
            # count CNV's that are directly adjacent (but technically don't overlap)
            $tree->insert($line_ref, $start - 1, $end + 1);

            #count the unique sample ID's
            $sampleIDs{$sample} = 1;

            $valid_lines++;
        }
        # print a warning that the line is bad (unless it's the header )
        elsif ($line_number != 0) {
            print "Warning:\n  $file_type file line " . ($line_number + 1) .
              " is missing data and was skipped.\n".
              "  None of these fields can be 0 or empty:\n".
              "  Sample:$sample, State:$state, Chr:$chr, Pos:$pos\n";
        }
    }
    return (\%hash, $valid_lines, scalar(keys %sampleIDs));
}

# ----------------------------------------------------------------------------
# @values = &findMatches( \%hash, $line, $file_type )
# ----------------------------------------------------------------------------
sub findMatches
{
    my ($start, $end, $testStart, $testEnd);
    my (@matchesForDumpfile, %identicalMatchIDs, $matchtype);
    my ($identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical, $coverageMap, $coverage) =
        (0, 0, 0, 0, 0, 0);

    my $hash_ref = $_[0];
    my $file_type = $_[2];
    my ($sample, $state, $chr, $position) = ParseLine($_[1]);
    if (!($sample && $state && $chr && $position)){
        return "badline";
    }
    ($start, $end) = split(/-/, $position);

    my $tree = $hash_ref->{$state}{$chr};
    $coverageMap = &getCoverageMap($start, $end, $tree);

    my $matching_nodes;
    if ($tree) {
        $matching_nodes = $tree->fetch($start, $end);
        $coverage = &getCoverage($sample, $start, $end, $matching_nodes);
    }

    for my $node (@$matching_nodes) {
        my ($testSampleID, $testStart, $testEnd, $testConf) =
            ($node->{'sample'}, $node->{'start'}, $node->{'end'}, $node->{'conf'});

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
        # check for an 'exact' match
        if ($match_case >= $PERCENT && $match_test >= $PERCENT){
            $identical++;
            $matchtype = "Identical";
            $identicalMatchIDs{$testSampleID} = 1;
        }
        # If not an 'exact' match, check for the three types of matches below:
        # For definitions of these matches, see comments at the top of this file.
        elsif ($match_case >= $PERCENT) {
            $nonIdentLarger++;
            $matchtype = "Non_Identical_Larger";
        }
        #elsif (($start <= $testStart && $end >= $testEnd) || ($match_test < $PERCENT && $match_test > 0 )) {
        elsif ($match_test >= $PERCENT) {
            $nonIdentSmaller++;
            $matchtype = "Non_Identical_Smaller";
        }
        else {
            $nonIdentical++;
            $matchtype = "Non_Identical";
            #print "$chr-$testStart:$testEnd $_[2]\n";
        }

        # store information about the match for the dumpfile
        my $dumpline = "$chr:$start-$end\t".($end - $start)."\t".
            "$state\t$sample\t$matchtype\t$file_type\t".
            "$chr:$testStart-$testEnd\t".($testEnd - $testStart)."\t".
            "$state\t$testSampleID\t$testConf\n";

        push (@matchesForDumpfile, $dumpline);
    }

    my @identicalMatchIDs = keys %identicalMatchIDs;

    my %hash =
    return
    return ($identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical,
        $coverageMap, $coverage, \@identicalMatchIDs, \@matchesForDumpfile);

}

# ----------------------------------------------------------------------------
# @parsed_values = &ParseLine( $line )
# ----------------------------------------------------------------------------
sub ParseLine
{
    # split the line into fields by whitespace
    my @fields = split (/\s+/, $_[0]);
    my ($sample, $state, $chr, $pos, $conf) = ('', '', '', '', '');

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
    if ($fields[7]) {
        $conf = $fields[7];
    }


    return ($sample, $state, $chr, $pos, $conf);
}

sub EllapsedTime()
{
    my $seconds_ellapsed = time() - shift;
    my $hour = int($seconds_ellapsed / 3600);
    my $min = int(($seconds_ellapsed % 3600) / 60);
    my $sec = ($seconds_ellapsed % 3600) % 60;
    return $hour, $min, $sec;
}

# ----------------------------------------------------------------------------
# int = getCoverageMap($start, $end, $tree);
# ----------------------------------------------------------------------------
sub getCoverageMap( $ ){
    my $INTERVALS = 10;
    my ($start, $end, $tree) = @_;
    my $length = ($end - $start) / $INTERVALS;
    my @counts;

    if (!$tree) {
        return join('.', (0) x $INTERVALS)
    }

    $start += $length / 2;
    for (my $i = 0; $i < $INTERVALS; $i++) {

        my $nodes = $tree->fetch($start, $start + 1);

        push(@counts, scalar @$nodes);
        $start += $length;
    }

    return join('.', @counts);
}

# ----------------------------------------------------------------------------
# float = getCoverage($sample, $start, $end, $matching_nodes)
# ----------------------------------------------------------------------------
sub getCoverage( $ ) {
    my ($sample, $start, $end, $matching_nodes) = @_;
    my $coverage_tree = Set::IntervalTree->new;
    my $nodeid = 1;

    for my $node (@$matching_nodes) {
        my ($nodeStart, $nodeEnd) = ($node->{'start'}, $node->{'end'});

        # Skip the node if the sampleID and position match (don't count self)
        next if ($sample eq $node->{'sample'} && $start == $nodeStart && $end == $nodeEnd);

        # if any node completely overlaps, then coverage is 1
        if($nodeStart <= $start && $nodeEnd >= $end){
            return 1;
        }

        #trim the node if it extends past the ends of the case CNV
        $nodeStart = $start if $nodeStart < $start;
        $nodeEnd = $end if $nodeEnd > $end;


        # if after trimming the match cnv, it only overlaps the case CNV by one bp,
        # add one so tree->fetch() doesn't throw an error
        if ($nodeEnd == $nodeStart) {$nodeEnd++}
        # check for any overlapping nodes
        my $coverage_nodes = $coverage_tree->fetch($nodeStart, $nodeEnd);

        #if no overlapping nodes, then add to the tree
        if (scalar @$coverage_nodes == 0) {
            $coverage_tree->insert(
                {id => $nodeid, start => $nodeStart, end => $nodeEnd},
                $nodeStart, $nodeEnd
            );
        }
        # If adding the node into our coverage tree would create overlapping nodes,
        # remove the overlapping nodes and consolidate them into one node.
        # Then add back in as a single consolidated node
        else {
            my $overlapping_nodes = $coverage_tree->remove($nodeStart, $nodeEnd);
            my ($newStart, $newEnd);
            for my $node (@$overlapping_nodes) {
                $newStart = $node->{'start'} if (!$newStart || $node->{'start'} < $newStart);
                $newEnd = $node->{'end'} if (!$newEnd || $node->{'end'} > $newEnd);
            }
            # if the new consolidated node overlaps the entire CNV, then coverage is 1
            if($nodeStart <= $start && $nodeEnd >= $end){
                return 1;
            }
            $coverage_tree->insert({id => $nodeid, start => $newStart, end => $newEnd},
                $newStart, $newEnd);
        }
        $nodeid++;
    }

    #sum the coverage for each node in the tree
    my $coverage_nodes = $coverage_tree->fetch($start, $end);
    my $sum_coverage = 0;
    for my $node (@$coverage_nodes) {
        $sum_coverage += ($node->{'end'} - $node->{'start'} );
    }
    #divide coverage by the length of the CNV to get coverage rate
    return sprintf ("%.2g", $sum_coverage / ($end - $start));

}

# ----------------------------------------------------------------------------
# bool = validateNumber($num, $low, $high)
# ----------------------------------------------------------------------------
sub validateNumber( $ ){
    my ($num, $low, $high) = @_;
    return 1 if ($num >= $low && $num <= $high );
    return 0;
}
#------------------------------------------------------------
sub usage( $ ){
    #the command pod2text calls the commando given between the =head and =cut
    print "$_[0]\n";
    system("pod2text $0");
    exit(1);
}

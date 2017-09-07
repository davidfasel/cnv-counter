#!/usr/bin/perl

=head1 USAGE CNV Counter
================================================================================

 Required:
 -c --case <path>     Case CNV file, white space delimited

 Optional:
 -n --control [file]    Control cnv file, white space delimited
 -m --match [1-100]     Percent for match to be considered "identical" (defaults to 90)
 -o --out [file]        Output file path (defaults to Casefile.matches.txt)
 -d --dump [file]       Dump file path (defaults to Output.matches_dump.zip)
 -f --freq [0-1]        Frequency filter used to prevent bloated output and dump files.
                        If the match frequency is greater than this value for a case CNV:
                          -Suppress listing ID's of identical matches in the last columns
                           of the output file.
                          -Suppress listing matches for the case CNV in the dump file.
 -a --freqcase [0-1]    Frequency filter applied to matches from the Case file only
 -b --freqcontrol [0-1] Frequency filter applied to matches from the Control file only
 -g --greedy            Include 'Non-Identical Larger' match types when counting
                        positive cases and controls.

 Basic usage:
 ============
 CNV_counter.pl --case CaseFile.txt

 Usage with optional parameters:
 ===============================
 CNV_counter.pl -c CaseFile.txt -n ControlFile.txt -o OutputPath -d DumpPath -m 85 -f 0.1

 See the comments at the top of this file for the format of
 input and output files and descriptions of match types.

=cut

# This script analyses copy number variation (CNV) files to identify rare CNV's.
# By David A Fasel daf2139<at>columbia.edu 11/2012
# Adapted from code by Roel Sterken 6/2009
#
# Version info:
# 0.20:
#  - optimized with "interval tree" data type
#  - zip dump file
# 0.21:
#  - fixed bug in Fisher p-val (because positive controls were always >= 1)
#  - fixed bug in Fisher p-val (Frequency filter changed positive case/control values)
#  - if extra columns in input file, put placeholder column headers in output file
#  - added Coverage Map and Coverage Freq
#  - option to set independent frequency filters for case and controls
#  - output file is tab delimited only, even if input file uses mix of tabs and spaces
#  - fixed bug some 'Non-Identical' Matches being classified as 'Non-Identical Smaller'
#  - broke dumpfile (doesn't output lines)
# 0.22
#  - fixed dumpfile from 0.21
#  - calculate reciprocal identical matches (deletions cn0,1 vs duplications cn3+)
#  - Pos case/ctrl and "identical" case/ctrl are now 2 different columns (since they
#    can be different values if <match> percent is less than 50)
#  - Added -greedy flag to count NI_larger as positive case/control
#
# todo:
# copy case info to controls if no controls file
# suppress self matching in dump file?
# count reciprocal matches for chrX
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
# -- "number of snps"; "# of bp's"; "starting SNP"; "ending SNP"
#     are not used but must contain a value
#
# -- "copy number type" must contain <cn=#>, any other information in this field is ignored
#
# -- "sample id" identifies the individual and can be any string, i.e. "7" or "person_1".
#     It cannot be blank or 0.
#
# Format of Output File (tab delimited):
# ======================
# Each line in the output file corresponds to a line in the Case file:
#
# * Case line:
#     The entire original line from the case file.  The CNV in this line is the
#     'reference' CNV.
# * Pos_Case (Identical case CNV's):
#     Number of individuals in the Case file that have at least one 'identical' CNV.
# * Neg_Case (Cases Non-match):
#     Number of individuals in the Case file that do not have at least one 'identical' CNV.
# * Freq_Case:
#     Frequency of cases (individuals in the case file with at least one 'identical' CNV
#     divided by total individuals in the Case file)
# * Pos_Ctrl:
#     Number of CNV's in the Control file that are 'identical' to the reference CNV.
# * Neg_Ctrl:
#     Number of individuals in the Control file that do not have at least one 'identical' CNV.
# * Freq_Ctrl:
#     Frequency of controls (individuals in the Control file with at least one
#     'identical' CNV divided by total individuals in the Controls file)
# * Fisher_p:
#     P-value calculated using Fisher
# * Ident_Case:
#     Number of CNV's in the Case file that are 'identical' to the reference CNV.
#     There could be more than one match per individual if <match> percent is set to
#     less than 50.
# * Reciprocal_Case:
#     Number of identical matches in reciprocal CNV types.  For example, if the
#     reference CNV is a deletion (copy number of 0 or 1), this is the count of
#     matching CNV's that are duplications (copy number > 2)
# * NI_Larger_Cases (Non-identical Bigger):
#     Number of case CNV's that contain <match> percent of the reference CNV
# * NI_Smaller_Cases (Non-identical Smaller):
#     Number of CNV's in the Case file that are contained within the reference CNV
# * NI_Cases (Non-identical):
#     Count of CNV's that do not meet any of the above matching criteria,
#     however, there is at least some overlap with the reference CNV
# * CoverageMap_Case:
#     Shows the number of matching CNV's at individual locations within the CNV
# * Coverage_Case:
#     The fraction of the CNV that is covered by any case CNVs
#
# * Control matches:
#     The previous 7 fields are repeated for matches found in the Control file.
#     ... etc.
#
# * ID_Matches_Cases:
#     Comma separated list of individual's from the Case file that contain at least one
#    'identical' match
# * ID_Matches_Ctrl:
#     Comma separated list of individual's from the Control file that contain at least one
#    'identical' match
#
# Format of Dump File:
# ====================
# Case CNV Position
# Case Number of Base Pairs in CNV
# Case CNV State
# Case Sample ID
# Match Type (defined below)
# Matched CNV - File (Case or Control input file)
# Matched CNV - Position
# Matched CNV - Number of Base Pairs in CNV
# Matched CNV - CNV State
# Matched CNV - Sample ID
# Matched CNV - Confidence Score
#
# Match Types:
# ============
#
# "Identical" - matches within the <match> percent specified.  Checks for a two-way match
#  where both the reference CNV and the tested CNV are both at least <match> percent
#  matches of each other
#
#     |----------------|         [test start/end]
#     TS               TE
#      |---------------|         [ref start/end]
#      S                E
#
# "Non-identical larger" - the reference CNV is contained within the test CNV
#
#     |----------------|         [test start/end]
#     TS               TE
#         |-------|              [ref start/end]
#         S       E
#
# "Non-identical smaller" - the test CNV is contained within the reference CNV
#
#     |----------|          [test start/end]
#     TS         TE
#   |-----------------|     [ref start/end]
#   S                 E
#
# "Non-identical" - none of the matches above apply,
#  however, there is at least some overlap between the CNV's
#
#     |----------|          [test start/end]
#     TS         TE
#           |---------|     [ref start/end]
#
# Fisher Test
# ===========
# Module:  Text::NSP::Measures::2D::Fisher::right
# The Fishers exact test per CNV is a right-sided test for the contingency table
#
#          Case  Control    | Totals
# Present   a     b         | a+b
# Absent    x     y         | x+y
# __________________________|________
# Totals    a+x   b+y       | a+b+x+y
#

use strict;
use warnings;
use IO::Compress::Zip qw(zip $ZipError) ;
use Text::NSP::Measures::2D::Fisher::right;
use Set::IntervalTree;
use Getopt::Long;
use feature 'say';

my $Counter = 0;
my $start_time;
my ($casePath, $controlPath, $outPath, $dumpPath);
my $FREQ_FILTER = 1;
my ($FREQ_FILTER_CASE, $FREQ_FILTER_CTRL);
my $MATCH_RATE =  90;
my $MAP_INTERVAL = 10;
my $GREEDY = 0; #false

# set options from paramaters
GetOptions (
    'case|c=s' => \$casePath,
    'control|n=s' => \$controlPath,
    'out|o=s' => \$outPath,
    'match|m=i' => \$MATCH_RATE,
    'freq|f=s' => \$FREQ_FILTER,
    'freqcase|a=s' => \$FREQ_FILTER_CASE,
    'freqcontrol|b=s' => \$FREQ_FILTER_CTRL,
    'greedy|g' => \$GREEDY,
);
# check that the case file path has been specified
$casePath or die &usage(
    "\n**Error: Please specify a path to the case file using --case or -c.\n");

$FREQ_FILTER_CASE = $FREQ_FILTER if not defined $FREQ_FILTER_CASE;
$FREQ_FILTER_CTRL = $FREQ_FILTER if not defined $FREQ_FILTER_CTRL;
&validateNumber($MATCH_RATE, 1, 100) or die "Match percent must be between 1 and 100";
&validateNumber($FREQ_FILTER, 0, 1) or die "Frequency filter must be between 0 and 1";
&validateNumber($FREQ_FILTER_CASE, 0, 1) or die "Case frequency filter must be between 0 and 1";
&validateNumber($FREQ_FILTER_CTRL, 0, 1) or die "Control frequency filter must be between 0 and 1";

$outPath = $outPath || "$casePath.matches";
$dumpPath = $outPath . "_dump";
$outPath .= ".txt";
# check if output files exist in the current directory
if(-e $outPath){
    say "\nWarning: output files already exist";
    print "Overwrite (y/n)?:";
    exit if (<STDIN> ne "y\n");
}

say "\nMatch percent: $MATCH_RATE% (Matches greater than $MATCH_RATE% will be " .
    "considered 'Identical')\n";
say "Max Frequency filters:";
printf "Case: %g", $FREQ_FILTER_CASE;
print " (no filter)" if ($FREQ_FILTER_CASE == 1);
say "";
printf "Control: %g", $FREQ_FILTER_CTRL;
print " (no filter)" if $FREQ_FILTER_CTRL == 1;
say "\n";

say "Greedy is on. ('Non-Identical Larger' matches will be counted" .
    " as positive cases/controls.) \n" if $GREEDY;

$MATCH_RATE = ($MATCH_RATE)/100;

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

open (my $OUTPUTFILE, "> $outPath") or
    die "Couldn't open output file $outPath. Check folder permissions.\n";
my $DUMPFILE = new IO::Compress::Zip "$dumpPath.zip", Name => "$dumpPath.txt"
    or die "zip failed: $ZipError\n";

# Headers for the output files
my @headers = qw(
    Position  Num_SNPs  BP_range  CNV_State  Sample_ID
    Start_SNP  End_SNP  Conf_Score
    Pos_Case  Neg_Case  Freq_Case
    Pos_Ctrl  Neg_Ctrl  Freq_Ctrl
    Fisher_p
    Ident_Cases  Recip_Cases  NI_Larger_Cases  NI_Smaller_Cases  NI_Cases
    CoverageMap_Cases  Coverage_Cases
    Ident_Ctrls  Recip_Ctrls  NI_Larger_Ctrls  NI_Smaller_Ctrls  NI_Ctrls
    CoverageMap_Ctrls  Coverage_Ctrls
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

# Make a hash for all the cases by CNV State => Chromosome => Interval Tree
my ($cases_hash, $caselines, $totalCases) = &HashCNV("Case", \@casesFile);
my ($controls_hash, $ctrllines, $totalCont) = &HashCNV("Control", \@controlsFile);

say "Found $caselines valid case lines. ";
say "Found $ctrllines valid control lines.";

#while (my ($i, $el) = each @Array) {
for my $caseline (@casesFile)
{
    my ($sample, $state, $chr, $start, $end) = ParseLine($caseline);
    if (!($sample && $state && $chr && $start && $end)){
        next;
    }

    # get the case and control matches
    my $caseMatches = &findMatches($cases_hash, $caseline, "Case");
    my $ctrlMatches = &findMatches($controls_hash, $caseline, "Control");

    # Count the cases and control samples
    my $identicalMatchIDs = $caseMatches->{'identicalMatchIDs'};
    my $identicalMatchIDsCont = $ctrlMatches->{'identicalMatchIDs'};
    my $cases = scalar @$identicalMatchIDs;
    my $controls = scalar @$identicalMatchIDsCont;
    if ($GREEDY) {
        $cases += $caseMatches->{nonIdentLarger};
        $controls += $ctrlMatches->{nonIdentLarger};
    }

    # frequency is determined by number of unique samples that have
    # at least one match, divided by total samples
    my $caseFreq = $cases / $totalCases;
    my $controlFreq = $controls / $totalCont;

    # nonMatches are samples that do not contain at least one match
    my $nonMatch = $totalCases - $cases;
    my $nonMatchCont = $totalCont - $controls;

    # get the fisher pvalue
    my $rightFisher = &getFisher($cases, $controls, $totalCases, $totalCont);

    # get the coverage rate and the coverage maps
    my $tree = $cases_hash->{$state}{$chr};
    my $coverage = &getCoverage($sample, $start, $end, $tree);
    my $coverageMap = &getCoverageMap($start, $end, $MAP_INTERVAL, $tree);

    $tree = $controls_hash->{$state}{$chr};
    my $coverageCont = &getCoverage($sample, $start, $end, $tree);
    my $coverageMapCont = &getCoverageMap($start, $end, $MAP_INTERVAL, $tree);

    # edit the list of identical matches if empty or clear it if freq filter is exceeded
    if (!$caseFreq)
        {@$identicalMatchIDs = ("none") }
    elsif ($caseFreq >= $FREQ_FILTER_CASE)
        {@$identicalMatchIDs = ("exceeds_freq_filter")}

    if (!$controlFreq)
        {@$identicalMatchIDsCont = ("none") }
    elsif ($controlFreq >= $FREQ_FILTER_CTRL)
         {@$identicalMatchIDsCont = ("exceeds_freq_filter")}

    # clean up the original case line by removing any spaces around tabs
    $caseline =~ s:\s*\t\s*:\t:g;

    # build the output line and print it to the output file
    my @output = (
        $caseline,
        $cases,
        $nonMatch,
        sprintf("%.4g", $caseFreq),
        $controls,
        $nonMatchCont,
        sprintf("%.4g", $controlFreq),
        sprintf("%.4g", $rightFisher),
        $caseMatches->{identical},
        $caseMatches->{recip_identical},
        $caseMatches->{nonIdentLarger},
        $caseMatches->{nonIdentSmaller},
        $caseMatches->{nonIdentical},
        $coverageMap,
        sprintf("%.2g", $coverage),
        $ctrlMatches->{identical},
        $ctrlMatches->{recip_identical},
        $ctrlMatches->{nonIdentLarger},
        $ctrlMatches->{nonIdentSmaller},
        $ctrlMatches->{nonIdentical},
        $coverageMapCont,
        sprintf("%.2g", $coverageCont),
        join(",", sort @$identicalMatchIDs),
        join(",", sort @$identicalMatchIDsCont)
    );
    say $OUTPUTFILE join("\t", @output);

    # output the dump file if frequency is less than the filter set by user
    my $dumplines = $caseMatches->{'dumplines'};
    my $dumplinesCont = $ctrlMatches->{'dumplines'};

    say $DUMPFILE join("\n", @$dumplines) if ($caseFreq <= $FREQ_FILTER);
    say $DUMPFILE join("\n", @$dumplinesCont) if (scalar @$dumplinesCont && $controlFreq <= $FREQ_FILTER);

    # update the status in the terminal and flush cache to force terminal to display text
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
    my ($line_number, $valid_lines) = (0, 0);
    my (%hash, %sampleIDs);
    chomp @lines;

    foreach my $CNVline (@lines)
    {
        $line_number++;

        #skip lines that only consist of whitespace characters;
        next if ($CNVline =~ /^\s*$/);

        # Parse the line
        my ($sample, $state, $chr, $start, $end, $conf) = ParseLine($CNVline);

        # if it's valid, put it in the hash of trees and count it
        if ($sample && $chr && $start && $end && $state){
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
        elsif ($line_number != 1) {
            print "Warning:\n  $file_type file line " . ($line_number) .
              " is missing data and was skipped.\n".
              "  None of these fields can be 0 or empty:\n".
              "  Sample:$sample, State:$state, Chr:$chr, Start:$start, End:$end\n";
        }
    }
    return (\%hash, $valid_lines, scalar(keys %sampleIDs));
}

# ----------------------------------------------------------------------------
# @values = &findMatches( \%hash, $line, $file_type )
# ----------------------------------------------------------------------------
sub findMatches
{
    my (%Matches, @matchesForDumpfile, %identicalMatchIDs, $matchtype);
    my ($identical, $recip_identical, $nonIdentSmaller, $nonIdentLarger, $nonIdentical) = (0, 0, 0, 0, 0);

    my ($hash_ref, $line, $file_type) = @_;
    my ($sample, $state, $chr, $start, $end) = ParseLine($line);

    for my $hash_state (keys %$hash_ref) {
        my $isReciprocal = &isReciprocalState($hash_state, $state, $chr);
        next if not ($isReciprocal || $hash_state eq $state);

        my $tree = $hash_ref->{$hash_state}{$chr};

        my $matching_nodes;
        if ($tree) {
            $matching_nodes = $tree->fetch($start, $end);
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

            # check for reciprocal and identical match
            if ($isReciprocal) {
                if ($match_case >= $MATCH_RATE && $match_test >= $MATCH_RATE) {
                    $recip_identical++;
                }
                # If reciprocal, we only care about identical matches, so skip the rest
                next;
            }
            # check for identical match
            if ($match_case >= $MATCH_RATE && $match_test >= $MATCH_RATE){
                # if we are looking at CNV's with the same state, then
                # count matches normally
                if($hash_state eq $state) {
                    $identical++;
                    $matchtype = "Identical";
                    $identicalMatchIDs{$testSampleID} = 1;
                }
            }
            # If not an 'identical' match, check for the three types of matches below:
            # For definitions of these matches, see comments at the top of this file.
            elsif ($match_case >= $MATCH_RATE) {
                $nonIdentLarger++;
                $matchtype = "Non_Identical_Larger";
            }
            elsif ($match_test >= $MATCH_RATE) {
                $nonIdentSmaller++;
                $matchtype = "Non_Identical_Smaller";
            }
            else {
                $nonIdentical++;
                $matchtype = "Non_Identical";
                #print "$chr-$testStart:$testEnd $_[2]\n";
            }

            # store information about the match for the dumpfile
            my @dumpline = (
                "$chr:$start-$end",
                ($end - $start),
                $state,
                $sample,
                $matchtype,
                $file_type,
                "$chr:$testStart-$testEnd",
                ($testEnd - $testStart),
                $state,
                $testSampleID,
                $testConf
            );

            push (@matchesForDumpfile, join("\t", @dumpline));
        }
    }

    my @identicalMatchIDs = keys %identicalMatchIDs;

    $Matches{identical} = $identical;
    $Matches{recip_identical} = $recip_identical;
    $Matches{nonIdentSmaller} = $nonIdentSmaller;
    $Matches{nonIdentLarger} = $nonIdentLarger;
    $Matches{nonIdentical} = $nonIdentical;
    $Matches{identicalMatchIDs} = \@identicalMatchIDs;
    $Matches{dumplines} = \@matchesForDumpfile;

    return \%Matches;

}

# ----------------------------------------------------------------------------
# @parsed_values = &ParseLine( $line )
# ----------------------------------------------------------------------------
sub ParseLine
{
    # split the line into fields by whitespace
    my @fields = split (/\s+/, $_[0]);
    my ($sample, $state, $chr, $pos, $start, $end, $conf) = ('', '', '', '', '', '', '');

    # the first field should have the structure ChrX:Start-Stop
    if ($fields[0]) {
        ($chr,$pos) = split(/[:]/, $fields[0]);
        if ($pos) {
            ($start, $end) = split(/-/, $pos);
        }
    }
    # fourth field should have cn=#
    if ($fields[3]) {
         $fields[3] =~ m/(cn=\d)/i;
         $state = $1;
    }
    # fifth field should be the sample ID
    if ($fields[4]) {
        $sample = $fields[4];
    }
    if ($fields[7]) {
        $conf = $fields[7];
    }

    return ($sample, $state, $chr, $start, $end, $conf);
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
# string = getCoverageMap($start, $end, $intervals, $tree);
# ----------------------------------------------------------------------------
sub getCoverageMap( $ ){
    my ($start, $end, $intervals, $tree) = @_;
    my $length = ($end - $start) / $intervals;
    my @counts;

    if (!$tree) {
        return join('.', (0) x $intervals)
    }

    $start += $length / 2;
    for (my $i = 0; $i < $intervals; $i++) {

        my $nodes = $tree->fetch($start, $start + 1);

        push(@counts, scalar @$nodes);
        $start += $length;
    }

    return join('.', @counts);
}

# ----------------------------------------------------------------------------
# float = getCoverage($sample, $start, $end, $tree)
# ----------------------------------------------------------------------------
sub getCoverage( $ ) {
    my ($sample, $start, $end, $tree) = @_;
    my $coverage_tree = Set::IntervalTree->new;
    my $nodeid = 1;
    my $matching_nodes;
    if ($tree) {
        $matching_nodes = $tree->fetch($start, $end);
    }
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
    return $sum_coverage / ($end - $start);

}

# ----------------------------------------------------------------------------
# float = getFisher($cases, $controls, $totalCases, $totalCont)
# ----------------------------------------------------------------------------
sub getFisher( $ ) {
    my ($cases, $controls, $totalCases, $totalCont) = @_;
    # Calculate the fisher exact test
    #            case     ctrl
    #  present    n11      n12 | n1p
    #  absent     n21      n22 | n2p
    #            --------------
    #             np1      np2   npp
    my $n11 = $cases;                                           #PresCase;
    my $n1p = $cases + $controls;                               #PresCase+PresCont;
    my $np1 = $totalCases;                                      #PresCase+AbsCase;
    my $npp = $totalCases + $totalCont;                         #PresCase+AbsCase+PresCont+AbsCont;
    #Text::NSP::Measures::2D::Fisher::right::calculateStatistic
    my $rightFisher = calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
    if( my $errorCode = getErrorCode() ){
        say STDERR "\n".$errorCode." - ".Text::NSP::Measures::getErrorMessage();
    }
    return $rightFisher;
}

# ----------------------------------------------------------------------------
# bool = isReciprocalState($hash_state, $state, $chr)
# ----------------------------------------------------------------------------
sub isReciprocalState( $ ){
    my ($hash_state, $state, $chr) = @_;
    $hash_state =~ /(\d+)/;
    $hash_state = $1;
    $state =~ /(\d+)/;
    $state = $1;
    if ( (lc($chr) ne "chrx") && ($state < 2 && $hash_state > 2) || ($state > 2 && $hash_state < 2) ) {
        return 1;
    }
    return 0;
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

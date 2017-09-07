#!/usr/bin/perl
use strict;
use warnings;
use IO::Compress::Zip qw(zip $ZipError) ;
use Text::NSP::Measures::2D::Fisher::right;
use Set::IntervalTree;
use Getopt::Long;  #process command line arguments
use Pod::Usage;
use feature 'say';

my $Counter = 0;
my ($casePath, $controlPath, $outPath, $dumpPath);
my $FREQ_FILTER = 1;
my ($FREQ_FILTER_CASE, $FREQ_FILTER_CTRL);
my $MATCH_RATE =  90;
my $MAP_INTERVAL = 10;
my $ShowHelp = 0; #false
my ($ShowMatchIDCols, $ShowMatchIDColsCase, $ShowMatchIDColsCtrl) = ("", "", "");
my ($PositiveMatches, $PositiveMatchesCase, $PositiveMatchesCtrl) = ("", "", "");
my (%ShowColumns, %PosMatches);
my $start_time = time();
my ($sec, $min, $hour) = localtime();

# set options from paramaters
GetOptions (
    'help|h' => \$ShowHelp,
    'case|c=s' => \$casePath,
    'control|n=s' => \$controlPath,
    'out|o=s' => \$outPath,
    'match|m=i' => \$MATCH_RATE,
    'freq|f=s' => \$FREQ_FILTER,
    'freqcase|a=s' => \$FREQ_FILTER_CASE,
    'freqcontrol|b=s' => \$FREQ_FILTER_CTRL,
    'col=s' => \$ShowMatchIDCols,
    'ccol=s' => \$ShowMatchIDColsCase,
    'ncol=s' => \$ShowMatchIDColsCtrl,
    'pos=s' => \$PositiveMatches,
    'cpos=s' => \$PositiveMatchesCase,
    'npos=s' => \$PositiveMatchesCtrl,
);

$ShowHelp and &help();

# check that the case file path has been specified
$casePath && -f $casePath or  &usage(
    "\n**Error: Please specify a path to the case file using --case or -c.\n");
-f $controlPath or die &usage(
    "\n**Error: Please specify a path to the control file using --control or -n.\n");
$outPath = $outPath || "$casePath.matches";
$dumpPath = $outPath . "_dump";
$outPath .= ".txt";

$FREQ_FILTER_CASE = $FREQ_FILTER if not defined $FREQ_FILTER_CASE;
$FREQ_FILTER_CTRL = $FREQ_FILTER if not defined $FREQ_FILTER_CTRL;
&validateNumber($MATCH_RATE, 1, 100) or &usage("Match percent must be between 1 and 100");
&validateNumber($FREQ_FILTER, 0, 1) or &usage("Frequency filter must be between 0 and 1");
&validateNumber($FREQ_FILTER_CASE, 0, 1) or &usage("Case frequency filter must be between 0 and 1");
&validateNumber($FREQ_FILTER_CTRL, 0, 1) or &usage("Control frequency filter must be between 0 and 1");

$ShowColumns{'NILCase'} = 1 if ($ShowMatchIDCols =~ /l/i || $ShowMatchIDColsCase =~ /l/i);
$ShowColumns{'NISCase'} = 1 if ($ShowMatchIDCols =~ /s/i || $ShowMatchIDColsCase =~ /s/i);
$ShowColumns{'NICase'} = 1  if ($ShowMatchIDCols =~ /n/i || $ShowMatchIDColsCase =~ /n/i);
$ShowColumns{'NILCtrl'} = 1 if ($ShowMatchIDCols =~ /l/i || $ShowMatchIDColsCtrl =~ /l/i);
$ShowColumns{'NISCtrl'} = 1 if ($ShowMatchIDCols =~ /s/i || $ShowMatchIDColsCtrl =~ /s/i);
$ShowColumns{'NICtrl'} = 1  if ($ShowMatchIDCols =~ /n/i || $ShowMatchIDColsCtrl =~ /n/i);

$PosMatches{'NILCase'} = 1 if ($PositiveMatches =~ /l/i || $PositiveMatchesCase =~ /l/i);
$PosMatches{'NISCase'} = 1 if ($PositiveMatches =~ /s/i || $PositiveMatchesCase =~ /s/i);
$PosMatches{'NICase'} = 1  if ($PositiveMatches =~ /n/i || $PositiveMatchesCase =~ /n/i);
$PosMatches{'NILCtrl'} = 1 if ($PositiveMatches =~ /l/i || $PositiveMatchesCtrl =~ /l/i);
$PosMatches{'NISCtrl'} = 1 if ($PositiveMatches =~ /s/i || $PositiveMatchesCtrl =~ /s/i);
$PosMatches{'NICtrl'} = 1  if ($PositiveMatches =~ /n/i || $PositiveMatchesCtrl =~ /n/i);

print  "\nMatch percent: $MATCH_RATE% (Matches greater than $MATCH_RATE% will be " .
       "considered 'Identical')\n\n".
       "Max Frequency filters:\n";
printf "Case: %g", $FREQ_FILTER_CASE;
print  " (no filter)" if ($FREQ_FILTER_CASE == 1);
printf "\nControl: %g", $FREQ_FILTER_CTRL;
print  " (no filter)" if $FREQ_FILTER_CTRL == 1;
printf "\n\nStarted: %d:%02d:%02d.\n", $hour, $min, $sec;
$MATCH_RATE = ($MATCH_RATE)/100;

open (CASE, $casePath) || &usage("Couldn't open $casePath\n");
open (CONTROL, $controlPath)|| &usage("Couldn't open $controlPath\n");
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
);

push @headers, 'ID_Matches_Cases';
push @headers, 'NILarger_ID_Matches_Cases' if $ShowColumns{'NILCase'};
push @headers, 'NISmaller_ID_Matches_Cases' if $ShowColumns{'NISCase'};
push @headers, 'NI_ID_Matches_Cases' if $ShowColumns{'NICase'};

push @headers, 'ID_Matches_Ctrls';
push @headers, 'NILarger_ID_Matches_Ctrls' if $ShowColumns{'NILCtrl'};
push @headers, 'NISmaller_ID_Matches_Ctrls' if $ShowColumns{'NISCtrl'};
push @headers, 'NI_ID_Matches_Ctrls' if $ShowColumns{'NICtrl'};

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

for my $caseline (@casesFile)
{
    my ($caseMatches, $ctrlMatches);
    my ($sample, $state, $chr, $start, $end) = ParseLine($caseline);
    if (!($sample && $state && $chr && $start && $end)){
        next;
    }

    # get the case and control matches
    $caseMatches = &findMatches($cases_hash, $caseline, "Case");
    $ctrlMatches = &findMatches($controls_hash, $caseline, "Control");
    
    # count the positive matches. A 'positive' match is an 'identical' match by 
    # default, but can include other match types if set by the user
    my (%posCases, %posCtrls);
    $posCases{$_} = 1 for @{$caseMatches->{'identMatchIDs'}};
    if ($PosMatches{'NILCase'}) {$posCases{$_} = 1 for @{$caseMatches->{'NILMatchIDs'}}};
    if ($PosMatches{'NISCase'}) {$posCases{$_} = 1 for @{$caseMatches->{'NISMatchIDs'}}} ;
    if ($PosMatches{'NICase'}) {$posCases{$_} = 1 for @{$caseMatches->{'NIMatchIDs'}}};
    
    $posCtrls{$_} = 1 for @{$ctrlMatches->{'identMatchIDs'}};
    if ($PosMatches{'NILCtrl'}) {$posCtrls{$_} = 1 for @{$ctrlMatches->{'NILMatchIDs'}}};
    if ($PosMatches{'NISCtrl'}) {$posCtrls{$_} = 1 for @{$ctrlMatches->{'NISMatchIDs'}}};
    if ($PosMatches{'NICtrl'}) {$posCtrls{$_} = 1 for @{$ctrlMatches->{'NIMatchIDs'}}};
    
    # frequency is determined by number of unique samples that have
    # at least one match, divided by total samples
    my $cases = scalar(keys %posCases);
    my $controls = scalar(keys %posCtrls);
    my $caseFreq = $cases / $totalCases;
    my $ctrlFreq = $controls / $totalCont;

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

    # edit each list of matches if empty or clear it if freq filter is exceeded
    # we have all of these matchIDs: identMatchIDs, NILMatchIDs, etc.
    for my $key (keys %$caseMatches) {
        if ($key =~ /MatchIDs/ && not @{$caseMatches->{$key}}) {$caseMatches->{$key} = ["none"]}
        elsif ($key =~ /MatchIDs/ && $caseFreq > $FREQ_FILTER_CASE) {$caseMatches->{$key} = ["filtered"]}
    }
    for my $key (keys %$ctrlMatches) {
        if ($key =~ /MatchIDs/ && not @{$ctrlMatches->{$key}}) {$ctrlMatches->{$key} = ["none"]}
        elsif ($key =~ /MatchIDs/ && $ctrlFreq > $FREQ_FILTER_CTRL) {$ctrlMatches->{$key} = ["filtered"]}
    }

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
        sprintf("%.4g", $ctrlFreq),
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
        sprintf("%.2g", $coverageCont)
    );
    push (@output, join(",", sort @{$caseMatches->{'identMatchIDs'}}));
    push (@output, join(",", sort @{$caseMatches->{'NILMatchIDs'}})) if $ShowColumns{'NILCase'};
    push (@output, join(",", sort @{$caseMatches->{'NISMatchIDs'}})) if $ShowColumns{'NISCase'};
    push (@output, join(",", sort @{$caseMatches->{'NIMatchIDs'}})) if $ShowColumns{'NICase'};
    
    push (@output, join(",", sort @{$ctrlMatches->{'identMatchIDs'}}));
    push (@output, join(",", sort @{$ctrlMatches->{'NILMatchIDs'}})) if $ShowColumns{'NILCtrl'};
    push (@output, join(",", sort @{$ctrlMatches->{'NISMatchIDs'}})) if $ShowColumns{'NISCtrl'};
    push (@output, join(",", sort @{$ctrlMatches->{'NIMatchIDs'}})) if $ShowColumns{'NICtrl'};

    say $OUTPUTFILE join("\t", @output);

    # output the dump file if frequency is less than the filter set by user
    my $dumplines = $caseMatches->{'dumplines'};
    my $dumplinesCont = $ctrlMatches->{'dumplines'};

    say $DUMPFILE join("\n", @$dumplines) if ($caseFreq <= $FREQ_FILTER);
    say $DUMPFILE join("\n", @$dumplinesCont) if (scalar @$dumplinesCont && $ctrlFreq <= $FREQ_FILTER);

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
    my (%Matches, @matchesForDumpfile, %MatchIDs, $matchtype);
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
                    $MatchIDs{'ident'}{$testSampleID} = 1;
                }
            }
            # If not an 'identical' match, check for the three types of matches below:
            # For definitions of these matches, see comments in the synopsis.
            elsif ($match_case >= $MATCH_RATE) {
                $nonIdentLarger++;
                $matchtype = "Non_Identical_Larger";
                $MatchIDs{'NIL'}{$testSampleID} = 1;
            }
            elsif ($match_test >= $MATCH_RATE) {
                $nonIdentSmaller++;
                $matchtype = "Non_Identical_Smaller";
                $MatchIDs{'NIS'}{$testSampleID} = 1;
            }
            else {
                $nonIdentical++;
                $matchtype = "Non_Identical";
                $MatchIDs{'NI'}{$testSampleID} = 1;
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
    
    $Matches{identical} = $identical;
    $Matches{recip_identical} = $recip_identical;
    $Matches{nonIdentSmaller} = $nonIdentSmaller;
    $Matches{nonIdentLarger} = $nonIdentLarger;
    $Matches{nonIdentical} = $nonIdentical;
    $Matches{dumplines} = \@matchesForDumpfile;
    
    $Matches{identMatchIDs} = [keys %{$MatchIDs{'ident'}}];
    $Matches{NILMatchIDs} = [keys %{$MatchIDs{'NIL'}}];
    $Matches{NISMatchIDs} = [keys %{$MatchIDs{'NIS'}}];
    $Matches{NIMatchIDs} = [keys %{$MatchIDs{'NI'}}];

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
    #the command pod2usage prints the text given between the =head and =cut
    pod2usage(-message=>"$_[0]", -verbose=>0, -exitval=>2);
}
sub help( $ ){
    pod2usage(-message=>"", -verbose=>2, -exitval=>0);
}

=head1 NAME

B<CNV_counter> analyzes copy number variation (CNV) files to identify rare CNV's.

=head1 SYNOPSIS

 CNV_counter.pl --case CaseFile --control ControlFile

 CNV_counter.pl -c CaseFile -n ControlFile -o OutputPath -d DumpPath -m 85 -f 0.1 -g
 
 See the help page for descriptions of the options, the format of input 
 and output files, and definitions of match types.
 
 CNV_counter.pl --help

=head1 OPTIONS

=head2 Required:

=over 8

=item -c --case F<file>

Case CNV file, white space delimited

=item -n --control F<file>
 
Control CNV file, white space delimited

=back

=head2 Optional:

=over 4 

=item -m, --match [1-100] 

Percent for match to be considered 'identical' (defaults to 90)  

=item  -o, --out F<file>

Output file path (defaults to <Casefile>.matches.txt)

=item  -d, --dump F<file>

Dump file path (defaults to <Outfile>.matches_dump.zip)

=item  -f, --freq [0-1]

Frequency filter used to prevent bloated output and dump files.
If the match frequency is greater than this value for a case CNV:

   -Suppress listing ID's of identical matches in the last columns
    of the output file.
   -Suppress listing matches for the case CNV in the dump file.

=item  -a, --freqcase [0-1]

Frequency filter applied to matches from the Case file only

=item  -b, --freqcontrol [0-1] 

Frequency filter applied to matches from the Control file only

=item  --col [lsn]

Show additional columns with IDs of individuals with matching CNVs. 
By default only the IDs of individual with 'identical' CNV's are listed. 
Options l,s,n are: non-identical-[l]arger, non-identical-[s]maller 
and [n]on-identical.   
 
=item  --ccol [lsn]

Same as --col but only applies to matched IDs from the Case file
  
=item  --ncol [lsn]

Same as --col but only applies to matched IDs from the Control file
 
=item  --pos [lsn]

Choose which additional match types to consider as 'positive' matches.
By default, only 'identical' matches are considered positive matches 
(e.g. when calculating p-values).  Options l,s,n are: 
non-identical-[l]arger, non-identical-[s]maller and [n]on-identical. 
                        
=item  --cpos [lsn]

Same as --pos but only applies to matches from the Case file  
 
=item  --npos [lsn]

Same as --pos but only applies to matches from the Control file

=item -h, --help

Displays this page.

=back

=head1 DESCRIPTION

=head2 Format of Input Files

The input files must be whitespace delimited and must include the following columns:

  <chr:start-end> 
  <# of snps> 
  <# of bp's> 
  <copy number type> 
  <sample id> 
  <first SNP> 
  <last SNP> 
  <confidence score>

Example:

  chr16:32438272-32505195 91 66924 state2,cn=1 7 cnvi0020587 cnvi0020628 86.717

  Notes:
  - "number of snps"; "# of bp's"; "starting SNP"; "ending SNP"
      are not used but the column must exist in the input file

  - "copy number type" must contain <cn=#>, any other information in this field is ignored

  - "sample id" identifies the individual and can be any string, i.e. "123" or "person_1".
      It cannot be blank or 0.

=head2 Format of Output File

  Each line in the output file is tab delimited and corresponds to one line in the Case file.  

  * Case line:
      The entire original line from the case file.  The CNV in this line is the
      'reference' CNV.
  * Pos_Case (Identical case CNV's):
      Number of individuals in the Case file that have at least one 'identical' CNV.
  * Neg_Case (Cases Non-match):
      Number of individuals in the Case file that do not have at least one 'identical' CNV.
  * Freq_Case:
      Frequency of cases (individuals in the case file with at least one 'identical' CNV
      divided by total individuals in the Case file)
  * Pos_Ctrl:
      Number of CNV's in the Control file that are 'identical' to the reference CNV.
  * Neg_Ctrl:
      Number of individuals in the Control file that do not have at least one 'identical' CNV.
  * Freq_Ctrl:
      Frequency of controls (individuals in the Control file with at least one
      'identical' CNV divided by total individuals in the Controls file)
  * Fisher_p:
      P-value calculated using Fisher
  * Ident_Case:
      Number of CNV's in the Case file that are 'identical' to the reference CNV.
      There could be more than one match per individual if <match> percent is set to
      less than 50.
  * Reciprocal_Case:
      Number of identical matches in reciprocal CNV types.  For example, if the
      reference CNV is a deletion (copy number of 0 or 1), this is the count of
      matching CNV's that are duplications (copy number > 2)
  * NI_Larger_Cases (Non-identical Bigger):
      Number of case CNV's that contain <match> percent of the reference CNV
  * NI_Smaller_Cases (Non-identical Smaller):
      Number of CNV's in the Case file that are contained within the reference CNV
  * NI_Cases (Non-identical):
      Count of CNV's that do not meet any of the above matching criteria,
      however, there is at least some overlap with the reference CNV
  * CoverageMap_Case:
      Shows the number of matching CNV's at individual locations within the CNV
  * Coverage_Case:
      The fraction of the CNV that is covered by any case CNVs
  * [Control matches]:
      The previous 7 fields are repeated for matches found in the Control file.
  * ID_Matches_Cases / _Ctrls:
      Comma separated list of individual's that have at least one
      'identical' match
  * NILarger, NISmaller, NI_ID_Matches_Cases / _Ctrls:
      If --col, ccol, ncol flags are set, then additional lists of individuals 
      that have at least one of the corresponding match types. 

=head2 Format of Dump File:

  * Case CNV Position
  * Case Number of Base Pairs in CNV
  * Case CNV State
  * Case Sample ID
  * Match Type (defined below)
  * Matched CNV: File (Case or Control input file)
  * Matched CNV: Position
  * Matched CNV: Number of Base Pairs in CNV
  * Matched CNV: CNV State
  * Matched CNV: Sample ID
  * Matched CNV: Confidence Score

=head2 Match Types:

"Identical" - matches within the <match> percent specified.  Checks for a two-way match
where both the reference CNV and the tested CNV are both at least <match> percent
matches of each other

    |----------------|         [test start/end]
    TS               TE
     |---------------|         [ref start/end]
     S                E

 "Non-identical larger" - the reference CNV is contained within the test CNV

    |----------------|         [test start/end]
    TS               TE
        |-------|              [ref start/end]
        S       E

 "Non-identical smaller" - the test CNV is contained within the reference CNV

    |----------|          [test start/end]
    TS         TE
  |-----------------|     [ref start/end]
  S                 E

 "Non-identical" - none of the matches above apply,
 however, there is at least some overlap between the CNV's

    |----------|          [test start/end]
    TS         TE
          |---------|     [ref start/end]

=head2 Fisher Test
q
The Fishers exact test per CNV is a right-sided test for the contingency table

         Case  Control    | Totals
Present   a     b         | a+b
Absent    x     y         | x+y
__________________________|________
Totals    a+x   b+y       | a+b+x+y

Module:  Text::NSP::Measures::2D::Fisher::right

=head1 AUTHOR

David Fasel daf2139<at>columbia.edu 05/2015.  Adapted from code by Roel Sterken 6/2009.

=head1 VERSION

 0.20:
  - optimized with "interval tree" data type
  - zip dump file
 0.21:
  - fixed bug in Fisher p-val (because positive controls were always >= 1)
  - fixed bug in Fisher p-val (Frequency filter could change positive case/control values)
  - if extra columns in input file, put placeholder column headers in output file
  - added Coverage Map and Coverage Freq
  - option to set independent frequency filters for case and controls
  - output file is tab delimited only, even if input file uses mix of tabs and spaces
  - fixed bug some 'Non-Identical' Matches being classified as 'Non-Identical Smaller'
  - broke dumpfile (doesn't output lines)
 0.22
  - fixed dumpfile from 0.21
  - calculate reciprocal identical matches (deletions cn0,1 vs duplications cn3+)
  - Pos case/ctrl and "identical" case/ctrl are now 2 different columns (since they
    can be different values if <match> percent is less than 50)
  - Added -Greedy flag to count NI_larger as positive case/control
 0.23
  - Removed hack to decrement each Case coverage map position by 1.  Now Case coverage
    map will always show at least one match (since self is counted).
  - Require a control file to be specified (rather than defaulting to using case file
    as the control file which may be confusing to the user)
  - fixed a bug where --Greedy flag sometimes broke match lists
 0.24
  - Output a column which lists 'Non-identical Larger' match IDs when --Greedy is on
 0.25 
  - Fixed a bug where the --Greedy flag could count an individual more than once
 0.26
  - Added ability to show which Match ID cols are shown and which match types are 
    counted as a 'positive' match.  Removed "Greedy" flag since it was superseded by
    this change.
 0.27
  - Fixed a bug where -b option was using the value from -a.
  - fixed bug: if freq was 1, then filtering was applied even if turned off
  
 To Do:
 - suppress self matching in dump file?
 - count reciprocal matches for chrX (currently chrX is ignored since definition of a
   reciprocal match is more complicated for chrX for males - ie cn=1 is not a deletion).

=cut

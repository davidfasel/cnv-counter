# CNV Counter

Analyzes copy number variation (CNV) files to identify rare CNV's.

Usage

    CNV_counter.pl --case CaseFile --control ControlFile
    CNV_counter.pl -c CaseFile -n ControlFile -o OutputPath -d DumpPath -m 85 -f 0.1 -g
    CNV_counter.pl --help
   
   See the help within the app for descriptions of the options, the format of input 
   and output files, and definitions of match types.


## Description:

There are 2 main use cases:  One is to identify rare CNV by comparing a list of
patient CNVs with a list of control CNVs.  There are various match types
which can help enable this analysis.  For example, if I want to find rare, 
potential pathogenic CNVs, I might look for CNVs that don't or only partially 
overlap known CNVs within healthy populations.  The other use case is to
compare patient CNVs with known pathogenic CNVs.

### Match Types:

#### Identical
Matches within the <match> percent specified.  Checks for a two-way match
where both the reference CNV and the tested CNV are both at least <match> percent
matches of each other

    |----------------|         [test start/end]
    TS               TE
     |---------------|         [ref start/end]
     S                E

#### Non-identical larger
the reference CNV is contained within the test CNV

    |----------------|         [test start/end]
    TS               TE
        |-------|              [ref start/end]
        S       E

#### Non-identical smaller
The test CNV is contained within the reference CNV

    |----------|          [test start/end]
    TS         TE
  |-----------------|     [ref start/end]
  S                 E

#### Non-identical
None of the matches above apply,
 however, there is at least some overlap between the CNV's

    |----------|          [test start/end]
    TS         TE
          |---------|     [ref start/end]




## AUTHOR

David Fasel daf2139<at>columbia.edu 05/2015.  Adapted from code by Roel Sterken 6/2009.


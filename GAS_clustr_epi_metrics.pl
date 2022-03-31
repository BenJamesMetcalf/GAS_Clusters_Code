#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

sub checkOptions {
    my %opts;
    getopts('h:f:n:o:', \%opts);
    my ($help, $file, $outName, $outDir);

    if($opts{h}) {
        $help = $opts{h};
        help();
    }

    if (! $opts{f}) {
        print "No input cluster file provided (-f argument)\n";
        help();
    } elsif (-s $opts{f}) {
        $file = $opts{f};
        print "Cluster File: $file\n";
    } else {
        print "The input cluster file does not exist or is empty\n";
        help();
    }

    $outName = "epiFeat_clustr_metrics.txt";
    if($opts{n}) {
        $outName = $opts{n};
        chomp($outName);
        print "The output file name: $outName\n";
    } else {
        print "The default output file name: $outName\n";
    }

    $outDir = "./";
    if($opts{o}) {
        if (-d $opts{o}) {
            $outDir = $opts{o};
            print "The output directory is: $outDir\n";
        } else {
            $outDir = $opts{o};
            mkdir $outDir;
            print "The output directory has been created: $outDir\n";
        }
    } else {
        print "The files will be output into the current directory.\n";
    }

    return ($help, $file, $outName, $outDir);
}

sub help
{

die <<EOF

USAGE
loop.pl -t <Streptococcus type: Numeric> -f <Assembly list file: String> -n <Output File Name: String> -o <Output Directory: String>

    -h   print usage
    -f   assembly list file
    -n   output file name
    -o   output directory

EOF
}

my ($help, $file, $outName, $outDir) = checkOptions( @ARGV );


###Subroutines###
###Start Doing Stuff###
my $fileOut = "$outDir/$outName";
open ( my $outFile, ">>", "$fileOut" ) or die "Could not open file '$fileOut': $!";
print "Emm_Type,Cluster_ID,Cluster_Num,PEH_Count,PWID_Count,PEH/PWID_Count,LTCF_Count\n";
print $outFile "Emm_Type,Cluster_ID,Cluster_Num,PEH_Count,PWID_Count,PEH/PWID_Count,LTCF_Count\n";
open ( my $inFile, "<", "$file" ) or die "Could not open file '$file': $!";
local $/ = "";
while ( my $clustr = <$inFile> ) {
    chomp($clustr);
    #print "LINE_STARTS\n";
    #print "$clustr\n";
    my @rowsArr = split(/\n/,$clustr);
    my $headr = shift(@rowsArr);
    my @headr = split(':',$headr); #grab total cluster number
    my $headata = $headr[1];
    $headata =~ s/\n//g;
    #print "headr data: $headata\n";
    my @headrInfo = split(',',$headata);
    my @clustArr;
    foreach (@rowsArr) {
        #print "Line: $_\n";
        my @line = split(',',$_);
        push(@clustArr,\@line);
        #print "Sample: $clustArr[0][6]\n";
    }

    my $peh = 0;
    my $pwid = 0;
    my $ltcf = 0;
    my $peh_pwid = 0;
    foreach (@clustArr) {
        #print "Line :@$_ || Sample: @$_[7]\n";
        if (@$_[4] eq "Prev_Res") {
	    next;
	} else {
	    if (@$_[4] == 4) {
		$peh++;
	    }
	    if (@$_[4] == 2) {
		$ltcf++;
	    }
	    if (@$_[5] == 1) {
		$pwid++;
	    }
	    if (@$_[4] == 4 || @$_[5] == 1) {
		$peh_pwid++;
	    }
	}
    }
    print "$headrInfo[0], $headrInfo[1], $headrInfo[2], $peh, $pwid, $peh_pwid, $ltcf\n";
    print $outFile "$headrInfo[0], $headrInfo[1], $headrInfo[2], $peh, $pwid, $peh_pwid, $ltcf\n";
    #print "LINE_ENDS\n";
}
close $outFile;

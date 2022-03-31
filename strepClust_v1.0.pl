#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use File::Basename;
use Env;
#use lib $ENV{MODULESHOME}."/init";
use lib "/usr/share/Modules/init/";
use perl;
use Scalar::Util qw(looks_like_number);
use Cwd;
use Math::Round;
use List::Util qw(sum);
use Statistics::Basic qw(:all nofill);

sub checkOptions {
    my %opts;
    getopts('ht:f:n:m:o:', \%opts);
    my ($help, $type, $file, $outName, $month, $outDir);

    if($opts{h}) {
        $help = $opts{h};
        help();
    }

    if ((! $opts{t} && ! $opts{f}) || ($opts{t} && $opts{f})) {
        print "Program needs either a GAS type (-t) or user provided file of assembly paths (-f) as input arguments.\n";
        help();
    } elsif ($opts{t}) {
        if (looks_like_number($opts{t})) {
            $type = $opts{t};
            print "Streptococcus type: $type\n";
        } else {
            print "GAS type argument needs to be a number\n";
            help();
        }
    } elsif ($opts{f}) {
        if (-s $opts{f}) {
            $type = "User";
            $file = $opts{f};
            print "User provided file of assembly paths: $file\n";
        } else {
            print "The assembly paths file does not exist or is empty.\n";
            help();
        }
    }

    $outName = "./mash_nrDist.txt";
    if($opts{n}) {
        $outName = $opts{n};
        chomp($outName);
        print "The output file name prefix: $outName\n";
    }

    $month = "18";
    if($opts{m}) {
        $month = $opts{m};
        chomp($month);
        print "The interval of time (months) to search for clusters: $month\n";
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

    return ($help, $type, $file, $outName, $month, $outDir);
}

sub help
{

die <<EOF

USAGE
strepClust_v1.0.pl -t <Streptococcus type: Numeric> -f <Assembly list file: String> -n <Output File Name: String> -m <retrospective interval: Numeric> -o <Output Directory: String>

    -h   print usage
    -t   streptococcus type
    -f   assembly list file
    -n   output file name
    -m   retrospective interval (months)
    -o   output directory

EOF
}

my ($help, $type, $file, $outName, $month, $outDir) = checkOptions( @ARGV );


###Subroutines###
###Start Doing Stuff###
chdir "$outDir";
my $mcl_fin = "FINAL_m${type}_MCL_Clusters.txt";
open (my $fin, ">", $mcl_fin ) or die "Could not open file $mcl_fin: $!";
module "load Mash/1.1";
if ($type eq "User") {
    `cat $file > m${type}_paths.txt`;
} else {
    `mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select L_Contig_path,L_LABID from Main_GAS_2015thru2018_Ben where L_WGS_emm_subtype REGEXP '^${type}\$\|^${type}[.][0-9]+\$' AND (CULT2 IS NULL OR CULT2>= date_sub(now(),interval ${month} month))" | grep -v 'NULL\\|MN[0-9]\+\\|^contigs_path\\|auto_' | tr '\t' ',' | sed '1d' > m${type}_paths_names.txt`;
    `cat m${type}_paths_names.txt | cut -d, -f1 > m${type}_paths.txt`;
}
`mash sketch -l m${type}_paths.txt -o m${type}-s50K-k21_genomes -s 50000 -k 21 -p 10`;
`mash dist -t -p 10 m${type}-s50K-k21_genomes.msh m${type}-s50K-k21_genomes.msh | sed 's:/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/GAS_Typing_Output/GAS_2015/18-B-429_Illumina-MiSeq-M032../.\\{34\\}/GAS_Typing_Output/::g' | sed 's:/velvet_output/contigs.fa::g' > m${type}-s50K-k21_dist.txt`;

my @mash;
my $input = "m$type-s50K-k21_dist.txt";
open ( my $f_seq, "<", $input ) or die "Could not open file '$input': $!";
while ( my $line = <$f_seq> ) {
    chomp($line);
    my @line = split('\t',$line);
    push(@mash,\@line);
}
close $f_seq;

my @mashInfo;
for my $i (0 .. $#mash) {
    my $closeMatch = 0;
    my $ref = $mash[$i];
    my @arr = @$ref;
    foreach my $dist (@arr) {
        if (looks_like_number($dist) && $dist < 0.00001) { #0.000031
            $closeMatch++;
        }
    }
    print "MASHINFO | index $i, Sample: $mash[$i][0], Matches: $closeMatch\n";
    push(@mashInfo,[$i, $mash[$i][0], $closeMatch]);
}

foreach my $row (@mashInfo) {
    my @arr = @$row;
    my $sample = $arr[1];
    my $matches = $arr[2];
    if ($matches >= 2) {
        `grep '${sample}' m${type}_paths_names.txt >> m${type}_nucmer_paths.txt`;
    }
}

=pop
my @sortInfo = sort { $b->[2] <=> $a->[2] } @mashInfo;
my @mashSort;
foreach my $row (@sortInfo) {
    my @arr = @$row;
    my $index = $arr[0];
    push(@mashSort,$mash[$index]);
}

my @mashFin;
for my $i (0 .. $#mashSort) {
    my $sortMatch = 0;
    for my $j ($i+1 .. $#mashSort) {
        my $dist = $mashSort[$i][$j];
        if (looks_like_number($dist) && $dist < 0.00001) { #0.000031
            $sortMatch++;
        }
    }
    if ($sortMatch >= 1) {
        #print "SORT-MASH | index $i, Sample: $mashSort[$i][0], Matches: $sortMatch\n";
        my $sample = $mashSort[$i][0];
        `grep '${sample}' m${type}_paths_names.txt >> m${type}_nucmer_paths.txt`;
    }
}
=cut

my $cur = getcwd();
my $tasks = `cat m${type}_nucmer_paths.txt | wc -l` - 1;
print "task value is: $tasks\n";
my $out_qsub = "$cur/qsub_output";
`mkdir ${out_qsub}`;
###Send the jobs out on the cluster with each sample running in parallel###
`qsub -sync y -q all.q -t 1-${tasks} -cwd -o "$out_qsub" -e "$out_qsub" /scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/StrepLab_Clustr_Detection/Clustr_Tools/nucmer_dist_qsub.sh $cur/m${type}_nucmer_paths.txt $cur/m${type}_paths_names.txt`;

my $nucmer_out = "m${type}_nucmer-dist_ADJ.txt";
#my $nucmer_out2 = "m${type}_nucmer-dist_ARR.txt";
for (my $i=1; $i <= $tasks; $i++) {
    `cat nucmer-"${i}"_dist.txt >> "${nucmer_out}"`;
    `rm nucmer-"${i}"_dist.txt`;
}

module "purge";
print "Beginning MCL Clustering\n";
module "load MCL-edge/14-137";
my $mcl_name = "m${type}_MCL-graph";
`mcxload -abc $nucmer_out -write-tab $mcl_name.tab -o $mcl_name.mci --stream-mirror -stream-tf 'gt(0.004)'`;
`mcl  $mcl_name.mci -I 1.3`; #e^-((x/2.5)-0.5)

#Generate pairwise dist matrix and print output#
my @mcl2sampl;
my %sampl2mcl;
open ( my $m_tab, "<", "m${type}_MCL-graph.tab" ) or die "Could not open file 'm${type}_MCL-graph.tab': $!";
while ( my $line = <$m_tab> ) {
    chomp($line);
    my @ind = split('\t',$line);
    $mcl2sampl[$ind[0]] = $ind[1];
    $sampl2mcl{$ind[1]} = $ind[0];
}
close $m_tab;

module "purge";
module "load MUMmer/3.9.4";
my @mcl_clust;
open ( my $m_out, "<", "out.m${type}_MCL-graph.mci.I13" ) or die "Could not open file 'out.m${type}_MCL-graph.mci.I13': $!";
while ( my $line = <$m_out> ) {
    chomp($line);
    if ($line =~ /begin/) {
        local $/ = "\$";
        while (<$m_out>) {
            my $clustr = $_;
            last if /\)/;
            my @clust = split(' ',$clustr);
            pop(@clust);
            my $clustNum = shift(@clust);
            my $len = scalar(@clust);
            if ($len >= 4) {
                my @distOut;
                my @snpDist;
                push @distOut, [("N")x($len+9)] for (0..$len);
                $distOut[0][0] = "State";
		$distOut[0][1] = "Zip_Code";
		$distOut[0][2] = "Emm_subtype";
                $distOut[0][3] = "Age";
                $distOut[0][4] = "Prev_Res";
                $distOut[0][5] = "IVDRUG";
                $distOut[0][6] = "Cult_Date";
                $distOut[0][7] = "Sample";
                $distOut[0][8] = "Index";
                $distOut[0][9] = 0;
                my $length = scalar(@distOut);
                my $length2 = scalar(@{$distOut[1]});
                #print "rows: $length || columns: $length2 || cluster0 var: $clust[0]\n";
                #Sort cluster samples by cluster date
                my @clustLabl;
                for my $i (0 .. $#clust) {
                    my $sampl = "'".$mcl2sampl[$clust[$i]]."'";
                    push(@clustLabl,$sampl);
                }
                my $clustLabls = join(",",@clustLabl);
                #print "Cluster Labels:\n$clustLabls\n";
		my $sortClustStr = `mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select L_LABID,CULT2,PREVRES,IVDRUG,AGEYRS,L_WGS_emm_subtype,STATE,ZIP from Main_GAS_2015thru2018_Ben where L_LABID In ($clustLabls) ORDER BY CULT2" | cut -f1-8 --output-delimiter=, | sed '1d'`;

                #print "sorted cluster string:\n$sortClustStr\n";
                my @sortClust = map { [split/,/] } (split/\n/,$sortClustStr);

		my $dateClustStr = `mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select CULT2 from Main_GAS_2015thru2018_Ben where L_LABID In ($clustLabls) ORDER BY CULT2" | sed '1d' | grep -v "NULL"`;
		my @dateClust = split(/\n/,$dateClustStr);
		my $dur;
		if (@dateClust) {
		    my $date1 = $dateClust[0];
		    my $date2 = $dateClust[-1];
		    if (scalar(@dateClust) == scalar(@sortClust)) {
			$dur = `mysql -u streplab -h xxxx -pxxxx -e "SELECT DATEDIFF('${date2}','${date1}')" | sed -n "2p"`;
		    } else {
			my $day = `mysql -u streplab -h xxxx -pxxxx -e "SELECT DATEDIFF('${date2}','${date1}')" | sed -n "2p"`;
			$dur = ">="."$day";
		    }
		} else {
		    $dur = "NULL";
		}
		$dur =~ s/\n//g;
                for my $i (0 .. $#sortClust) {
                    #print "sortClust\n";
                    #print "$sortClust[$i][0] || $sortClust[$i][1] || $sortClust[$i][2] || $sortClust[$i][3]\n";
                    $distOut[$i+1][0] = $sortClust[$i][6];#state
		    $distOut[$i+1][1] = $sortClust[$i][7];#zip
		    $distOut[$i+1][2] = $sortClust[$i][5];#emm_subtype
                    $distOut[$i+1][3] = $sortClust[$i][4];#age
                    $distOut[$i+1][4] = $sortClust[$i][2];#prev_res
                    $distOut[$i+1][5] = $sortClust[$i][3];#ivduc
                    $distOut[$i+1][6] = $sortClust[$i][1];#cult_date
                    $distOut[$i+1][7] = $sortClust[$i][0];#sample
                    $distOut[$i+1][8] = $i;
                    #print "distout\n$distOut[$i+1][0] || $distOut[$i+1][1] || $distOut[$i+1][2] || $distOut[$i+1][3] || $distOut[$i+1][4] || $distOut[$i+1][5] || $distOut[$i+1][6]\n";
                    for my $j ($i+1 .. $#sortClust) {
                        $distOut[0][$j+9] = $j;
                        my $sampl1 = $sortClust[$i][0];
                        my $sampl2 = $sortClust[$j][0];
                        my $dist = `grep "$sampl1" m${type}_nucmer-dist_ADJ.txt | grep "$sampl2" | cut -d' ' -f3 | head -n1`;
                        $dist =~ s/\n//g;
                        if (looks_like_number($dist) && $dist > 0) {
                            my $snps = ((-log($dist)+1)*2.3);
                            my $snps2 = round($snps);
                            $distOut[$i+1][$j+9] = $snps2;
                            push(@snpDist,$snps2);
                        } else {
                            ##Calculate nucmer distance##
                            print "Need to calculate distance: $sampl1 vs $sampl2\n";
                            my $genome1 = `grep '${sampl1}' m${type}_paths_names.txt | cut -d, -f1`; #`grep '${sampl1}' m${type}_paths_names.txt | cut -d, -f1`
                            my $genome2 = `grep '${sampl2}' m${type}_paths_names.txt | cut -d, -f1`;
                            $genome1 =~ s/\n//g;
                            $genome2 =~ s/\n//g;
                            print "genome 1: $genome1 || genome 2: $genome2\n";
                            `nucmer --mum --prefix=TEMP_${type}_nucDist ${genome1} ${genome2}`;
                            `delta-filter -1 TEMP_${type}_nucDist.delta > TEMP_${type}_nucDist-Fg.delta`;
                            if (-s "TEMP_${type}_nucDist-Fg.delta") {
                                #my $snps2 = `show-snps -CIH TEMP_${type}_nucDist-Fg.delta | wc -l`;
                                #my $cmd = "show-snps -CHI TEMP_${type}_nucDist-Fg.delta | awk '$6 >= 100 {print $0}' | wc -l";
                                #print "$cmd\n";
				my $snps = `show-snps -CHI TEMP_${type}_nucDist-Fg.delta | awk '\$6 >= 100 {print \$0}' | wc -l`; #Changed from 1000
                                $snps =~ s/\n//g;
                                push(@snpDist,$snps);
                                $distOut[$i+1][$j+9] = $snps;
                                `rm TEMP_${type}_nucDist.delta TEMP_${type}_nucDist-Fg.delta`;
                            } else {
                                $distOut[$i+1][$j+9] = "N";
                            }
                        }
                    }
                }

                #my $avgDist = sum(@snpDist)/@snpDist;
		my $avgDist = mean(@snpDist);
		my $stdDev = stddev(@snpDist);
                #print "Emm_Type: $type, Cluster: $clustNum, Samples: $len, Avg Pairwise: $avgDist, Std Dev: $stdDev, Duration: $dur\n";
                #print $fin "Emm_Type: $type, Cluster: $clustNum, Samples: $len, Avg Pairwise: $avgDist, Std Dev: $stdDev, Duration: $dur\n";
		print "Emm_Type,Cluster,Sample_Num,Avg_Pair,Std_Dev,Dur:$type,$clustNum,$len,$avgDist,$stdDev,$dur\n";
		print $fin "Emm_Type,Cluster,Sample_Num,Avg_Pair,Std_Dev,Dur:$type,$clustNum,$len,$avgDist,$stdDev,$dur\n";
                for my $row (@distOut) {
                    print join(",",@{$row}), "\n";
                    print $fin join(",",@{$row}), "\n";
                }
                print "\n";
                print $fin "\n";
            }
        }
    }
}
=cut

#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

#emmTotal = total count (cluster and non cluster) for focal emm type
#emmClust = count of cluster isolates for focal emm type
#emmBackground = count of non-cluster isolates for focal emm type
#countClust = total count of cluster isolates for all emm types
#countGAS = total count of GAS

#emmOdds = count of focal cluster isolates / count of focal non-cluster isolates
#refOdds = total count of cluster isolates for all emm types minus focal emm type / total count of GAS minus total cluster samples minus focal non-cluster isolates
#emmOR = focal emm type odds ratio / reference odds ratio

sub checkOptions {
    my %opts;
    getopts('hc:n:o:', \%opts);
    my ($help, $clustFile, $outName, $outDir);

    if($opts{h}) {
        $help = $opts{h};
        help();
    }

    if($opts{c}) {
        $clustFile = $opts{c};
        if (-e $clustFile) {
            print "File containing the predicted GAS clusters: $clustFile\n";
        } else {
            print "The location given for the file containing GAS clusters is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/cluster_file).\n";
            help();
        }
    } else {
        print "The location of the GAS clusters file (including full path) has not been given.\n";
        help();
    }

    $outName = "FINAL_emmIOR_emmPrev.txt";
    if($opts{n}) {
        $outName = $opts{n};
        chomp($outName);
        print "The output file name prefix: $outName\n";
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

    return ($help, $clustFile, $outName, $outDir);
}

sub help
{

die <<EOF

USAGE
clustrOR_creatr.pl -c <Cluster file: String> -n <Output File Name: String> -o <Output Directory: String>

    -h   print usage
    -c   cluster file
    -n   output file name
    -o   output directory

EOF
}

my ($help, $clustFile, $outName, $outDir) = checkOptions( @ARGV );


###Subroutines###
###Start Doing Stuff###
chdir "$outDir";
my $fileOut = "$outDir/$outName";
open ( my $outFile, ">>", "$fileOut" ) or die "Could not open file '$fileOut': $!";
print "Emm_Type,Emm_Cluster_Num,Emm_Total_Num,Emm_NonCluster_Num,Emm_Odds,Emm_OR,Emm_Prevalence,AGEYRS_AVG,WEIGHTLB_AVG,BLACK_SUM,WHITE_SUM,Capsule_SUM,SDA1_SUM,SLAA_SUM,SIC_SUM,ROCA_SUM,PNGA_SUM,NADase_D330G_SUM,PWID_SUM,LTCF_SUM,PEH_SUM,VULN_SUM,ExoA_SUM,ExoC_SUM,ExoG_SUM,ExoH_SUM,ExoI_SUM,ExoJ_SUM,ExoK_SUM,ExoL_SUM,ExoM_SUM,ExoS_SUM,ExoZ_SUM,PRTF2_SUM,R28_SUM,SFB1_SUM,SOF_SUM\n";
print $outFile "Emm_Type,Emm_Cluster_Num,Emm_Total_Num,Emm_NonCluster_Num,Emm_Odds,Emm_OR,Emm_Prevalence,AGEYRS_AVG,WEIGHTLB_AVG,BLACK_SUM,WHITE_SUM,Capsule_SUM,SDA1_SUM,SLAA_SUM,SIC_SUM,ROCA_SUM,PNGA_SUM,NADase_D330G_SUM,PWID_SUM,LTCF_SUM,PEH_SUM,VULN_SUM,ExoA_SUM,ExoC_SUM,ExoG_SUM,ExoH_SUM,ExoI_SUM,ExoJ_SUM,ExoK_SUM,ExoL_SUM,ExoM_SUM,ExoS_SUM,ExoZ_SUM,PRTF2_SUM,R28_SUM,SFB1_SUM,SOF_SUM\n";

my $countGAS = `mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select count(*) as Total from Main_GAS_2015thru2018_Ben where L_WGS_emm_subtype IS NOT NULL and L_WGS_emm_subtype NOT REGEXP '[A-Za-z]'" | sed -n '2p'`;
chomp($countGAS);
my $countClust = `cat "$clustFile" | grep -E -v '^Emm_Type|^State' | sed '/^\$/d' | wc -l`;
chomp($countClust);
#print "GAS isolates: $countGAS | Cluster isolates: $countClust\n";
#my @emmArr = `cat "$clustFile" | grep -E -v '^Emm_Type|^State' | sed '/^\$/d' | cut -d, -f3 | sort | uniq`;
#my @emmArr = `cat "$clustFile" | grep -E -v '^Emm_Type|^State' | sed '/^\$/d' | cut -d, -f3 | sed 's/\\([0-9]\\+\\)\\.[0-9]\\+/\\1/g' | sort | uniq`;
my @emmArr = `mysql -Ns -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select distinct L_WGS_emm_subtype from Main_GAS_2015thru2018_Ben where L_WGS_emm_subtype IS NOT NULL and L_WGS_emm_subtype NOT REGEXP '[A-Za-z]'" | sed 's/\\([0-9]\\+\\)\\.[0-9]\\+/\\1/g' | sort | uniq`;
my @states = ('MD', 'TN', 'MN', 'CO', 'CT', 'CA', 'GA', 'NM', 'NY', 'OR');

foreach my $emmType (@emmArr) {
    chomp($emmType);
    #Calculate Cluster Odds Ratios and Emm Subtype Prevalence Metrics for Each Emm Subtype
    my $emmTotal = `mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select count(*) as Total from Main_GAS_2015thru2018_Ben where L_WGS_emm_subtype REGEXP '^$emmType\[.\]' and L_WGS_emm_subtype NOT REGEXP '[A-Za-z]'" | sed -n '2p'`;
    chomp($emmTotal);
    if ($emmTotal >= 10) {
	my $emmClust = `cat "$clustFile" | grep -E -v '^Emm_Type|^State' | sed '/^\$/d' | awk -F"," -v pat="^$emmType\[.\]" '\$3~pat {print \$0}' | wc -l`;
	chomp($emmClust);
        if (! $emmClust) {
            $emmClust = 0;
            #print "emm Type: $emmType | No Cluster Isolates\n";
        } else {
            #print "emm Type: $emmType | Cluster: $emmClust\n";
        }
        #Calculate Cluster Odds Ratios and Emm Subtype Prevalence Metrics for Each Emm Type
        my %stateNum;
        foreach my $site (@states) {
            my $emmState = `mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select count(*) as Total from Main_GAS_2015thru2018_Ben where STATE='$site' and  L_WGS_emm_subtype REGEXP '^$emmType\[.\]' and L_WGS_emm_subtype NOT REGEXP '[A-Za-z]'" | sed -n '2p'`;
            chomp($emmState);
	    $stateNum{$site} = $emmState;
        }
	#my $emmPrev = ($stateNum{MD}[0]*0.0871 + $stateNum{TN}[0]*0.0567 + $stateNum{MN}[0]*0.1748 + $stateNum{CO}[0]*0.0864 + $stateNum{CT}[0]*0.1132 + $stateNum{CA}[0]*0.1157 + $stateNum{GA}[0]*0.1777 + $stateNum{NM}[0]*0.0659 + $stateNum{NY}[0]*0.0659 + $stateNum{OR}[0]*0.0567);
	my $emmPrev = ($stateNum{MD}*0.0871 + $stateNum{TN}*0.0567 + $stateNum{MN}*0.1748 + $stateNum{CO}*0.0864 + $stateNum{CT}*0.1132 + $stateNum{CA}*0.1157 + $stateNum{GA}*0.1777 + $stateNum{NM}*0.0659 + $stateNum{NY}*0.0659 + $stateNum{OR}*0.0567);
        #print "Emm Type: $emmType | Emm Prevalence: $emmPrev\n";
        my $emmBackground = $emmTotal - $emmClust;

        my $emmOdds;
        my $emmOR;
        if ($emmClust > 0 and $emmBackground > 0) {
            $emmOdds = $emmClust / $emmBackground;
            my $refOdds = ($countClust - $emmClust) / ($countGAS - $countClust - $emmBackground);
            $emmOR = $emmOdds / $refOdds;
        } else {
            $emmOdds = ($emmClust + 0.5) / ($emmBackground + 0.5);
            my $refOdds = ($countClust - $emmClust + 0.5) / ($countGAS - $countClust - $emmBackground + 0.5);
            $emmOR = $emmOdds / $refOdds;
        }
        print "$emmType,$emmClust,$emmTotal,$emmBackground,$emmOdds,$emmOR,$emmPrev,";
	print $outFile "$emmType,$emmClust,$emmTotal,$emmBackground,$emmOdds,$emmOR,$emmPrev,";

	my $attr1 = `mysql -Ns -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select AVG(AGEYRS), AVG(WEIGHTLB), SUM(case when BLACK = 1 then 1 else 0 end) as BLACK, SUM(case when WHITE = 1 then 1 else 0 end) as WHITE, SUM(case when L_CAPSULE = 'HASA' then 1 else 0 end) as Capsule, SUM(case when L_SDA1 = 'SDA1' then 1 else 0 end) as SDA1, SUM(case when L_SLAA = 'SLAA' then 1 else 0 end) as SLAA, SUM(case when L_SIC = 'SIC' then 1 else 0 end) as SIC, SUM(case when L_ROCA REGEXP 'ROCAM*' then 1 else 0 end) as ROCA, SUM(case when L_PNGA3 = 'PNGA3' then 1 else 0 end) as PNGA, SUM(case when L_NADase_D330G = '330G' then 1 else 0 end) as NADase_D330G, SUM(IVDRUG),SUM(case when PREVRES = 2 then 1 else 0 end) as LTCF,SUM(case when PREVRES = 4 then 1 else 0 end) as PEH, SUM(case when (PREVRES = 4 OR IVDRUG = 1) then 1 else 0 end) as VULN from Main_GAS_2015thru2018_Ben where L_WGS_emm_subtype REGEXP '^$emmType\[.\]' and L_WGS_emm_subtype NOT REGEXP '[A-Za-z]'" | tr '\t' ','`;
	chomp($attr1);
	print "$attr1,";
	print $outFile "$attr1,";
       
	#Exotoxin Proteins
	my $exo_all = `mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select L_EXOTOXINS from Main_GAS_2015thru2018_Ben where L_WGS_emm_subtype REGEXP '^$emmType\[.\]' and L_WGS_emm_subtype NOT REGEXP '[A-Za-z]'" | tr ':' '\n' | grep -v 'L_EXOTOXINS\\\|neg\\\|NULL' | sort -k2,2 | uniq -c | sed 's/^[ ]*//' | awk -F" " '{print \$2","\$1}' | tr '\n' ':'`;
	chomp($exo_all);
	#print "EXOTOXINS: $exo_all\n";
	my %exoHash = split /[:,]/, $exo_all;
	my @exoTypes = ("A", "C", "G", "H", "I", "J", "K", "L", "M", "S", "Z");
	my %exoOut;
	foreach my $tox (@exoTypes) {
	    if (exists($exoHash{$tox})) {
		$exoOut{$tox} = $exoHash{$tox};
		#print "Exotoxin $tox Exists: $exoHash{$tox}\n";
	    } else {
		$exoOut{$tox} = 0;
		#print "Exotoxin $tox Does Not Exist\n";
	    }
	}
	#print map { "$_ => $exoOut{$_}\n" } keys %exoOut;
	print "$exoOut{'A'},$exoOut{'C'},$exoOut{'G'},$exoOut{'H'},$exoOut{'I'},$exoOut{'J'},$exoOut{'K'},$exoOut{'L'},$exoOut{'M'},$exoOut{'S'},$exoOut{'Z'},";
	print $outFile "$exoOut{'A'},$exoOut{'C'},$exoOut{'G'},$exoOut{'H'},$exoOut{'I'},$exoOut{'J'},$exoOut{'K'},$exoOut{'L'},$exoOut{'M'},$exoOut{'S'},$exoOut{'Z'},";

	#Other Surface Proteins
	my $miscSurf = `mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select L_OTHER_SURFACE_PROTEINS from Main_GAS_2015thru2018_Ben where L_WGS_emm_subtype REGEXP '^$emmType\[.\]' and L_WGS_emm_subtype NOT REGEXP '[A-Za-z]'" | tr ':' '\n' | grep -v 'L_OTHER_SURFACE_PROTEINS\\\|neg\\\|NULL' | sort -k2,2 | uniq -c | sed 's/^[ ]*//' | awk -F" " '{print \$2","\$1}' | tr '\n' ':'`;
	chomp($miscSurf);
	#print "Misc Surface: $miscSurf\n";
	my %surfHash = split /[:,]/, $miscSurf;
	my @surfTypes = ("PRTF2", "R28", "SFB1", "SOF");
	my %surfOut;
	foreach my $prot (@surfTypes) {
	    if (exists($surfHash{$prot})) {
		$surfOut{$prot} = $surfHash{$prot};
		#print "Surface Protein $prot Exists: $surfHash{$prot}\n";
	    } else {
		$surfOut{$prot} = 0;
		#print "Surface Protein $prot Does Not Exist\n";
	    }
	}
	#print map { "$_ => $surfOut{$_}\n" } keys %surfOut;
	print "$surfOut{'PRTF2'},$surfOut{'R28'},$surfOut{'SFB1'},$surfOut{'SOF'}\n";
	print $outFile "$surfOut{'PRTF2'},$surfOut{'R28'},$surfOut{'SFB1'},$surfOut{'SOF'}\n";
    }
}

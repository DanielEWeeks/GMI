#
#   GMI: Genetic Map Interpolator
#   Copyright (C) 2010 Nandita Mukhopadhyay, Xinyu Tang, Daniel E. Weeks
# 
#   This file is part of the GMI program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option) any later
#   version.
# 
#   GMI is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
#
# function to create a genetic map for a user-provided 
# list of SNPs. It uses merged SNP physical and genetic positions
# from Ensembl and Rutgers combined maps, using linear interpolation
# for SNPs that are not already mapped.
# 
# Usage: gmi_bioc_query.pl snp_list_file [additional args]
# Arguments:
#           --namecol [col], default 1
#           --poscol [col], default = -1, meaning absent
#           --chrcol [col], default = -1 meaning absent
#           --header (no arguments), first line will be skipped if specified
#           --noplots (no arguments), plots will be drawn if specified
#           --newsnps (no arguments), new SNP names will be used if specified
#           --outfile [output-file-prefix], default snp_mapped
#
# Configuration parameters:
# $mpath
# $wget_prog
# $gmi_version
# ========================insert parameters here=============


$| = 1;

if (scalar(@ARGV) == 0) {
    print_usage();
    exit(0);
}

use warnings;
no warnings qw(once numeric);
use Getopt::Long;
use File::Copy;
use Fcntl  qw(:flock);

open(UTILS, "$mpath/gmi_bioc_utils.pl")  || 
    die("Unable to find $mpath/gmi_bioc_utils.pl");
close(UTILS);

require "$mpath/gmi_bioc_utils.pl";

my $output_file = 'snp_mapped';
my $snpcol = 1;
my $has_header = 0;
my $chrcol = -1;
my $poscol = -1;
my $skip_plots = 0;
my $new_snps = 0;

$snpfile = $ARGV[0];

# Check if the snplist file exists;

open(SNPDATA, $snpfile) || die "Could not open $snpfile!\n";
flock(SNPDATA, LOCK_EX);

my @snplines = <SNPDATA>;
flock(SNPDATA, LOCK_UN);
close(SNPDATA);

GetOptions('namecol=i' => \$snpcol,
           'chrcol=i' => \$chrcol,
           'poscol=i' => \$poscol,
	   'header'=> \$has_header,
	   'outfile=s'=> \$output_file,
	   'noplots'=> \$skip_plots,
	   'newsnps'=> \$new_snps
    );

if ($poscol > -1 && $chrcol == -1) {
    print "Warning: both position and chromosome must be provided,\n";
    print "         ignoring positions given in input file.\n";
    $poscol = -1;
}

if ($chrcol > -1 && $poscol == -1) {
    print "Warning: both position and chromosome must be provided,\n";
    print "         ignoring chromosome numbers given in input file.\n";
    $chrcol = -1;
}

if ($has_header && scalar @snplines < 2) {
    print "Error: $snpfile does not contain data (only has a header line).\n";
    print("Terminating GMI.");
    exit(-1);
}

# print "Checking wget program\n";

# if (&check_installed_utils("wget") == 0) {
#     die "Error: wget program not found\n";
# }


########################################################################
# Check Ensembl current release version
########################################################################
# Currently used Ensembl release;
# NM - copy the global files into the user's local space, then 
# use those

# set global path and registry

print "GMI version $gmi_version\n";
print "Checking latest Ensembl release version against installed version:\n";
if (-f "version.txt") {
    unlink("version.txt");
}
system("R --slave < $mpath/gmi.bioc.checkver.R");
if (!-f "version.txt") {
    print("Unable to connect to Biomart.\n");
    print("Terminating GMI.");
    exit(-1);
}

open(VH, "version.txt") || die("Unable to open version.txt");

$ver = <VH>;
chomp($ver);
$ver = $ver * 1.0;
close(VH);

print "The latest Ensembl release version is $ver \n";

open(OVER, "$mpath/version.txt") || die("Unable to open $mpath/version.txt");

@ov1 = <OVER>;
$over = $ov1[0];
chomp($over);

print "The installed Ensembl version is $over \n";

if($ver != $over) {
    print "Local databases need to be updated,\n";
    print "continue with GMI query? [answer \"yes\" or \"no\"] (default no) > ";
    $ans = <STDIN>;
    chomp($ans);

    if ($ans ne "yes") {
	print("Error: up-to-date check failed for Ensembl version!\n");
        print("To update, run gmi_bioc_install.pl inside the GMI install folder again.\n");
	print("Terminating gmi query.\n");
	exit(-1);
    }
}

########################################################################
# Print some run-related information
########################################################################

my $wfile = $output_file."_log.txt";

# Print the first few lines of the input file to show
# how the data is being read in.

open(WAR, ">$wfile") || die("Cannot open $wfile");

printf WAR "---------------------------------------------------\n";
printf WAR "   GMI version $gmi_version\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday, @rest) = localtime(time);
$year += 1900; ## $year contains no. of years since 1900, to add 1900 to make Y2K compliant
my @month_abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );

printf WAR "   Run date:        $wday $month_abbr[$mon] $mday $hour:$min:$sec $year\n";
printf WAR "   Input file:      $snpfile\n";
printf WAR "   GMI options:\n";

if ($new_snps) {
    printf WAR "      Output map files contain merged NCBI SNP names.\n";
}
else {
    printf WAR "      Output map files contain user's SNP names.\n";
}
  
if (!$skip_plots) {
    printf WAR "      Graphical output files created.\n";
}
 
printf WAR "---------------------------------------------------\n\n";

if (scalar @snplines <= 5) {
    if (! $has_header) {
	@flines = @snplines;
    }
    else {
	@flines = @snplines[1 .. scalar @snplines - 1];
    }
}
else {
    if (! $has_header) {
	@flines = @snplines[0 .. 4];
    }
    else {
	@flines = @snplines[1 .. 5];
    }
}

$flines_num = scalar @flines;

@filehandles=qw/STDOUT WAR/;


foreach $fl (@filehandles) {

    printf $fl "Reading in data as per options to GMI:\n";
    if ($has_header) {
	printf $fl  "   Input SNP file has a header line.\n";
    }
    else {
	printf $fl  "   Input SNP file does not have a header line.\n";
    }
    printf $fl "   Reading SNP IDs from column $snpcol";

    if ($chrcol > -1 && $poscol > -1) {
	printf $fl 
	    ", chromosomes from column $chrcol, user positions from column $poscol";
    }
    printf $fl "\n";
    printf $fl "---------------------------------------------------\n\n";    
    printf $fl  "First $flines_num lines of data read in from $snpfile:\n";
    printf $fl  "   SNP-name(column $snpcol)";

    if ($chrcol > -1 && $poscol > -1) {
	printf $fl ", Chromosome(column $chrcol)\tPosition(column $poscol)";
    }
    printf $fl "\n";

    foreach $line (@flines) {
	chomp($line);
	@word = split(/\s+/, $line);
	if (scalar @word < $snpcol) {
	    next;
	}
	printf $fl "   ".$word[$snpcol-1];
	if ($chrcol > -1 && $poscol > -1) {
	    printf $fl "\t\t".$word[$chrcol - 1]."\t\t\t".$word[$poscol - 1];
	}
	printf $fl "\n";
    }
    printf $fl "---------------------------------------------------\n\n";
}

close(WAR);

########################################################################
# Obtain the list of snp ids
########################################################################

# First create just the list of SNPs

my @snpid = ();

# This is not necessary, we are not writing anything here,
# just processing lists

my $i=1;

my %seen = (); 
my %user_chromosomes = ();
my %user_positions = ();

# Modified to check for duplicates as well as store snps,
# chromosomes and positions
my @duplicates = ();
my %chr_pos_undef_dict = ();
my @chr_pos_undef = ();
my %lowered = ();
my @found_RS = ();

$lnum = 1;
foreach $line (@snplines)
{
    if ($i && $has_header) {
	$i = 0;
    }
    else {
	chomp($line);
	@word = split(/\s+/, $line);
	if (scalar @word < $snpcol) {
	    print "WARNING: Line $lnum has fewer than $snpcol columns, skipping.\n";
	    $lnum++;
	    next;
	}
	$snp1 = $word[$snpcol-1];
	$snpname = $snp1;

	$snpname=~ s/^[R]/r/ ; $snpname=~ s/^rS/rs/ ;
	$lnum++;
	
	if($snp1 ne $snpname) {
	    push(@found_RS, $snp1);
	    $lowered{$snp1} = $snpname;
	}

	if (!exists($seen{$snpname}))  {
	    # First time we've seen this one
	    $seen{$snpname} = 0;
	    push(@snpid, $snpname);
	    
	    if ($chrcol > -1 ) {
		if ($word[$chrcol-1] eq ($word[$chrcol - 1] + 0)) {
		    # This is a test for numeric values
		    if ($word[$chrcol-1] eq "23") {
			$user_chromosomes{$snpname}  = "X";
		    }
		    elsif ($word[$chrcol-1] eq "24") {
			$user_chromosomes{$snpname}  = "Y";
		    }
		    elsif ($word[$chrcol-1] eq "25") {
			$user_chromosomes{$snpname}  = "XY";
		    }
#		    elsif ($word[$chrcol-1] eq "24") {
#			$user_chromosomes{$snpname}  = "XY";
#		    }
#		    elsif ($word[$chrcol-1] eq "25") {
#			$user_chromosomes{$snpname}  = "Y";
#		    }
		    elsif ($word[$chrcol-1] eq "26") {
			$user_chromosomes{$snpname}  = "MT";
		    }
		    elsif ($word[$chrcol-1] eq "999") {
			$user_chromosomes{$snpname}  = "U";
		    }
		    else {
			if (($word[$chrcol - 1] + 0) >= 1 &&
			    ($word[$chrcol - 1] + 0) <= 22) {
			    $user_chromosomes{$snpname} = $word[$chrcol-1];
			}
			else {
			    $user_chromosomes{$snpname}  = "U";
			}
		    }
		}
		else {
		    if ($word[$chrcol-1] eq "U" ||
			$word[$chrcol-1] eq "X" ||
			$word[$chrcol-1] eq "Y" ||
			$word[$chrcol-1] eq "XY" ||
			$word[$chrcol-1] eq "MT") {
			$user_chromosomes{$snpname} = $word[$chrcol-1];
		    }
		    else {
			$user_chromosomes{$snpname}  = "U";
		    }		
		}

		if ($user_chromosomes{$snpname} eq "U" &&
		    ! exists($chr_pos_undef_dict{$snpname})) {
		    push(@chr_pos_undef, $snpname);
		    $chr_pos_undef_dict{$snpname} = 1;
		}
	    }
	    	    
	    if ($poscol > -1) {
		# Test if position is numeric or NA
		if ($word[$poscol - 1] eq "NA") {
		    $user_positions{$snpname} = $word[$poscol-1];
		}
		elsif ($word[$poscol - 1] eq "-99.00") {
		    $user_positions{$snpname} = "NA";
		}
		elsif (($word[$poscol - 1] + 0) eq $word[$poscol - 1]) {
		    $user_positions{$snpname} = $word[$poscol-1];
		}
		else {
		    $user_positions{$snpname} = "NA";
		    if (! exists($chr_pos_undef_dict{$snpname})) {
			push(@chr_pos_undef, $snpname);
			$chr_pos_undef_dict{$snpname} = 1;
		    }
		}		    
	    }
	    
	} elsif ($seen{$snpname} == 1) {
	    # We've seen this one before and reported
	    $seen{$snpname}++;
	    
	} else {
	    # Second time, so report the duplicate
	    $seen{$snpname} = 1;
	    push(@duplicates, $snpname);
	}	
    }
	
}

# If there are duplicate SNP ids in user's input file, list these
# inside duplicates file.

$dupl = $output_file."_duplicates.txt";
unlink($dupl);

if (@duplicates > 0) { 
    print("WARNING: Found duplicate SNPs in input file, \n         see ". 
	   "$dupl for complete list.\n");
    open(DH, ">$dupl") || die("Unable to open duplicate SNPs file $dupl");

    flock(DH, LOCK_EX);
    foreach $dupl (@duplicates) {
	print DH $dupl, "\n";
    }

    flock(DH, LOCK_UN);
    close(FH);
}
# If SNPs were renamed from RS to rs, log them now
if (scalar @found_RS > 0) {

    print("WARNING: Found uppercase rs IDs in input file, \n         see ". 
	  "$wfile for complete list.\n");

    open(WAR, ">>$wfile") || die("Cannot open $wfile");

    printf WAR "SNPs with uppercase rs IDs that were changed to lowercase\n";
    printf WAR "---------------------------------------------------\n\n";
    foreach $snp (@found_RS) {
	print WAR $snp, " => ", $lowered{$snp}, "\n";
    }
    printf WAR "---------------------------------------------------\n\n";
    close(WAR);
}
    

# If there were SNPs with unrecognized chromosomes and positions
if (scalar @chr_pos_undef > 0) {
    open(WAR, ">>$wfile") || die("Cannot open $wfile");
    # If every snp is undefined, then don't pass in a 3 column snplist
    if (scalar(@chr_pos_undef) == scalar(@snpid)) {
	print WAR "All SNPs have unrecognized chromosome numbers and/or positions\n";
	print WAR "chromosomes and positions will be ignored.\n";
	printf WAR "---------------------------------------------------\n\n";
	$chrcol = -1;
	$poscol = -1;
    }
    else {
	printf WAR "SNPs with unrecognized chromosome numbers and/or positions\n";
	printf WAR "(chromosome numbers set to \"U\" and positions set to \"NA\")\n";
	printf WAR "---------------------------------------------------\n\n";
	foreach $snp (@chr_pos_undef) {
	    print WAR $snp, "\n";
	}
    }
    printf WAR "---------------------------------------------------\n\n";
    close(WAR);
}

# removed code that fetched SNP info from ensembl
my $user_supplied_snps = "snp_user";

open(UFH, ">$user_supplied_snps") ||
    die "Could not open temporary file $user_supplied_snps!";

flock(UFH, LOCK_EX);	       

foreach $snp (@snpid)
{    
    print UFH $snp; # print the SNP name
    if ($chrcol > -1) {
	print UFH "\t", $user_chromosomes{$snp};
    }
	
    if ($poscol > -1) {
	print UFH "\t", $user_positions{$snp};
    }
    print UFH "\n";
}

flock(UFH, LOCK_UN);	       
close(UFH);

########################## CREATE R BATCH ###########################

$skip_plots = 1 - $skip_plots;
print "Creating genetic maps\n";

my $rsour = "gmi.bioc.mapping.R";
my $rbatch = "batch.R";

open(BAT, ">$rbatch") || die "Could not create R batch file $rbatch!\n";
flock(BAT, LOCK_EX);	       

print BAT "source('$mpath/$rsour')", "\n";
print BAT 
    "snp.bioc.mapping('$user_supplied_snps', '$output_file', '$over', \
'$mpath/rutgers_merged_$over', $skip_plots, $new_snps, '$maintainer_email')\n";

flock(BAT, LOCK_UN);	       
close(BAT);

#####################################################################

############################ RUN R BATCH ############################
if (-f "$rbatch.out") {
    unlink("$rbatch.out");
}

system("R --slave < $rbatch | tee -a -i $rbatch.out");

#files created by batch.R
#  output_annot.txt
#  output_nonannot.txt
#  output.pdf
#  output_log.txt


######################## PRINT OUT WARNINGS ############################

# Check if the warnings file exists;
if (-f $wfile) {
    open(WAR, "$wfile") || die "can't open $wfile!\n";
    my @warn = <WAR>;
    close(WAR);
    
    $show=0; 
    $linenum=0;
    foreach $line (@warn)    {
	if ($line =~ "SUMMARY OF WARNINGS FROM INTERPOLATION") {
	    $show=1;
	    print $warn[$linenum-1];
	}
	if ($show)	{
	    print $line;
	}
	$linenum++;
    }

    if (scalar @duplicates > 0 || scalar @chr_pos_undef > 0) {
	open(WAR, ">>$wfile") || die "can't open $wfile!\n";
	foreach $filehandle (qw/STDOUT WAR/) {
	    print $filehandle "In addition, input file contained problem SNPs:\n";
	    if (scalar @found_RS > 0) {
		print $filehandle 
		    "  -> Uppercase rs IDs found (listed at the top of this file)\n";
		print $filehandle 
		    "     (IDs have been changed to lowercase IDs in output map files.)\n";
	    }
	    if (scalar @duplicates > 0) {
		print $filehandle 
		    "  -> Duplicate SNP names (listed in ".$output_file."_duplicates.txt),\n";
		print $filehandle 
		    "     (duplicate SNP names have been omitted from output map files.)\n";
	    }
	    if (scalar @chr_pos_undef > 0) {
		if (scalar @chr_pos_undef == scalar @snpid) {
		    print $filehandle 
			"  -> All SNPs have unrecognized chromosomes or physical positions\n";
		    print $filehandle 
			"    chromosome numbers and positions have been ignored.\n";
		}
		else {
		    print $filehandle 
			"  -> Unrecognized chromosomes or physical positions (see top of log file)\n";
		    print $filehandle 
			"    (chromosome numbers and positions have been set to unknown)\n";
		}
	    }
	    print $filehandle "---------------------------------------------------\n\n";
	}
	close(WAR);
    }
    print "Please look at log file for more detailed information!\n";
}

#####################################################################
print("\n");

if (-f $output_file."_annot.txt" && -f $output_file."_nonannot.txt") {
    print("GMI query successful.\n");
    print "The following files were created:\n";
    print "  ".$output_file."_annot.txt       Map file in Mega2's annotated format\n";
    print "  ".$output_file."_nonannot.txt    Map file in Mega2's unannotated format\n";

    if (-f $output_file.".pdf") {
	print "  $output_file.pdf             PDF plots\n";
    }

    if ( -f $wfile) {
	print "  $wfile         Run summary\n";
    }

    if ( -f $dupl) {
	print "  $dupl  Duplicate SNP names\n";
    }

} 

else {
    print("GMI query failed!\n");
}

#####################################################################

# print "Cleaning up\n";
##r unlink($user_supplied_snps);

# print "Completed!\n";


sub print_usage {

    print <<USAGE 
      GMI version 1.0	
      Usage: gmi.pl snp_list_file [additional args]

      Arguments:  --namecol [col], default 1
                  --poscol [col], default = -1, meaning absent
                  --chrcol [col], default = -1 meaning absent
                  --header (no arguments), first line will be skipped if specified
                  --outfile [output-file-prefix], default snp_mapped
                  --noplots (no arguments), plotting will be skipped if specified
		  --newsnps (no arguments), new SNP names will be used if specified

USAGE
;
}

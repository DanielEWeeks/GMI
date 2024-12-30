#!/usr/bin/env perl -w
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
# ===========================================================================

use File::Basename;
require File::Spec;
use Cwd;
require "./gmi_bioc_utils.pl";

$RUTGERS_MAP = "smooth_map_b36";
$RUTGERS_URL = "http://compgen.rutgers.edu/downloads/" . $RUTGERS_MAP . ".zip";

sub download_rutgers_maps ($)
{
    
    $wget_prog = shift;
    
    # exec command;
    if (-e "$RUTGERS_MAP.zip") {
	system("/bin/rm -f $RUTGERS_MAP.zip");
    }
    $exec_command = $wget_prog . " " . $RUTGERS_URL;
    system($exec_command);

    system("unzip -t $RUTGERS_MAP.zip");
    if ($? == 0) {
	# rpath should already have been created
	system("unzip $RUTGERS_MAP.zip");
	if ($? == 0) {
	    system("/bin/rm -f $RUTGERS_MAP.zip");
	    return($rpath);
	}
	else {
	    print("Failed to unzip Rutgers maps.\n");
	    system("/bin/rm -rf $RUTGERS_MAP");
	    return "Failed";
	}
    }
    else {
	#clean up
	print("Corrupt or incomplete file $RUTGERS_MAP.zip, please run install again.\n");
	return "Failed";
    }
}

sub get_config ($)
    
{
    
    my $config_file = shift;    
    my %config_params = ();
    
    open(CONFIG, $config_file) || die("Could not open $config_file for reading.");
    
  LINE: while (<CONFIG>) {
      
      if (/^#/) {
	  next LINE;
      }	
      @items = split(/=/);
      if ($#items >= 1) {
	  $param=lc(trim($items[0]));
	  $value=lc(trim($items[1]));
	  $config_params{$param}=$value;
      }
      elsif ($#items == 0) {
	  $param=lc(trim($items[0]));
	  $config_params{$param}="";
      }
  }
    close(CONFIG);
    
    return %config_params;
}

sub print_config

{
    my %config_params = @_;
    
    open(CONF, ">CONFIG.txt") ||
	die("Could not open CONFIG.txt for writing");
    
    my @keys = ("rutgers-map-enclosing-folder", 
		"applications-folder", 
		"htdocs-folder",
		"setup-merged-maps",
		"install-perl-modules",
		"maintainer-email");
    
    foreach $param (@keys) {	
	if (exists($config_params{$param})) {
	    print CONF $param, " = ", $config_params{$param}, "\n";
	}
    }
    
    close(CONF);
    
    return;
}

sub print_log 
    
{
    my $log = shift;
    my $message = shift;
    
    print $log $message, "\n";
    print $message, "\n";
}


my %config_params=();

if (! $ARGV[0]) {
    $config_file = "";
}
else {
    $config_file=$ARGV[0];
    %config_params = &get_config($config_file);
}


######################################################################

# Check for necessary Perl modules

@required = ("Cwd", "Fcntl", "warnings", "Getopt::Long", 
	     "File::Basename", "File::Copy");

@install_perl = ();

foreach $module (@required) {		 
    if (system("perl -M".$module." -e 1 2>&1") > 0) {
	push(@install_perl, $module);
    }
    else {
	print "Found module $module\n";
    }
}

if (scalar(@install_perl) > 0 ) {  
    print("Missing Perl modules " . join(", ", @install_perl) . "\n");
    if (exists($config_params{"install-perl-modules"})) {
	$setup = $config_params{"install-perl-modules"};
	($setup eq "no") ? 
	    print("Not installing Perl modules.\n") :
	    print("Installing Perl modules ... \n");
    }
    else {
	print("\n");
	print "Call cpan to install missing Perl module(s)? [yes/no](default yes) > ";
	$setup = <STDIN>;
	chomp($setup);
	$config_params{"install-perl-modules"} = $setup;
    }
    
    if ($setup eq "") {
	$setup = "yes";
    }
    
    if ($setup ne "no") {
	system("/usr/bin/env cpan -fi " . join(" ", @install_perl));
	print("\n");
	
	print("Installed " . join(" ", @install_perl) . " in your Perl libraries\n.");
    }
    else {
	print "Missing Perl modules need to be installed first.\n";
	print "Terminating GMI install.\n";
	exit(-1);
    }
}

$current_path = getcwd();

if (! -w $current_path) {
    print("Error: You do not have write permissions on this folder.\n");
    die("Installation requires such permissions, terminating GMI installation.\n");
}

if (-f "install_log.txt" && ! -w "install_log.txt") {
    print("Error: Unable to overwrite install_log.txt.\n");
    print("Terminating GMI installation.\n");
    exit(-1);
}

open(LOG, ">install_log.txt") || 
    die("Could not open install_log.txt for writing.");
my $log=LOG;

print "Configuring GMI ...\n";


$wget_prog = `which ftp`;
chomp($wget_prog);

if ($wget_prog eq "") {
    print_log($log, "Ftp not found, looking for wget ...");    
    $wget_prog = `which wget`;
    chomp($wget_prog);
    if ($wget_prog eq "") {
	die("Could not find ftp or wget, please install one of these programs first!");
    }
}
else {
    print_log($log, "Found $wget_prog.");
}

print_log($log, "Looking for R ...");
$r_prog = `which R`;
chomp($r_prog);

if ($r_prog eq "") {
    die("Could not find R, please install R first!");
}
else {
    print_log($log, "Found $r_prog.");
}

print_log($log, " ");

$rpath = $current_path;

if ( ! -d $RUTGERS_MAP ) {
    print("Download Rutgers maps now (necessary for installing GMI)? (yes/no) [yes] > ");
    $create = <STDIN>; chomp($create);
    if ($create ne "no") {	
	$installed_rutgers = &download_rutgers_maps($wget_prog);
	if ($installed_rutgers eq "Failed") {
	    print_log($log, 
		      "ERROR: Failed to create Rutgers combined map folder.");
	    die("Terminating GMI install.\n");
	}       
    }
    else {
	print_log($log, 
		  "ERROR: Cannot proceed with installing GMI (please download Rutgers maps first).");
	die("Terminating GMI install.\n");
    }
}
else {
    print("Rutgers map folder already exists, re-install now? (yes/no) [yes] > ");
    $reinstall = <STDIN>; chomp($reinstall);
    if ($reinstall ne "no") {
	if (! -w $RUTGERS_MAP) {
	    print_log($log, 
		      "Error: do not have write permissions to $RUTGERS_MAP.");
	    die("Terminating GMI install.\n");
	}
	system("/bin/rm -rf $RUTGERS_MAP $RUTGERS_MAP.backup");
	$installed_rutgers = &download_rutgers_maps($wget_prog);
	if ($installed_rutgers eq "Failed") {
	    print_log($log, 
		      "ERROR: Failed to create Rutgers combined map folder.");
	    system("/bin/mv $RUTGERS_MAP.backup $RUTGERS_MAP");
	    print_log($log, "Restored old Rutgers files.\n");
	    die("Terminating GMI install.\n");
	}       
	else {
	    print_log($log, "Successfully downloaded Rutgers maps.");
	    system("/bin/rm -rf $RUTGERS_MAP.backup");
	}
    }
    else {
	print_log($log, "Rutgers maps will not be re-installed.\n");		  
	$installed_rutgers = "not re-installed";
    }
}

print_log($log, " ");

if (exists($config_params{"applications-folder"})) {
    $epath = $config_params{"applications-folder"};
}
else {
    print "Folder in which to install gmi.pl [/usr/local/bin] > ";
    $epath = <STDIN>;

    chomp($epath);
}

if ($epath eq "") {
    $epath = "/usr/local/bin";
}

else {
    $epath = File::Spec->rel2abs($epath) ;
}
    
if ( ! -d  $epath ) {
    die("ERROR: Specified folder for installing scripts does not exist or is unreadable.\n");
}
else {
    print_log($log, "Found executables folder: " . $epath);
    $config_params{"applications-folder"} = $epath;
}

if (exists($config_params{"maintainer-email"})) {
    $maintainer_email = $config_params{"maintainer-email"};
}
else {
    print "Your e-mail address (necessary for Entrez queries) > ";
    $maintainer_email = <STDIN>;
    chomp($maintainer_email);
}

if ($maintainer_email eq "") {
    print_log($log, "ERROR: Invalid e-mail address.");
    die("Terminating GMI install.\n");
}
else {
    $config_params{"maintainer-email"} = $maintainer_email;
}

$maintainer_email =~ s/\@/\\\@/ ;

print_log($log, " ");

######################################################################
# Set up path for mapping
######################################################################

# Set the mpath

$mpath = $current_path;
$perl_path = `which perl`;
chomp($perl_path);

open(QUE, "gmi_bioc_query.pl") ||
    die("Could not open gmi_bioc_query.pl for reading.");
my @oque = <QUE>;
close(QUE);

open(NQUE, ">gmi_bioc_query_new.pl") ||
    die("Could not open temporary file gmi_bioc_query_new.pl for writing.");
print NQUE "#!".$perl_path."\n";

# get rid of the first block of comments

$line = shift(@oque);

while ($line !~ /insert/) {
    print NQUE $line;
    $line = shift(@oque);
}

print NQUE $line;

#1) set variable $mpath
print NQUE "\n# Path where merged files are located.\n";
print NQUE "\$mpath=\"".$mpath."\";\n";

#2) wget program
print NQUE "\n# Path to wget executable.\n";
print NQUE "\$wget_prog=\"".$wget_prog."\";\n";

#3) version number
print NQUE "\n# GMI version number.\n";
print NQUE "\$gmi_version=\"1.5\";\n";

#4) user's e-mail
print NQUE "\n# GMI user e-mail.\n";
print NQUE "\$maintainer_email=\"".$maintainer_email."\";\n";

#Rest of the script need not be changed
foreach $line (@oque)
{
    print NQUE $line;
}

close(NQUE);

# Move gmi_bioc_query_new.pl to gmi_bioc_query.pl
# Not a good idea to overwrite, if we are not going to install
# Check for argv[1] is adequate

if (-w $epath) {
	# Move gmi_bioc_query.pl inside environment path;

    if (-f "$epath/gmi.pl") {
	print_log($log, "Existing gmi.pl will be overwritten.");
    }
    system("mv gmi_bioc_query_new.pl $epath/gmi.pl");
    system("chmod +x $epath/gmi.pl");

    $installed_query = "$epath/gmi.pl";
    print_log($log, "gmi.pl written to $epath.");
}
else {
    print_log($log, "ERROR: Cannot write gmi.pl in $epath, no write permissions.");
    $installed_query = "not installed";
    system("mv gmi_bioc_query_new.pl gmi.pl");
    system("chmod +x gmi.pl");    
}

print_log($log, " ");

######################################################################
# Removed Perl biomart code, 
# use Bioc setdata function instead.
######################################################################

# Set the directory path

my $path_raw = "$mpath/smooth_map_b36";

print "Check if Rutgers map files are in the right position\n";

open(CHE, "$path_raw/chr1.sm.map2") || die "can't locate Rutgers map files!\n";
close(CHE);

######################################################################

# Call R to run the setdata function

print_log($log, "Creating R batch file gmi.bioc.setdata.batch.R");

my $rbatch = "$mpath/gmi.bioc.setdata.batch.R";

my $rsour = "gmi.bioc.setdata.R";
open(BAT, ">$rbatch") || die("Could not create R batch file $rbatch");
print BAT "source('$mpath/$rsour')\n";
print BAT "status <- gmi.bioc.setdata('$mpath','$path_raw')\n";
print BAT "write(file='install_status.txt', status)\n";
print BAT "q()\n";
close(BAT);

if ($config_file eq "") {
    print_config(%config_params);
    print_log($log, "Wrote configuration parameters to CONFIG.txt");

}

print_log($log, "Running $rbatch, this may take a while ...");
close(LOG);

system("R --slave < $rbatch | tee -i -a install_log.txt");

open(IL, "install_status.txt") || die("Could not open install_status.txt");
$setdata_status = <IL>;
chomp($setdata_status);
$setdata_status = 1 * $setdata_status;
close(IL);

open(LOG, ">>install_log.txt") || 
    die("Could not append to install_log.txt.");
$log=LOG;

print_log($log, " ");
if ($setdata_status <= 0) {
     print_log($log, "Failed to create merged maps.");
     print_log($log, "Please try again.\n");
     $installed_merged = "not installed";
}
else {
    $installed_merged = "rutgers_merged_$setdata_status";
}

print_log($log,
	  "To re-configure, run gmi_bioc_install.pl again from within this folder.");
print_log($log, " ");
print_log($log, "Installed the following GMI components:");
print_log($log, "  GMI query program:   $installed_query");
print_log($log, "  Rutgers map files:   $installed_rutgers");
print_log($log, "  Merged map files:    $installed_merged");
print_log($log, "  Configuration file:  $mpath/CONFIG.txt");

if ($installed_query eq "not installed") {
    print_log($log, " ");
    print_log($log, "REMINDER: gmi.pl was not installed in $epath.");
    print_log($log, "          To run GMI, please copy gmi.pl to $epath.");
}
else {
    print_log($log, "To use GMI, run gmi.pl.");
}

print_log($log, "\nFinished installing GMI.");

# }
# print  "Run $rbatch now? [yes/no](default yes) > ";
# $ans=<>;
# chomp($ans);
# if ($ans ne "no") {
#     print_log($log, "Running $rbatch, this may take a while ...");
#     system("R CMD BATCH --no-save --no-restore --no-readline $rbatch $rbatch.out");
#     print_log($log, "Finished setting up GMI data files.");
# }
# else {
# }

# print_log($log,"Please remember to run $rbatch using R before using GMI query.\n");


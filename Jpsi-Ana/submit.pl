#!/usr/bin/perl -w

use File::Copy;
use Cwd;

open (FILE, "<list_pp_allfile.lst") || die "File Open Error!!\n";
@fullpathlist = <FILE>;
close(FILE);

open (FILE, "<list_pp_runnum.lst") || die "File Open Error!!\n";
@runlist = <FILE>;
close(FILE);

$maindir = getcwd();
$jobdir = "${maindir}/jobdir_Run15pp200_02";

mkdir $jobdir;

foreach $runnum (@runlist){

	chomp($runnum);
	
	$wrkdir = "${jobdir}/run${runnum}";
	mkdir $wrkdir;
	chdir $wrkdir;

	open (FILE, ">file.lst");
	foreach $filename (@fullpathlist){
		if ( $filename =~ /$runnum/ ){
			print FILE $filename;
		}
	}
	close (FILE);

	open (FILE, ">condor");
	print FILE "Universe      = vanilla\n";
	print FILE "Notification  = Never\n";
#	print FILE "Requirements  = CPU_Speed >= 2\n";
#	print FILE "Rank          = CPU_Speed\n";
	print FILE "Getenv        = true\n";
	print FILE "Priority      = +50\n";
	print FILE "Executable    = jobscript\n";
	print FILE "Output        = jobscript.out\n";
	print FILE "Error         = jobscript.err\n";
	print FILE "Log           = jobscript.log\n";
#	print FILE "Notify_user   = shlim\@rcf.rhic.bnl.gov\n";
#	print FILE "+Experiment   = \"phenix\"\n";
#	print FILE "+Job_Type     = \"cas\"\n";
	print FILE "Queue";
	close (FILE);

	system "cp ${maindir}/ScanEventTree.C .";

	open (FILE, ">jobscript");
	print FILE "#!/bin/bash -f\n";
#	print FILE "source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/setup.csh\n\n";
#	print FILE "setenv ODBCINI /opt/phenix/etc/odbc.ini.phnxdbrcf2\n";

	$outputname = "${jobdir}/Run15pp200_ScanEventTree_${runnum}.root";

	print FILE "root -l -b -q 'ScanEventTree.C+g'\n";

	print FILE "mv outfileScanEventTree.root ${outputname}\n";
	print FILE "rm -rf *.root\n";

	close (FILE);
	chmod 0755, "jobscript";

	system "condor_submit condor";
}

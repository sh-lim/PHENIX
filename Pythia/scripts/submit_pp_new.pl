#!/usr/bin/perl -w

use Cwd;

$package = "pp13TeV_set00";
$maindir = getcwd();

$groupnum = 0;

$rundir = "${maindir}/run_${package}_grp${groupnum}";
mkdir $rundir;

for ($irun=0; $irun<100; $irun++){

#	sleep 1;

	$wrkdir = "${rundir}/wrk_${irun}";
	mkdir $wrkdir;

	chdir $wrkdir;
	open(FILE, ">condor");
	print FILE "Universe = vanilla\n";
	print FILE "Notification = Never\n";
	print FILE "Requirements = CPU_Speed>=1\n";
	print FILE "Rank = CPU_Speed\n";
	print FILE "Getenv = true\n";
	print FILE "Priority = +20\n";
	print FILE "Executable = jobscript\n";
	print FILE "Log = jobscript.log\n";
	print FILE "Output = jobscript.out\n";
	print FILE "Error = jobscript.err\n";
	print FILE "Notify_user = shlim\@rcf.rhic.bnl.gov\n";
#	print FILE "+Experiment = \"phenix\"\n";
#	print FILE "+Job_Type = \"cas\"\n";
	print FILE "Queue\n";
	close(FILE);

#	$seednum = int(rand(100000000));

#	$grpdir = sprintf("%s/%s_grp%03d",$maindir,$package,$groupnum);
#	$outputdir = sprintf("%s/%s_grp%03d","/gpfs/mnt/gpfs02/sphenix/user/shlim/99.Tmp/02.ClosureSample",$package,$groupnum);
	$outputdir = sprintf("%s/%s_grp%03d","/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/shlim/simulation/ClosureSample",$package,$groupnum);

	$randnum = int(rand(60));

	open(FILE, ">jobscript");
	print FILE "#!/bin/csh -f\n";
#	print FILE "sleep ${randnum}; source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/setup.csh\n\n\n";
	print FILE "sleep ${randnum}; source /opt/phenix/bin/phenix_setup.csh\n\n\n";

#	print FILE "mkdir -p $grpdir\n\n";
	print FILE "mkdir -p $outputdir\n\n";

#	$condor_wrk_dir = "grp${groupnum}_run${irun}";
	$condor_wrk_dir = "grp${groupnum}_run${irun}";
	print FILE "mkdir -p \$_CONDOR_SCRATCH_DIR/${condor_wrk_dir}\n";
	print FILE "cd \$_CONDOR_SCRATCH_DIR/${condor_wrk_dir}\n";

	print FILE "cp -av ${maindir}/codes/mainEx00.cc .\n";
	print FILE "cp -av ${maindir}/codes/mainEx00.cfg .\n";
	print FILE "cp -av ${maindir}/codes/Makefile .\n";
	print FILE "cp -av ${maindir}/codes/Makefile.inc .\n\n";

	print FILE "make mainEx00\n";
	print FILE "./mainEx00\n";

#	$outputfile = sprintf("%s/outfile_%s_grp%03d_%05d.root",$grpdir,$package,$groupnum,$irun);
	$outputfile = sprintf("%s/outfile_%s_grp%03d_%05d.root",$outputdir,$package,$groupnum,$irun);
	print FILE "mv Pythia8_event.root $outputfile\n";

	print FILE "rm -rf mainEx00\n\n";

	close(FILE);
	chmod 0755, "jobscript";
	system "condor_submit condor";
}


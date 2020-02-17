#!/usr/bin/perl -w

use Cwd;

#$package = "dAu200GeV";
#$package = "pAu200GeV";
#$package = "pp13TeV";
#$package = "pp13TeV_PYTHIA";
$package = "pp13TeV_nondiffractive";
#$package = "pp200GeV";
#$package = "pp200GeV_nondiffractive";
#$package = "pp200GeV_all";
#$package = "pp200GeV_inelastic";
#$package = "pPb8000GeV";
#$package = "pPb8160GeV_PYTHIA";
#$package = "pPb200GeV_PYTHIA";
#$package = "pAu200GeV_PYTHIA";
$maindir = getcwd();

$filedir = "/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/shlim/simulation/ClosureSample";
#$filedir = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/02.PythiaClosure/sHijing";
#$filedir = "/gpfs/mnt/gpfs02/sphenix/user/shlim/99.Tmp/02.ClosureSample";

$groupnum = 106;
$try = 1000;

#$files_per_job = 50; #pAu200GeV hijing, 10000 files
#$files_per_job = 20; #pp200GeV hijing, 2000 files
#$files_per_job = 10; #pp200GeV pythia, 1000 files
#$files_per_job = 20; #pp200GeV pythia softQCD:all, 2000 files
#$files_per_job = 20; #pp200GeV pythia softQCD:inelastic, 2000 files
#$files_per_job = 10; #pp13TeV pythia, 5000 files
#$files_per_job = 10; #pp13TeV hijing, 5000 files
#$files_per_job = 5; #pPb8000GeV hijing, 5000 files
#$files_per_job = 10; #pPb8000GeV hijing jetpt15, 5000 files
#$files_per_job = 15; #pPb8000GeV hijing jetpt30, 15000 files
#$files_per_job = 40; #pPb8000GeV hijing jetpt50, 8000 files
#$files_per_job = 10; #pPb8160GeV pythia pPb, 8000 files
#$files_per_job = 40; #pPb8160GeV pythia pPb jetpt30 R04, 10000 files
#$files_per_job = 30; #pPb8160GeV pythia pPb jetpt30 R07, 6000 files
#$files_per_job = 200; #pPb8160GeV pythia pPb jetpt50 R04, 50000 files
#$files_per_job = 20; #pp13TeV Pythia pp, 2000 files
#$files_per_job = 40; #pp13TeV Pythia pp, 4000 files
#$files_per_job = 2000; #pp13TeV Pythia pp for HFmu, 10000 files per group, group 102-106
#$files_per_job = 10; #pp13TeV pythia (string shoving), 1000 files
#$files_per_job = 30; #pPb200GeV pythia, 5000 files
$files_per_job = 50; #pp13TeV pythia (string shoving), 1000 files

$rundir = "${maindir}/runhist_${package}_grp${groupnum}_try${try}";
mkdir $rundir;

for ($irun=0; $irun<20; $irun++){

	$wrkdir = "${rundir}/wrk_${irun}";
	mkdir $wrkdir;

	chdir $wrkdir;
	open(FILE, ">condor");
	print FILE "Universe = vanilla\n";
	print FILE "Notification = Never\n";
	print FILE "Requirements = CPU_Speed>=2\n";
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

	open(FILE, ">jobscript");
	print FILE "#!/bin/csh -f\n";
	print FILE "source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/setup.csh\n\n\n";

	print FILE "cp -av ${maindir}/Run.C .\n";
	print FILE "cp -av ${maindir}/make_hist_dAu200GeV.C .\n";
	print FILE "cp -av ${maindir}/make_hist_dAu200GeV_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pAu200GeV.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pAu200GeV_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pAu200GeV_pT_phenix.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pPb8TeV_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pPb8TeV_jet_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp200GeV.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp200GeV_4pc_mult.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp200GeV_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp200GeV_all_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp200GeV_inelastic_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp200GeV_hijing.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp200GeV_hijing_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_ineff.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_jet_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_jetonly_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_hijing_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_4pc_pT.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_4pc_mult.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_jet_4pc_mult.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pPb8TeV_pT_pythia.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pPb8TeV_pT_pythia_jet.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pPb8TeV_pT_hijing.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_HFmu.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pPb200GeV_pT_pythia.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pAu200GeV_pT_pythia.C .\n";
	print FILE "cp -av ${maindir}/make_hist_pp13TeV_alice.C .\n";

#	print FILE "root -l -b -q make_hist_pp13TeV_HFmu.C+g\n";
#	print FILE "root -l -b -q make_hist_pp13TeV_ineff.C+g\n";
#	print FILE "root -l -b -q make_hist_pp13TeV.C+g\n";
#	print FILE "root -l -b -q make_hist_pp13TeV_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pp13TeV_jet_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pp13TeV_jetonly_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pp13TeV_hijing_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pp13TeV_4pc_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pp13TeV_4pc_mult.C+g\n";
#	print FILE "root -l -b -q 'make_hist_pp13TeV_jet_4pc_mult.C+g(\"file.lst\",50.0)'\n";
#	print FILE "root -l -b -q make_hist_pp200GeV.C+g\n";
#	print FILE "root -l -b -q make_hist_pp200GeV_4pc_mult.C+g\n";
#	print FILE "root -l -b -q make_hist_pp200GeV_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pp200GeV_all_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pp200GeV_inelastic_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pAu200GeV.C+g\n";
#	print FILE "root -l -b -q make_hist_pAu200GeV_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pAu200GeV_pT_phenix.C+g\n";
#	print FILE "root -l -b -q make_hist_pp200GeV_hijing.C+g\n";
#	print FILE "root -l -b -q make_hist_pp200GeV_hijing_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pPb8TeV_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pPb8TeV_jet_pT.C+g\n";
#	print FILE "root -l -b -q make_hist_pPb8TeV_pT_pythia.C+g\n";
#	print FILE "root -l -b -q make_hist_pPb8TeV_pT_pythia_jet.C+g\n";
#	print FILE "root -l -b -q make_hist_pPb8TeV_pT_hijing.C+g\n";
#	print FILE "root -l -b -q make_hist_pPb200GeV_pT_pythia.C+g\n";
#	print FILE "root -l -b -q make_hist_pAu200GeV_pT_pythia.C+g\n";
	print FILE "root -l -b -q 'make_hist_pp13TeV_alice.C+g(\"file.lst\",${groupnum})'\n";

	$filename = sprintf("%s/outfile_%s_grp%03d_%05d.root",$rundir,$package,$groupnum,$irun);
	print FILE "mv -v outfile_hist.root $filename\n\n";

	print FILE "rm -rf *.C\n\n";

	close(FILE);
	chmod 0755, "jobscript";

	$grpdir = sprintf("%s/%s_grp%03d",$filedir,$package,$groupnum);

	open(FILE,">file.lst");
	for ($iseg=0; $iseg<$files_per_job; $iseg++){
		$filename = sprintf("%s/outfile_%s_grp%03d_%05d.root",$grpdir,$package,$groupnum,$files_per_job*$irun+$iseg);
		if ( -e $filename ){
			print FILE $filename."\n";
		}
	}
	close(FILE);

	system "condor_submit condor";

}

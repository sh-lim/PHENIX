#!/usr/bin/perl -w

use File::Copy;
use Cwd;

#open (FILE, "<list_Run15pAu200MU_pro105_picoDST_dimu_fullpath_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200MU_pro105_pro107_picoDST_dimu_fullpath_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200MU_pro105_pro108_picoDST_dimu_fullpath_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pp200MU_pro108_picoDST_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pp200MU_pro108_picoDST_dimu_fullpath_100.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pp200MU_pro108_picoDST_dimu_fullpath_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pp200MB_ERT_MPC_pro108_picoDST_fullpath.lst") || die "File Open Error!!\n";
open (FILE, "<list_Run15pAu200MU_pro105_noVTX_pro108_picoDST_dimu_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run14AuAu200MB_pro109_picoDST_dimu_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run14HeAu200MU_pro109_picoDST_dimu_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run14HeAu200MU_pro109_picoDST_dimu_fullpath_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run14HeAu200MB_pro109_picoDST_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAl200MU_pro105_picoDST_dimu_fullpath_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAl200MU_pro105_picoDST_dimu_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAl200ALL_pro105_picoDST_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAl200OT_pro105_picoDST_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAl200MU_pro105_picoDST_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200MU_pro108_picoDST_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200MU_pro105_pro107_picoDST_fullpath.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run16dAu200MB_pro107_picoDST_fullpath.lst") || die "File Open Error!!\n";
@fullpathlist = <FILE>;
close(FILE);

#open (FILE, "<list_Run16dAu200MB_pro107_picoDST_runnum.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAl200MU_pro105_picoDST_runnum.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAl200MU_pro105_picoDST_dimu_runnum_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAl200MU_pro105_picoDST_dimu_runnum_Matt.lst") || die "File Open Error!!\n";
#
#open (FILE, "<list_Run15pp200MU_pro108_picoDST_runnum.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pp200MB_ERT_MPC_pro108_picoDST_runnum.lst") || die "File Open Error!!\n";
#
#open (FILE, "<list_Run15pAu200MU_pro105_picoDST_dimu_runnum_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200MU_pro105_pro108_picoDST_dimu_runnum_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200MB_ERT_MPC_FVTX_pro105_pro108_picoDST_runnum.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200MU_pro105_noVTX_pro108_picoDST_runnum.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200ALL_pro108_picoDST_runnum.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run15pAu200MU_pro108_picoDST_runnum.lst") || die "File Open Error!!\n";
open (FILE, "<list_Run15pAu200MU_pro105_noVTX_pro108_picoDST_runnum.lst") || die "File Open Error!!\n";
#
#open (FILE, "<list_Run14HeAu200MU_pro103_picoDST_dimu_runnum_Matt.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run14HeAu200MB_pro109_picoDST_runnum.lst") || die "File Open Error!!\n";
#open (FILE, "<list_Run14AuAu200MB_pro109_picoDST_dimu_runnum.lst") || die "File Open Error!!\n";
#open (FILE, "<missing.lst") || die "File Open Error!!\n";
@runlist = <FILE>;
close(FILE);

$maindir = getcwd();
#$jobdir = "${maindir}/jobdir_Run15pp200MU_pro108_sngmu_00";
#$jobdir = "${maindir}/jobdir_Run15pp200MU_pro108_post_aligned_swap_sngmu_20th";
#$jobdir = "${maindir}/jobdir_Run15pp200MB_pro106_post_aligned_swap_sngmu_1st";
#$jobdir = "${maindir}/jobdir_Run15pp200MU_pro106_dimu_swap_10th";
#$jobdir = "${maindir}/jobdir_Run15pp200MU_pro106_dimu_mixed_bbcz";
#$jobdir = "${maindir}/jobdir_Run15pp200MU_pro108_FVTX_track";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro108_FVTX_track";
#$jobdir = "${maindir}/jobdir_Run15pp200MU_pro106_inclHF";
#$jobdir = "${maindir}/jobdir_Run15pp200MU_pro108_hadron_FVTX_3rd";
#$jobdir = "${maindir}/jobdir_Run15pp200MB_ERT_MPC_pro108_hadron_FVTX_00";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro107_sngmu_post_aligned_swap_5th";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro107_dimu_swap_4th";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro107_dimu_mixed_bbcz";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro107_FVTX_track";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro108_hadron_FVTX_00";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro105_pro107_inclHF";
#$jobdir = "${maindir}/jobdir_Run15pAl200MU_OT_pro105_FVTX_track";
#$jobdir = "${maindir}/jobdir_Run15pAl200OT_pro105_hadron_FVTX_00";
#$jobdir = "${maindir}/jobdir_Run15pAl200MU_pro105_hadron_FVTX_00";
#$jobdir = "${maindir}/jobdir_Run15pAl200OT_pro105_sngmu_01";
#$jobdir = "${maindir}/jobdir_Run15pAl200MU_pro105_sngmu_01";
#$jobdir = "${maindir}/jobdir_Run15pAl200MB_ERT_FVTX_pro105_sngmu_00";
#$jobdir = "${maindir}/jobdir_Run15pAu200MB_ERT_MPC_FVTX_pro105_pro108_sngmu_00";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro108_sngmu_01";
#$jobdir = "${maindir}/jobdir_Run14HeAu200MB_pro109_sngmu_00";
#$jobdir = "${maindir}/jobdir_Run14HeAu200MB_pro109_hadron_FVTX_00";
#$jobdir = "${maindir}/jobdir_Run16dAu200MB_pro107_hadron_FVTX_01";
#$jobdir = "${maindir}/jobdir_Run15pp200MU_pro108_inclusive_jpsi_100";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro105_pro107_inclusive_jpsi_01";
$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro105_pro108_inclusive_jpsi_100";
#$jobdir = "${maindir}/jobdir_Run14HeAu200MU_pro109_inclusive_jpsi_05";
#$jobdir = "${maindir}/jobdir_Run15pAl200MU_pro105_inclusive_jpsi_100";
#$jobdir = "${maindir}/jobdir_Run14AuAu200MB_pro109_dimuon_02";
#$jobdir = "${maindir}/jobdir_Run15pAu200MU_pro105_dimuon_02";

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
#			$thisname = $filename;
		}
	}
	close (FILE);

	open (FILE, ">condor");
	print FILE "Universe      = vanilla\n";
	print FILE "Notification  = Never\n";
	print FILE "Requirements  = CPU_Speed >= 2\n";
	print FILE "Rank          = CPU_Speed\n";
	print FILE "Getenv        = true\n";
	print FILE "Priority      = +50\n";
	print FILE "Executable    = jobscript\n";
	print FILE "Output        = jobscript.out\n";
	print FILE "Error         = jobscript.err\n";
	print FILE "Log           = jobscript.log\n";
	print FILE "Notify_user   = shlim\@rcf.rhic.bnl.gov\n";
#	print FILE "+Experiment   = \"phenix\"\n";
#	print FILE "+Job_Type     = \"cas\"\n";
	print FILE "Queue";
	close (FILE);

	system "cp ${maindir}/Run_fill_ana_tree.C .";
	system "cp ${maindir}/Run_fill_ana_tree_dimu.C .";
	system "cp ${maindir}/Run_fill_fvtx_tracklet.C .";

	open (FILE, ">jobscript");
	print FILE "#!/bin/csh -f\n";
	print FILE "source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/setup.csh\n\n";
#	print FILE "setenv ODBCINI /opt/phenix/etc/odbc.ini.phnxdbrcf2\n";
#	print FILE "ln -sf /direct/phenix+hhj/shlim/work/15.run15/10.runQA/09.Golden_dimu/Run_fill_ana_tree.C\n";
#	print FILE "ln -sf /direct/phenix+hhj/shlim/work/15.run15/10.runQA/09.Golden_dimu/make_picodstobj.C\n";
#	print FILE "ln -sf /direct/phenix+hhj/shlim/work/15.run15/09.re_association/make_picodstobj_dimu.C\n";
#	print FILE "ln -sf /direct/phenix+hhj/shlim/work/15.run15/09.re_association/make_picodstobj.C\n";

#	$outputname = "${jobdir}/${dataset}${trig}_pro105_beamXY_hist_${runnum}.root";
#	$outputname = "${jobdir}/Run15pp200_pro108_femtoDST_dimu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pp200_pro108_femtoDST_sngmu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAu200MU_pro107_femtoDST_dimu_mixed_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAu200MU_pro107_femtoDST_dimu_swapped_${runnum}.root";
#	$outputname = "${jobdir}/Run15pp200MB_ERT_MPC_pro108_femtoDST_sngmu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAu200_pro108_femtoDST_sngmu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAu200MU_pro105_pro107_femtoDST_sngmu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAu200_pro105_pro108_femtoDST_sngmu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAl200_pro105_femtoDST_sngmu_${runnum}.root";
#	$outputname = "${jobdir}/Run14HeAu200_pro109_femtoDST_sngmu_${runnum}.root";
#	$outputname = "${jobdir}/Run16dAu200_pro107_femtoDST_sngmu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAu200MU_pro105_pro107_femtoDST_dimu_${runnum}.root";
	$outputname = "${jobdir}/Run15pAu200MU_pro105_pro108_femtoDST_dimu_${runnum}.root";
#	$outputname = "${jobdir}/Run14HeAu200_pro109_femtoDST_dimu_${runnum}.root";
#	$outputname = "${jobdir}/Run14AuAu200_pro109_femtoDST_dimu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAl200_pro105_femtoDST_dimu_${runnum}.root";
#	$outputname = "${jobdir}/Run15pAu200_pro105_femtoDST_dimu_${runnum}.root";

#	print FILE "root -l -b -q 'make_picodstobj.C(0,\"file.lst\")'\n";
#	print FILE "root -l -b -q 'make_picodstobj_dimu.C(0,\"file.lst\")'\n";
#	print FILE "ls picodstobj_local.root > file.lst\n";

#	print FILE "root -l -b -q 'Run_fill_ana_tree.C(\"file.lst\",\"femtoDST.root\",0)'\n";
	print FILE "root -l -b -q 'Run_fill_ana_tree_dimu.C(\"file.lst\",\"femtoDST.root\",0)'\n";
#	print FILE "root -l -b -q 'Run_fill_fvtx_tracklet.C(\"femtoDST.root\",\"file.lst\",0)'\n";
	print FILE "mv femtoDST.root ${outputname}\n";
#	print FILE "mv histogram_out.root ${outputname}\n";
	print FILE "rm -rf *.root\n";

	close (FILE);
	chmod 0755, "jobscript";

	system "condor_submit condor";
}

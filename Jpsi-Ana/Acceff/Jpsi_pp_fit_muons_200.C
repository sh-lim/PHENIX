/// p+p -> J/psi-> mu mu data from PHENIX muon data at 200 GeV

//use this to figure out a weighting function for the pT dependence of the dimuon pairs I will use for the acc*eff correction

TH1D* hfit;

#include "/phenix/u/shlim/Style.h"

void Jpsi_pp_fit_muons_200()
{

	gRandom = new TRandom3(0);
  
  //data from PPG
  double pt[27], y[27], e1[27],e2[27],m[27];

  pt[0] =    0.125; y[0] =  3.48; e1[0] =  0.16; e2[0] =  0.23; m[0] = 1;
  pt[1] =    0.375; y[1] =  3.24; e1[1] =  0.10; e2[1] =  0.22; m[1] = 1;
  pt[2] =    0.625; y[2] =  2.82; e1[2] =  0.08; e2[2] =  0.19; m[2] = 1;
  pt[3] =    0.875; y[3] =  2.52; e1[3] =  0.06; e2[3] =  0.16; m[3] = 1;
  pt[4] =    1.125; y[4] =  2.07; e1[4] =  0.05; e2[4] =  0.13; m[4] = 1;
  pt[5] =    1.375; y[5] =  1.62; e1[5] =  0.04; e2[5] =  0.10; m[5] = 1;
  pt[6] =    1.625; y[6] =  1.20; e1[6] =  0.03; e2[6] =  0.08; m[6] = 1;
  pt[7] =    1.875; y[7] =  8.92; e1[7] =  0.21; e2[7] =  0.58; m[7] = 0.1;
  pt[8] =    2.125; y[8] =  6.50; e1[8] =  0.17; e2[8] =  0.42; m[8] = 0.1;
  pt[9] =    2.375; y[9] =  4.56; e1[9] =  0.13; e2[9] =  0.30; m[9] = 0.1;
  pt[10] =   2.625; y[10] =  3.16; e1[10] =  0.10; e2[10] =  0.21; m[10] = 0.1;
  pt[11] =   2.875; y[11] =  2.20; e1[11] =  0.08; e2[11] =  0.14; m[11] = 0.1;
  pt[12] =   3.125; y[12] =  1.51; e1[12] =  0.06; e2[12] =  0.10; m[12] = 0.1;
  pt[13] =   3.375; y[13] =  1.09; e1[13] =  0.05; e2[13] =  0.07; m[13] = 0.1;
  pt[14] =   3.625; y[14] =  7.71; e1[14] =  0.38; e2[14] =  0.50; m[14] = 0.01;
  pt[15] =   3.875; y[15] =  5.38; e1[15] =  0.31; e2[15] =  0.35; m[15] = 0.01;
  pt[16] =   4.125; y[16] =  3.30; e1[16] =  0.23; e2[16] =  0.22; m[16] = 0.01;
  pt[17] =   4.375; y[17] =  2.23; e1[17] =  0.18; e2[17] =  0.15; m[17] = 0.01;
  pt[18] =   4.625; y[18] =  1.44; e1[18] =  0.14; e2[18] =  0.09; m[18] = 0.01;
  pt[19] =   4.875; y[19] =  1.07; e1[19] =  0.12; e2[19] =  0.07; m[19] = 0.01;
  pt[20] =   5.125; y[20] =  9.19; e1[20] =  1.07; e2[20] =  0.60; m[20] = 0.001;
  pt[21] =   5.375; y[21] =  4.68; e1[21] =  0.72; e2[21] =  0.31; m[21] = 0.001;
  pt[22] =   5.625; y[22] =  3.28; e1[22] =  0.57; e2[22] =  0.22; m[22] = 0.001;
  pt[23] =   5.875; y[23] =  1.96; e1[23] =  0.44; e2[23] =  0.13; m[23] = 0.001;
  pt[24] =   6.250; y[24] =  1.17; e1[24] =  0.22; e2[24] =  0.08; m[24] = 0.001;
  pt[25] =   6.750; y[25] =  7.98; e1[25] =  1.76; e2[25] =  0.52; m[25] = 0.0001;
  pt[26] =   7.500; y[26] =  1.91; e1[26] =  0.73; e2[26] =  0.13; m[26] = 0.0001;

	double sum = 0;

	for(int i=0; i<27; i++)
	{
		float dpt = 0.25;
		if ( i==24 || i==25 ) dpt = 0.50;
		else if ( i==26 ) dpt = 1.0;
		sum += y[i]*m[i]*pt[i]*dpt*1.0*2*TMath::Pi();
		y[i] = y[i]*m[i]/42.0*1e-6;
		e1[i] = sqrt(e1[i]*e1[i] + e2[i]*e2[i]);
		e1[i] = e1[i]*m[i]/42.0*1e-6;
	} 

	cout << sum << endl;

	double pt_dAu_bwd[27], y_dAu_bwd[27], e1_bwd[27], e2_bwd[27], e3_bwd[27];
	double pt_dAu_fwd[27], y_dAu_fwd[27], e1_fwd[27], e2_fwd[27], e3_fwd[27];
	pt_dAu_bwd[0] = 0.125, y_dAu_bwd[0] = 4.770e-07, e1_bwd[0]= 5.920e-08, e2_bwd[0]= 4.190e-08, e3_bwd[0]= 5.370e-10;
	pt_dAu_bwd[1] = 0.375, y_dAu_bwd[1] = 4.520e-07, e1_bwd[1]= 5.340e-08, e2_bwd[1]= 3.960e-08, e3_bwd[1]= 5.080e-10;
	pt_dAu_bwd[2] = 0.625, y_dAu_bwd[2] = 4.350e-07, e1_bwd[2]= 4.530e-08, e2_bwd[2]= 3.820e-08, e3_bwd[2]= 4.890e-10;
	pt_dAu_bwd[3] = 0.875, y_dAu_bwd[3] = 4.090e-07, e1_bwd[3]= 2.540e-08, e2_bwd[3]= 3.590e-08, e3_bwd[3]= 4.600e-10;
	pt_dAu_bwd[4] = 1.125, y_dAu_bwd[4] = 3.570e-07, e1_bwd[4]= 1.610e-08, e2_bwd[4]= 3.130e-08, e3_bwd[4]= 4.020e-10;
	pt_dAu_bwd[5] = 1.375, y_dAu_bwd[5] = 2.740e-07, e1_bwd[5]= 1.270e-08, e2_bwd[5]= 2.400e-08, e3_bwd[5]= 3.080e-10;
	pt_dAu_bwd[6] = 1.625, y_dAu_bwd[6] = 2.140e-07, e1_bwd[6]= 9.180e-09, e2_bwd[6]= 1.880e-08, e3_bwd[6]= 2.410e-10;
	pt_dAu_bwd[7] = 1.875, y_dAu_bwd[7] = 1.610e-07, e1_bwd[7]= 7.020e-09, e2_bwd[7]= 1.410e-08, e3_bwd[7]= 1.810e-10;
	pt_dAu_bwd[8] = 2.215, y_dAu_bwd[8] = 1.220e-07, e1_bwd[8]= 5.420e-09, e2_bwd[8]= 1.070e-08, e3_bwd[8]= 1.370e-10;
	pt_dAu_bwd[9] = 2.375, y_dAu_bwd[9] = 9.740e-08, e1_bwd[9]= 4.180e-09, e2_bwd[9]= 8.540e-09, e3_bwd[9]= 1.100e-10;
	pt_dAu_bwd[10] = 2.625, y_dAu_bwd[10] = 6.620e-08, e1_bwd[10]= 2.640e-09, e2_bwd[10]= 5.810e-09, e3_bwd[10]= 7.450e-11;
	pt_dAu_bwd[11] = 	2.875, y_dAu_bwd[11] = 4.630e-08, e1_bwd[11]= 2.080e-09, e2_bwd[11]= 4.060e-09, e3_bwd[11]= 5.210e-11;
	pt_dAu_bwd[12] = 	3.125, y_dAu_bwd[12] = 3.140e-08, e1_bwd[12]= 1.600e-09, e2_bwd[12]= 2.750e-09, e3_bwd[12]= 3.530e-11;
	pt_dAu_bwd[13] = 	3.375, y_dAu_bwd[13] = 2.090e-08, e1_bwd[13]= 1.170e-09, e2_bwd[13]= 1.840e-09, e3_bwd[13]= 2.350e-11;
	pt_dAu_bwd[14] = 	3.625, y_dAu_bwd[14] = 1.720e-08, e1_bwd[14]= 1.130e-09, e2_bwd[14]= 1.510e-09, e3_bwd[14]= 1.930e-11;
	pt_dAu_bwd[15] = 	3.875, y_dAu_bwd[15] = 1.010e-08, e1_bwd[15]= 7.590e-10, e2_bwd[15]= 8.910e-10, e3_bwd[15]= 1.140e-11;
	pt_dAu_bwd[16] = 	4.125, y_dAu_bwd[16] = 7.650e-09, e1_bwd[16]= 7.190e-10, e2_bwd[16]= 6.720e-10, e3_bwd[16]= 8.610e-12;
	pt_dAu_bwd[17] = 	4.375, y_dAu_bwd[17] = 4.560e-09, e1_bwd[17]= 4.890e-10, e2_bwd[17]= 4.000e-10, e3_bwd[17]= 5.130e-12;
	pt_dAu_bwd[18] = 	4.625, y_dAu_bwd[18] = 4.050e-09, e1_bwd[18]= 5.240e-10, e2_bwd[18]= 3.550e-10, e3_bwd[18]= 4.560e-12;
	pt_dAu_bwd[19] = 	4.875, y_dAu_bwd[19] = 2.450e-09, e1_bwd[19]= 4.100e-10, e2_bwd[19]= 2.150e-10, e3_bwd[19]= 2.760e-12;
	pt_dAu_bwd[20] = 	5.125, y_dAu_bwd[20] = 1.730e-09, e1_bwd[20]= 3.430e-10, e2_bwd[20]= 1.520e-10, e3_bwd[20]= 1.950e-12;
	pt_dAu_bwd[21] = 	5.375, y_dAu_bwd[21] = 1.230e-09, e1_bwd[21]= 3.150e-10, e2_bwd[21]= 1.080e-10, e3_bwd[21]= 1.380e-12;
	pt_dAu_bwd[22] = 	5.625, y_dAu_bwd[22] = 5.820e-10, e1_bwd[22]= 1.750e-10, e2_bwd[22]= 5.110e-11, e3_bwd[22]= 6.550e-13;
	pt_dAu_bwd[23] = 	5.875, y_dAu_bwd[23] = 5.740e-10, e1_bwd[23]= 2.510e-10, e2_bwd[23]= 5.030e-11, e3_bwd[23]= 6.460e-13;
	pt_dAu_bwd[24] = 	6.250, y_dAu_bwd[24] = 2.940e-10, e1_bwd[24]= 1.090e-10, e2_bwd[24]= 2.580e-11, e3_bwd[24]= 3.310e-13;
	pt_dAu_bwd[25] = 	6.750, y_dAu_bwd[25] = 1.450e-10, e1_bwd[25]= 5.410e-11, e2_bwd[25]= 1.270e-11, e3_bwd[25]= 1.630e-13;
	pt_dAu_bwd[26] = 	7.500, y_dAu_bwd[26] = 5.660e-11, e1_bwd[26]= 2.830e-11, e2_bwd[26]= 4.970e-12, e3_bwd[26]= 6.370e-14;

	y_dAu_fwd[ 0] = 4.360e-07, e1_fwd[ 0] =  3.530e-08, e2_fwd[ 0] =  4.020e-08, e3_fwd[ 0] =  4.900e-10;
	y_dAu_fwd[ 1] = 4.030e-07, e1_fwd[ 1] =  3.210e-08, e2_fwd[ 1] =  3.720e-08, e3_fwd[ 1] =  4.530e-10;
	y_dAu_fwd[ 2] = 3.390e-07, e1_fwd[ 2] =  2.230e-08, e2_fwd[ 2] =  3.120e-08, e3_fwd[ 2] =  3.810e-10;
	y_dAu_fwd[ 3] = 3.000e-07, e1_fwd[ 3] =  1.560e-08, e2_fwd[ 3] =  2.760e-08, e3_fwd[ 3] =  3.370e-10;
	y_dAu_fwd[ 4] = 2.440e-07, e1_fwd[ 4] =  1.080e-08, e2_fwd[ 4] =  2.250e-08, e3_fwd[ 4] =  2.740e-10;
	y_dAu_fwd[ 5] = 1.960e-07, e1_fwd[ 5] =  8.110e-09, e2_fwd[ 5] =  1.810e-08, e3_fwd[ 5] =  2.200e-10;
	y_dAu_fwd[ 6] = 1.600e-07, e1_fwd[ 6] =  5.340e-09, e2_fwd[ 6] =  1.480e-08, e3_fwd[ 6] =  1.800e-10;
	y_dAu_fwd[ 7] = 1.130e-07, e1_fwd[ 7] =  3.950e-09, e2_fwd[ 7] =  1.040e-08, e3_fwd[ 7] =  1.270e-10;
	y_dAu_fwd[ 8] = 8.600e-08, e1_fwd[ 8] =  2.750e-09, e2_fwd[ 8] =  7.930e-09, e3_fwd[ 8] =  9.670e-11;
	y_dAu_fwd[ 9] = 6.370e-08, e1_fwd[ 9] =  2.190e-09, e2_fwd[ 9] =  5.870e-09, e3_fwd[ 9] =  7.170e-11;
	y_dAu_fwd[10] = 4.360e-08, e1_fwd[10] =  1.580e-09, e2_fwd[10] =  4.020e-09, e3_fwd[10] =  4.900e-11;
	y_dAu_fwd[11] = 3.270e-08, e1_fwd[11] =  1.240e-09, e2_fwd[11] =  3.010e-09, e3_fwd[11] =  3.680e-11;
	y_dAu_fwd[12] = 2.300e-08, e1_fwd[12] =  9.690e-10, e2_fwd[12] =  2.120e-09, e3_fwd[12] =  2.590e-11;
	y_dAu_fwd[13] = 1.420e-08, e1_fwd[13] =  6.860e-10, e2_fwd[13] =  1.310e-09, e3_fwd[13] =  1.600e-11;
	y_dAu_fwd[14] = 1.070e-08, e1_fwd[14] =  5.890e-10, e2_fwd[14] =  9.830e-10, e3_fwd[14] =  1.200e-11;
	y_dAu_fwd[15] = 7.650e-09, e1_fwd[15] =  5.320e-10, e2_fwd[15] =  7.050e-10, e3_fwd[15] =  8.610e-12;
	y_dAu_fwd[16] = 6.150e-09, e1_fwd[16] =  4.740e-10, e2_fwd[16] =  5.670e-10, e3_fwd[16] =  6.920e-12;
	y_dAu_fwd[17] = 4.150e-09, e1_fwd[17] =  3.990e-10, e2_fwd[17] =  3.820e-10, e3_fwd[17] =  4.670e-12;
	y_dAu_fwd[18] = 2.910e-09, e1_fwd[18] =  3.250e-10, e2_fwd[18] =  2.680e-10, e3_fwd[18] =  3.270e-12;
	y_dAu_fwd[19] = 2.030e-09, e1_fwd[19] =  2.770e-10, e2_fwd[19] =  1.870e-10, e3_fwd[19] =  2.280e-12;
	y_dAu_fwd[20] = 1.230e-09, e1_fwd[20] =  2.140e-10, e2_fwd[20] =  1.140e-10, e3_fwd[20] =  1.380e-12;
	y_dAu_fwd[21] = 8.540e-10, e1_fwd[21] =  1.920e-10, e2_fwd[21] =  7.870e-11, e3_fwd[21] =  9.610e-13;
	y_dAu_fwd[22] = 5.410e-10, e1_fwd[22] =  1.200e-10, e2_fwd[22] =  4.990e-11, e3_fwd[22] =  6.090e-13;
	y_dAu_fwd[23] = 3.890e-10, e1_fwd[23] =  1.000e-10, e2_fwd[23] =  3.580e-11, e3_fwd[23] =  4.380e-13;
	y_dAu_fwd[24] = 2.600e-10, e1_fwd[24] =  6.590e-11, e2_fwd[24] =  2.400e-11, e3_fwd[24] =  2.920e-13;
	y_dAu_fwd[25] = 1.410e-10, e1_fwd[25] =  4.580e-11, e2_fwd[25] =  1.300e-11, e3_fwd[25] =  1.590e-13;
	y_dAu_fwd[26] = 4.390e-11, e1_fwd[26] =  2.310e-11, e2_fwd[26] =  4.050e-12, e3_fwd[26] =  4.940e-14;

	for (int i=0; i<27; i++){
		e1_bwd[i] = sqrt(e1_bwd[i]*e1_bwd[i] + e2_bwd[i]*e2_bwd[i] + e3_bwd[i]*e3_bwd[i]);

		pt_dAu_fwd[i] = pt_dAu_bwd[i];
		e1_fwd[i] = sqrt(e1_fwd[i]*e1_fwd[i] + e2_fwd[i]*e2_fwd[i] + e3_fwd[i]*e3_fwd[i]);
	}


	TGraphErrors *pt_Jpsi = new TGraphErrors(27,pt,y,0,e1);
	TGraphErrors *pt_Jpsi_dAu_bwd = new TGraphErrors(27, pt_dAu_bwd, y_dAu_bwd, 0, e1_bwd);
	TGraphErrors *pt_Jpsi_dAu_fwd = new TGraphErrors(27, pt_dAu_fwd, y_dAu_fwd, 0, e1_fwd);


  ///////////Fit function
	pt_Jpsi_fit = new TF1("pt_Jpsi_fit"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",0,10);
  //pt_Jpsi_fit->SetParameter(0,0.36);
  //pt_Jpsi_fit->SetParameter(1,-0.009);
  //pt_Jpsi_fit->SetParameter(2,1000);
  //pt_Jpsi_fit->SetParameter(3,3.24);
  //pt_Jpsi_fit->SetParameter(4,12.69766);
  pt_Jpsi_fit->SetParameter(0,0.1);
  pt_Jpsi_fit->SetParameter(1,-0.01);
  pt_Jpsi_fit->SetParameter(2,1e-7);
  pt_Jpsi_fit->SetParameter(3,1);
  pt_Jpsi_fit->SetParameter(4,10);

	pt_Jpsi_fit_bwd = new TF1("pt_Jpsi_fit_bwd"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",0,10);
  pt_Jpsi_fit_bwd->SetParameter(0,0.1);
  pt_Jpsi_fit_bwd->SetParameter(1,-0.01);
  pt_Jpsi_fit_bwd->SetParameter(2,5e-7);
  pt_Jpsi_fit_bwd->SetParameter(3,1);
  pt_Jpsi_fit_bwd->SetParameter(4,10);

	pt_Jpsi_fit_fwd = new TF1("pt_Jpsi_fit_fwd"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",0,10);
  pt_Jpsi_fit_fwd->SetParameter(0,0.1);
  pt_Jpsi_fit_fwd->SetParameter(1,-0.01);
  pt_Jpsi_fit_fwd->SetParameter(2,5e-7);
  pt_Jpsi_fit_fwd->SetParameter(3,1);
  pt_Jpsi_fit_fwd->SetParameter(4,10);
  
	TCanvas *c1 = new TCanvas("c1","c1",1.1*800,400);
	c1->Divide(2,1);

	c1->cd(1);
	SetPadStyle();
	gPad->SetLogy();
	htmp = (TH1F*)gPad->DrawFrame(0,5e-13,8,1e-5);
	SetHistoStyle("p_{T} (GeV/c)","B_{#mu#mu}/2#pip_{T} d^{2}N/dydp_{T} [(GeV/c)^{-2}]");
	pt_Jpsi->SetMarkerStyle(20);
	pt_Jpsi->SetLineColor(1);
	pt_Jpsi->SetLineWidth(2);
	pt_Jpsi->Draw("P");

	pt_Jpsi_dAu_bwd->SetMarkerStyle(20);
	pt_Jpsi_dAu_bwd->SetMarkerColor(4);
	pt_Jpsi_dAu_bwd->SetLineColor(4);
	pt_Jpsi_dAu_bwd->SetLineWidth(2);
	pt_Jpsi_dAu_bwd->Draw("p");

	pt_Jpsi_dAu_fwd->SetMarkerStyle(20);
	pt_Jpsi_dAu_fwd->SetMarkerColor(2);
	pt_Jpsi_dAu_fwd->SetLineColor(2);
	pt_Jpsi_dAu_fwd->SetLineWidth(2);
	pt_Jpsi_dAu_fwd->Draw("p");

  pt_Jpsi_dAu_bwd->Fit(pt_Jpsi_fit_bwd,"R0Q",0,10);
	pt_Jpsi_fit_bwd->SetLineColor(4);
	pt_Jpsi_fit_bwd->SetLineWidth(3);
  pt_Jpsi_fit_bwd->Draw("same");

  pt_Jpsi_dAu_fwd->Fit(pt_Jpsi_fit_fwd,"R0Q",0,10);
	pt_Jpsi_fit_fwd->SetLineColor(2);
	pt_Jpsi_fit_fwd->SetLineWidth(3);
  pt_Jpsi_fit_fwd->Draw("same");

  pt_Jpsi->Fit(pt_Jpsi_fit,"R0Q",0,10);
	pt_Jpsi_fit->SetLineColor(1);
	pt_Jpsi_fit->SetLineWidth(3);
  pt_Jpsi_fit->Draw("same");

	TLegend *leg = new TLegend(0.3,0.75,0.9,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(pt_Jpsi,"p+p, 1.2<|y|<2.2, PPG104","P");
	leg->AddEntry(pt_Jpsi_dAu_fwd,"d+Au, 1.2<y<2.2, PPG125","P");
	leg->AddEntry(pt_Jpsi_dAu_bwd,"d+Au, -2.2<y<-1.2, PPG125","P");
	leg->Draw();

	//return;

	TF1 *pt_Jpsi_gen = new TF1("pt_Jpsi_gen"," [2]*x/pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",0,20);

  // double A = pt_Jpsi_fit->GetParameter(0);
  // double B = pt_Jpsi_fit->GetParameter(1);
  // double n = pt_Jpsi_fit->GetParameter(2);
  cout<<endl; cout<<endl;
  cout<<"pt_wt = new TF1(\"pt_Jpsi_fit\",\" [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])\",0,20);"<<endl;
	for(int i = 0;i<5;i++)
	{
		cout<<"pt_wt->SetParameter("<<i<<", "<<pt_Jpsi_fit->GetParameter(i)<<" );"<<endl;  
		pt_Jpsi_gen->SetParameter(i,pt_Jpsi_fit->GetParameter(i));
	}
  cout<<endl; cout<<endl;

  double chisq = pt_Jpsi_fit->GetChisquare();
  double ndf =  pt_Jpsi_fit->GetNDF();
  double chisqdf = chisq/ndf;
  cout<<" Chisquare/ndf for pt fit function: " << chisq <<"/"<<ndf<<" = "<< chisqdf << endl;
  
  //break;

	double rap[12], rap_y[12], rap_e1[12], rap_e2[12];
	double rap_dAu[12], rap_y_dAu[12], rap_e1_dAu[12], rap_e2_dAu[12], rap_e3_dAu[12], rap_e4_dAu[12];

	rap[0]	= -2.075;	rap_y[0]	= 17.6; rap_e1[0]		= 0.5; rap_e2[0]	= 1.5;
	rap[1]	= -1.825;	rap_y[1]	= 24.4; rap_e1[1]		= 0.4; rap_e2[1]	= 1.9;
	rap[2]	= -1.575;	rap_y[2]	= 31.5; rap_e1[2]		= 0.5; rap_e2[2]	= 2.2;
	rap[3]	= -1.325;	rap_y[3]	= 41.2; rap_e1[3]		= 1.1; rap_e2[3]	= 5.3;
	rap[4]	= -0.3;		rap_y[4]	= 49.0; rap_e1[4]		= 2.1; rap_e2[4]	= 5.4;
	rap[5]	= 0.0;		rap_y[5]	= 45.6; rap_e1[5]		= 1.6; rap_e2[5]	= 5.0;
	rap[6]	= 0.35;		rap_y[6]	= 46.1; rap_e1[6]		= 1.9; rap_e2[6]	= 5.1;
	rap[7]	= 1.325;	rap_y[7]	= 40.7; rap_e1[7]		= 1.2; rap_e2[7]	= 5.7;
	rap[8]	= 1.575;	rap_y[8]	= 33.6; rap_e1[8]		= 0.7; rap_e2[8]	= 3.0;
	rap[9]	= 1.825;	rap_y[9]	= 25.6; rap_e1[9]		= 0.4; rap_e2[9]	= 2.5;
	rap[10]	= 2.075;	rap_y[10]	= 18.9; rap_e1[10]	= 0.4; rap_e2[10]	= 1.9;
	rap[11]	= 2.325;	rap_y[11]	= 13.9; rap_e1[11]	= 0.9; rap_e2[11]	= 1.4;

	rap_dAu[0] = -2.075, rap_y_dAu[0] =  2.50E-06, rap_e1_dAu[0] =  5.70E-08, rap_e2_dAu[0] =  5.10E-08, rap_e3_dAu[0] =  2.20E-07, rap_e4_dAu[0] =  5.30E-08;
	rap_dAu[1] = -1.825, rap_y_dAu[1] =  4.00E-06, rap_e1_dAu[1] =  4.90E-08, rap_e2_dAu[1] =  4.20E-08, rap_e3_dAu[1] =  3.50E-07, rap_e4_dAu[1] =  8.50E-08;
	rap_dAu[2] = -1.575, rap_y_dAu[2] =  5.40E-06, rap_e1_dAu[2] =  6.00E-08, rap_e2_dAu[2] =  5.00E-08, rap_e3_dAu[2] =  4.80E-07, rap_e4_dAu[2] =  1.10E-07;
	rap_dAu[3] = -1.325, rap_y_dAu[3] =  7.20E-06, rap_e1_dAu[3] =  1.40E-07, rap_e2_dAu[3] =  1.10E-07, rap_e3_dAu[3] =  6.30E-07, rap_e4_dAu[3] =  1.50E-07;
	rap_dAu[4] = +1.325, rap_y_dAu[4] =  5.50E-06, rap_e1_dAu[4] =  1.10E-07, rap_e2_dAu[4] =  8.50E-08, rap_e3_dAu[4] =  5.10E-07, rap_e4_dAu[4] =  1.20E-07;
	rap_dAu[5] = +1.575, rap_y_dAu[5] =  4.20E-06, rap_e1_dAu[5] =  4.60E-08, rap_e2_dAu[5] =  3.80E-08, rap_e3_dAu[5] =  3.90E-07, rap_e4_dAu[5] =  8.90E-08;
	rap_dAu[6] = +1.825, rap_y_dAu[6] =  3.20E-06, rap_e1_dAu[6] =  3.60E-08, rap_e2_dAu[6] =  3.10E-08, rap_e3_dAu[6] =  3.00E-07, rap_e4_dAu[6] =  6.80E-08;
	rap_dAu[7] = +2.075, rap_y_dAu[7] =  2.50E-06, rap_e1_dAu[7] =  3.60E-08, rap_e2_dAu[7] =  3.10E-08, rap_e3_dAu[7] =  2.30E-07, rap_e4_dAu[7] =  5.20E-08;
	rap_dAu[8] = +2.325, rap_y_dAu[8] =  1.50E-06, rap_e1_dAu[8] =  5.60E-08, rap_e2_dAu[8] =  4.90E-08, rap_e3_dAu[8] =  1.40E-07, rap_e4_dAu[8] =  3.20E-08;
	rap_dAu[ 9] = -0.3, rap_y_dAu[ 9] = 8.93e-07, rap_e1_dAu[ 9] =  2.2e-08, rap_e2_dAu[ 9] =  6.0e-08, rap_e3_dAu[ 9] =  3.1e-08;
	rap_dAu[10] = +0.0, rap_y_dAu[10] = 8.30e-07, rap_e1_dAu[10] =  1.8e-08, rap_e2_dAu[10] =  5.6e-08, rap_e3_dAu[10] =  2.8e-08;
	rap_dAu[11] = +0.3, rap_y_dAu[11] = 8.61e-07, rap_e1_dAu[11] =  2.2e-08, rap_e2_dAu[11] =  5.8e-08, rap_e3_dAu[11] =  2.9e-08;

	sum = 0;

	for (int i=0; i<12; i++){

		if ( fabs(rap[i])>1.2 && fabs(rap[i])<2.2 ){
			sum += rap_y[i]*0.25;
		}

		rap_y[i] /= 42.0;
		rap_e1[i] = sqrt(rap_e1[i]*rap_e1[i] + rap_e2[i]*rap_e2[i]);
		rap_e1[i] /= 42.0;

	}

	cout << sum/2.0 << endl;

	for (int i=0; i<9; i++){
		rap_y_dAu[i] = rap_y_dAu[i]*1e6/7.6;
		rap_e1_dAu[i] = sqrt(rap_e1_dAu[i]*rap_e1_dAu[i] + rap_e2_dAu[i]*rap_e2_dAu[i] + rap_e3_dAu[i]*rap_e3_dAu[i] + rap_e4_dAu[i]*rap_e4_dAu[i]);
		rap_e1_dAu[i] = rap_e1_dAu[i]*1e6/7.6;
	}

	for (int i=9; i<12; i++){
		rap_y_dAu[i] *= 1e6;
		rap_e1_dAu[i] = sqrt(rap_e1_dAu[i]*rap_e1_dAu[i] + rap_e2_dAu[i]*rap_e2_dAu[i] + rap_e3_dAu[i]*rap_e3_dAu[i]);
		rap_e1_dAu[i] *= 1e6;
	}

  //TF1 *rap_Jpsi_fit = new TF1("rap_Jpsi_fit","[0]*exp(-0.5*((x-[1])/[2])**2)",1.0,2.4);

	TGraphErrors *rap_Jpsi = new TGraphErrors(12, rap, rap_y, 0, rap_e1);
  rap_Jpsi->SetLineWidth(2);
	rap_Jpsi->SetMarkerStyle(20);

	TGraphErrors *rap_Jpsi_dAu = new TGraphErrors(12, rap_dAu, rap_y_dAu, 0, rap_e1_dAu);
  rap_Jpsi_dAu->SetLineWidth(2);
	rap_Jpsi_dAu->SetMarkerStyle(20);
	rap_Jpsi_dAu->SetMarkerColor(2);
	rap_Jpsi_dAu->SetLineColor(2);
  
  TF1 *rap_Jpsi_fit = new TF1("rap_Jpsi_fit","[0]*exp(-0.5*((x-[1])/[2])**2) + [0]*exp(-0.5*((x+[1])/[2])**2)",-2.4,2.4);
  rap_Jpsi_fit->SetParameter(0,1);
  rap_Jpsi_fit->SetParameter(1,-1);
  rap_Jpsi_fit->SetParameter(2,0.5);
	rap_Jpsi_fit->SetLineColor(1);
	rap_Jpsi_fit->SetLineWidth(3);

  TF1 *rap_Jpsi_fit_dAu = new TF1("rap_Jpsi_fit_dAu","[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)",-2.4,2.4);
  rap_Jpsi_fit_dAu->SetParameter(0,1);
  rap_Jpsi_fit_dAu->SetParameter(1,-1);
  rap_Jpsi_fit_dAu->SetParameter(2,0.5);
  rap_Jpsi_fit_dAu->SetParameter(3,1);
  rap_Jpsi_fit_dAu->SetParameter(4,1);
  rap_Jpsi_fit_dAu->SetParameter(5,0.5);
	rap_Jpsi_fit_dAu->SetLineColor(2);
	rap_Jpsi_fit_dAu->SetLineWidth(3);
  
	c1->cd(2);
	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(-3, 0, 3, 1.5);
	SetHistoStyle("y","B_{ll} dN/dy (A.U.)");
  rap_Jpsi->Draw("P");
  rap_Jpsi->Fit(rap_Jpsi_fit,"R0Q");
	rap_Jpsi_fit->SetRange(-3,3);
  rap_Jpsi_fit->Draw("same");

	rap_Jpsi_dAu->Draw("p");
  rap_Jpsi_dAu->Fit(rap_Jpsi_fit_dAu,"R0Q");
	rap_Jpsi_fit_dAu->SetRange(-3,3);
  rap_Jpsi_fit_dAu->Draw("same");

	TLegend *leg = new TLegend(0.2,0.75,0.6,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(rap_Jpsi,"p+p, PPG104","PL");
	leg->AddEntry(rap_Jpsi_dAu,"d+Au, PPG109","PL");
	leg->Draw();

  cout<<endl; cout<<endl;

  cout<<"rap_wt = new TF1(\"rapt_Jpsi_fit\",\"[0]*exp(-0.5*((x-[1])/[2])**2) + [0]*exp(-0.5*((x+[1])/[2])**2)\",-2.5,2.5);"<<endl;
	for(int i = 0;i<3;i++)
	{
		cout<<"rap_wt->SetParameter("<<i<<", "<<rap_Jpsi_fit->GetParameter(i)<<" );"<<endl;  
	}

	/*
	TFile *outfile = new TFile("outfile_Jpsi_fit_func.root","recreate");
	rap_Jpsi_fit->Write("fit_func_rap_pp");
	rap_Jpsi_fit_dAu->Write("fit_func_rap_dAu");
	pt_Jpsi_fit->Write("fit_func_pT_pp");
	pt_Jpsi_fit_bwd->Write("fit_func_pT_dAu_bwd");
	pt_Jpsi_fit_fwd->Write("fit_func_pT_dAu_fwd");
	outfile->Close();
	*/

	/*
	TH1F *hpT = new TH1F("hpT","hpT",100,0,10);
	TH1F *hy = new TH1F("hy","hy",100,1.2,2.2);
	for (int ii=0; ii<1e6; ii++){
		hpT->Fill(pt_Jpsi_gen->GetRandom(0,10));
		hy->Fill(rap_Jpsi_fit->GetRandom(1.2,2.2));
	}

	TCanvas *2 = new TCanvas("c2","c2",800,400);
	c2->Divide(2,1);
	c2->cd(1);
	SetPadStyle();
	htmp = (TH1F*)hpT;
	SetHistoStyle();
	hpT->Draw("");

	c2->cd(2);
	SetPadStyle();
	htmp = (TH1F*)hy;
	SetHistoStyle();
	hy->Draw("");
	*/

}

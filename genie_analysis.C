#define GENIE_ANALYSIS_C
#include "genie_analysis.h"
#include "Constants.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TF1.h>

#include <exception>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;

// -----------------------------------------------------------------------------------------------------------------------------------------------

double DeltaAlphaTFunction(TVector3 MuonVector,TVector3 ProtonVector) {
			
	// STV Calculation		
			
	TVector3 MuonVectorTrans;
	MuonVectorTrans.SetXYZ(MuonVector.X(),MuonVector.Y(),0.);
	double MuonVectorTransMag = MuonVectorTrans.Mag();
			
	TVector3 ProtonVectorTrans;
	ProtonVectorTrans.SetXYZ(ProtonVector.X(),ProtonVector.Y(),0.);
	double ProtonVectorTransMag = ProtonVectorTrans.Mag();

	TVector3 PtVector = MuonVectorTrans + ProtonVectorTrans;

	double fPt = PtVector.Mag();

	double fDeltaAlphaT = TMath::ACos( (- MuonVectorTrans * PtVector) / ( MuonVectorTransMag * fPt ) ) * 180./TMath::Pi();
	if (fDeltaAlphaT > 180.) { fDeltaAlphaT -= 180.; }
	if (fDeltaAlphaT < 0.) { fDeltaAlphaT += 180.; }
	
	return fDeltaAlphaT;

}

// -----------------------------------------------------------------------------------------------------------------------------------------------

double DeltaPhiTFunction(TVector3 MuonVector,TVector3 ProtonVector) {
			
	// STV Calculation		
			
	TVector3 MuonVectorTrans;
	MuonVectorTrans.SetXYZ(MuonVector.X(),MuonVector.Y(),0.);
	double MuonVectorTransMag = MuonVectorTrans.Mag();
			
	TVector3 ProtonVectorTrans;
	ProtonVectorTrans.SetXYZ(ProtonVector.X(),ProtonVector.Y(),0.);
	double ProtonVectorTransMag = ProtonVectorTrans.Mag();

	double fDeltaPhiT = TMath::ACos( (- MuonVectorTrans * ProtonVectorTrans) / ( MuonVectorTransMag * ProtonVectorTransMag ) ) * 180./TMath::Pi();
	if (fDeltaPhiT > 180.) { fDeltaPhiT -= 180.; }
	if (fDeltaPhiT < 0.) { fDeltaPhiT += 180.; }
	
	return fDeltaPhiT;

}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Loading all the constants from Constant.h (e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565,
// H3_bind_en, He4_bind_en, C12_bind_en, B_bind_en, He3_bind_en, D2_bind_en, Fe_bind_en, Mn_bind_en

void genie_analysis::Loop() {

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	int NSectors = 6;
	int NDeltaAlphaTBins = 36; double DeltaAlphaTMin = 0.; double DeltaAlphaTMax = 180.;
	int NDeltaPhiTBins = 36; double DeltaPhiTMin = 0.; double DeltaPhiTMax = 180.;

	const int NInt = 6; // All Interactions = 0, QE = 1, MEC = 2, RES = 3, DIS = 4, COH = 5

	// ---------------------------------------------------------------------------------------------------------------

	// Do we want to apply fiducials & the acceptance map weights
	// Do we want a truth level study ? if so, stop ditching sectors

	bool TruthLevel0piSignalStudy = true;
	bool TruthLevel1p0piSignalStudy = false;

	bool ApplyThetaSlice = false; double MinThetaSlice = 36, MaxThetaSlice = 39;

	// ---------------------------------------------------------------------------------------------------------------

	std::map<std::string,double> bind_en;
	std::map<std::string,double> target_mass;
	std::map<std::string,double> residual_target_mass;
	std::map<std::string, double> Ecal_offset; // that might not be necessary for simulation data

	target_name = ftarget; //std string for target name
	en_beam["1161"]=1.161;
	en_beam["2261"]=2.261;
	en_beam["4461"]=4.461;

	en_beam_Ecal["1161"]=1.161;
	en_beam_Ecal["2261"]=2.261;
	en_beam_Ecal["4461"]=4.461;

	en_beam_Eqe["1161"]=1.161;
	en_beam_Eqe["2261"]=2.261;
	en_beam_Eqe["4461"]=4.461;

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	double Wcut = 2; //cut for all beam energies < 2
	double Q2cut = 0; // cut for 1.1 GeV > 0.1, for 2.2 GeV > 0.4 and 4.4 GeV > 0.8

	const int n_slice=3; // Stick to the 3 slices
	const double pperp_min[n_slice]={0.,0.2,0.4};
	const double pperp_max[n_slice]={0.2,0.4,10.};

	TVector3 V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) //1.1 GeV Configuration parameters and cuts
	{
		Q2cut = 0.1;
	}


	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) //2.2 GeV Configuration parameters and cuts
	{
		Q2cut = 0.4;
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) //4.4 GeV Configuration parameters and cuts
	{
		Q2cut = 0.8;
	}

	// Further constants for binding energies and target masses

	Ecal_offset["3He"]  = 0.004;
	Ecal_offset["4He"]  = 0.005;
	Ecal_offset["C12"]  = 0.005;
	Ecal_offset["56Fe"] = 0.011;

	bind_en["3He"]  = He3_bind_en-D2_bind_en + Ecal_offset["3He"]; // the offset is used to shift the peak to be at 0
	bind_en["4He"]  = He4_bind_en-H3_bind_en + Ecal_offset["4He"];
	bind_en["C12"]  = C12_bind_en-B_bind_en	+ Ecal_offset["C12"];
	bind_en["56Fe"] = Fe_bind_en-Mn_bind_en	+ Ecal_offset["56Fe"];
	bind_en["CH2"]  = C12_bind_en-B_bind_en;

	target_mass["3He"] = 2*m_prot+m_neut-He3_bind_en;
	target_mass["4He"] = 2*m_prot+2*m_neut-He4_bind_en;
	target_mass["C12"] = 6*m_prot+6*m_neut-C12_bind_en;
	target_mass["56Fe"]= 26*m_prot+30*m_neut-Fe_bind_en;
	target_mass["CH2"] = 6*m_prot+6*m_neut-C12_bind_en;

	residual_target_mass["3He"] = m_prot+m_neut-D2_bind_en;
	residual_target_mass["4He"] = m_prot+2*m_neut-H3_bind_en;
	residual_target_mass["C12"] = 5*m_prot+6*m_neut-B_bind_en;
	residual_target_mass["56Fe"]= 25*m_prot+30*m_neut-Mn_bind_en;
	residual_target_mass["CH2"] = 25*m_prot+30*m_neut-Mn_bind_en;

	// ----------------------------------------------------------------------------------------

	TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);
	TLorentzVector V4_target(0,0,0,target_mass[ftarget]);

	// ---------------------------------------------------------------------------------------------------------------

	//Output file definition

	TString FileName = Form("genie_e2a_ep_%s_%s.root",ftarget.c_str(),fbeam_en.c_str());
	TFile *file_out = new TFile(FileName, "Recreate");

	// ---------------------------------------------------------------------------------------------------------------

	fiducialcut->InitPiMinusFit(fbeam_en);

	// ---------------------------------------------------------------------------------------------------------------

	// Binning for energy reconstruction histograms && feeddown

	int n_bins;
	double *x_values;
	double *x_qe;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=17;i++) { x_values[i]=0.4+i*0.04; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en]; }
		for (int i=0;i<=20;i++) { x_values[i+18]=1.08+(i+1)*0.02; x_qe[i+18] = (x_values[i+18] - en_beam[fbeam_en]) / en_beam[fbeam_en]; }
	}

	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
		n_bins=54;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=23;i++) { x_values[i]=i*0.09; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=30;i++) { x_values[i+24]=2.07+(i+1)*0.03; x_qe[i+24] = (x_values[i+24] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=21;i++)	{ x_values[i]=i*0.2; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=16;i++)	{ x_values[i+22]=4.2+(i+1)*0.05; x_qe[i+22] = (x_values[i+22] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}	

	// ---------------------------------------------------------------------------------------------------------------

	TH1F *h1_Etot_Npi0[n_slice],*h1_Erec_Npi0[n_slice];

	for(int h = 0; h < n_slice; h++){
		
		h1_Etot_Npi0[h] = new TH1F(Form("h1_Etot_Npi0_%d",h+1),"",n_bins,x_values);
		h1_Erec_Npi0[h] = new TH1F(Form("h1_Erec_Npi0_%d",h+1),"",n_bins,x_values);

	}	

	// ---------------------------------------------------------------------------------------------------------------		

	TH1F *h1_Nprot=new TH1F("h1_Nprot","",5,-0.5,4.5);
	TH1F *h1_Nprot_NonZeroProt=new TH1F("h1_Nprot_NonZeroProt","",4,0.5,4.5);
	TH1F *h1_Nphot=new TH1F("h1_Nphot","",10,-0.5,4.5);
	TH1F *h1_Npiphot=new TH1F("h1_Npiphot","",10,-0.5,4.5);
	TH1F *h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,-0.5,4.5);
	TH1F *h1_Npi=new TH1F("h1_Npi","",5,-0.5,4.5);
	TH1F *h1_Npi_NonZeroProt=new TH1F("h1_Npi_NonZeroProt","",4,-0.5,4.5);
	TH1F *h1_Npipl=new TH1F("h1_Npipl","",10,-0.5,4.5);
	TH1F *h1_Npimi=new TH1F("h1_Npimi","",10,-0.5,4.5);
	TH1F *h1_MissMomentum = new TH1F("MissMomentum","",80,0.,1.);

	TH1F *h1_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector[NSectors];
	TH1F *h1_InteractionBreakDown_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector[NInt][NSectors];

	for (int WhichSector = 0; WhichSector < NSectors; WhichSector++) {

		h1_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector[WhichSector]  = new TH1F("h1_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_"+TString(std::to_string(WhichSector)),"",6000,0.,6.);

		for (int WhichInt = 1; WhichInt < NInt; WhichInt++) {

			h1_InteractionBreakDown_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector[WhichInt][WhichSector]  = new TH1F("h1_"+TString(std::to_string(WhichInt))+"_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_"+TString(std::to_string(WhichSector)),"",6000,0.,6.);

		}

	}

	// -------------------------------------------------------------------------------------------------------

	// 0pi plots

	TH1F *h1_E_rec_0pi_frac_feed=new TH1F("h1_E_rec_0pi_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_rec_0pi = new TH1F("h_Erec_subtruct_piplpimi_noprot_3pi","",n_bins,x_values);

	// 1p0pi plots	

	TH1F *h1_E_tot_cut2_fracfeed = new TH1F("h1_E_tot_cut2_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_cut2 = new TH1F("epRecoEnergy_slice_0","",n_bins,x_values);

	// ---------------------------------------------------------------------------------------------------

	// Plots for interaction break down for GENIE samples

	TH1D* ECal_BreakDown[NInt];
	TH1D* EQE_BreakDown[NInt];

	TH1D* ECal_LowPmiss_BreakDown[NInt];
	TH1D* EQE_LowPmiss_BreakDown[NInt];

	TH1D* ECal_MidPmiss_BreakDown[NInt];
	TH1D* EQE_MidPmiss_BreakDown[NInt];

	TH1D* ECal_HighPmiss_BreakDown[NInt];
	TH1D* EQE_HighPmiss_BreakDown[NInt];

	TH1D* InclusiveEQE_BreakDown[NInt];
	TH1D* Pmiss_BreakDown[NInt];
		 
	TH1D* DeltaPhiT_BreakDown[NInt];
	TH1D* DeltaAlphaT_BreakDown[NInt];

	for (int WhichInt = 0; WhichInt < NInt; WhichInt++) {

		ECal_BreakDown[WhichInt] = new TH1D(Form("ECal_Int_%d",WhichInt),";E^{Cal} (GeV)",n_bins,x_values);
		EQE_BreakDown[WhichInt] = new TH1D(Form("EQE_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);

		ECal_LowPmiss_BreakDown[WhichInt] = new TH1D(Form("ECal_LowPmiss_Int_%d",WhichInt),";E^{Cal} (GeV)",n_bins,x_values);
		EQE_LowPmiss_BreakDown[WhichInt] = new TH1D(Form("EQE_LowPmiss_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);

		ECal_MidPmiss_BreakDown[WhichInt] = new TH1D(Form("ECal_MidPmiss_Int_%d",WhichInt),";E^{Cal} (GeV)",n_bins,x_values);
		EQE_MidPmiss_BreakDown[WhichInt] = new TH1D(Form("EQE_MidPmiss_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);

		ECal_HighPmiss_BreakDown[WhichInt] = new TH1D(Form("ECal_HighPmiss_Int_%d",WhichInt),";E^{Cal} (GeV)",n_bins,x_values);
		EQE_HighPmiss_BreakDown[WhichInt] = new TH1D(Form("EQE_HighPmiss_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);

		InclusiveEQE_BreakDown[WhichInt] = new TH1D(Form("InclusiveEQE_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);
		Pmiss_BreakDown[WhichInt] = new TH1D(Form("Pmiss_Int_%d",WhichInt),";P_{miss}^{#perp} [GeV/c]",80,0.,1.);
		
		DeltaPhiT_BreakDown[WhichInt] = new TH1D(Form("DeltaPhiT_Int_%d",WhichInt),\
								";#delta#phi_{T} [deg]",NDeltaPhiTBins,DeltaPhiTMin,DeltaPhiTMax);
								
		DeltaAlphaT_BreakDown[WhichInt] = new TH1D(Form("DeltaAlphaT_Int_%d",WhichInt),\
								";#delta#alpha_{T} [deg]",NDeltaAlphaTBins,DeltaAlphaTMin,DeltaAlphaTMax);	
				
	}	  						  

	// Signal Event Counter -> 1e1p0pi events
	int SignalEvents = 0;
	int QESignalEvents = 0;
	int MECSignalEvents = 0;
	int RESSignalEvents = 0;
	int DISSignalEvents = 0;
	int OtherSignalEvents = 0;

	// ---------------------------------------------------------------------------------------------------------------

	TF1 *myElectronFit = new TF1("myElectronFit","[0]+[1]/x",0.,5.);

	if (en_beam[fbeam_en] == 1.161) { myElectronFit->SetParameters(17,7); }
	if (en_beam[fbeam_en] == 2.261) { myElectronFit->SetParameters(16,10.5); }
	if (en_beam[fbeam_en] == 4.461) { myElectronFit->SetParameters(13.5,15); }

	// ---------------------------------------------------------------------------------------------------------------

	int TotalCounter = 0;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		//Read Entry
		int nb = GetEntry(jentry);
		if (nb == 0) { std::cout <<"Event loop: 0 byte read for entry " << jentry << ". Indicate failure in reading the file" <<	std::endl;}

		if (jentry%1000 == 0) {std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/double(fChain->GetEntries())*100. << " %"<< std::endl;}

		if( jentry%200000 == 0 )
		{
			gDirectory->Write("hist_Files", TObject::kOverwrite);
			//cout<<jentry<<endl;
		}

		TotalCounter ++;

		// ---------------------------------------------------------------------------------------------------------------

		std::string StoreEnergy = fbeam_en;

		// -------------------------------------------------------------------------------------------------------------------------

		// For GENIE samples, identify the interaction type

		int Interaction = -1;

		if (qel) { Interaction = 1; }
		if (mec) { Interaction = 2; }
		if (res) { Interaction = 3; }
		if (dis) { Interaction = 4; }

		// ---------------------------------------------------------------------------------------------------------------

		if(jentry == 0){ //first entry to initialize TorusCurrent, Fiducials and Subtraction classes

			//The TorusField has to be set before the Fiducialcut parameters are initialized
			if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. ) //1.1 GeV, we are not using the 1.1 GeV data with 1500 current field
			{
				 fTorusCurrent = 750;
			}
			else if( (en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) || (en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) ) //2.2 GeV	or 4.4 GeV
			{
				 fTorusCurrent = 2250;
			}
			else { std::cout << "genie_analysis::Loop(): fTorusCurrent could not be assigned" << std::endl;}

			fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
			fiducialcut->SetFiducialCutParameters(fbeam_en);

		}

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Counters for truth level studies

		int TrueElectronsAboveThreshold = 0;
		int TrueProtonsAboveThreshold = 0;
		int TrueChargedPionsAboveThreshold = 0;
		int TruePiPlusAboveThreshold = 0;
		int TruePiMinusAboveThreshold = 0;
		int TrueGammasAboveThreshold = 0;

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		double SmearedPe;
		double SmearedEe;

		// Outgoing e',	Uncorr and corrected are the same read from root file.
		//V4_el and V3_el will be changed by smearing for GENIE simulation data
		TLorentzVector V4_el(pxl,pyl,pzl,El);
		TLorentzVector V4_el_uncorr(pxl,pyl,pzl,El);
		TVector3 V3_el(pxl,pyl,pzl);

		double el_momentum = V3_el.Mag();
		double el_theta = V3_el.Theta();

		// ----------------------------------------------------------------------------------------------------------------------	

		//Smearing of Electron Vector from Simulation
		SmearedPe = pl;
		SmearedEe = sqrt( SmearedPe*SmearedPe + e_mass * e_mass );
		V3_el.SetXYZ(pxl,pyl,pzl);
		V4_el.SetPxPyPzE(V3_el.X(),V3_el.Y(),V3_el.Z(),SmearedEe);
		double phi_ElectronOut = V3_el.Phi(); //in Radians

		V3_el.SetPhi(phi_ElectronOut + TMath::Pi() ); // Vec.Phi() is between (-180,180), GENIE coordinate system flipped with respect to CLAS

		phi_ElectronOut += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
		el_momentum = V3_el.Mag(); //Momentum after smearing
		el_theta = V3_el.Theta(); //Angle after smearing	

		// ----------------------------------------------------------------------------------------------------------------------

		double theta_min = myElectronFit->Eval(el_momentum); // in deg
		if (el_theta*180./TMath::Pi() < theta_min) { continue; }

		// ----------------------------------------------------------------------------------------------------------------------

		// Explicit cuts on electron momentum

		if (fbeam_en=="1161" && el_momentum < 0.4) { continue; }
		if (fbeam_en=="2261" && el_momentum < 0.55) { continue; }
		if (fbeam_en=="4461" && el_momentum < 1.1) { continue; }

		//Definition as for data. It is also correct for GENIE simulation data since V3_el is rotated above by 180 degree in phi
		double el_phi_mod = V3_el.Phi()*TMath::RadToDeg()  + 30; //Add 30 degree for plotting and photon phi cut
		if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree

		double reco_q3 = (V4_el-V4_beam).Rho();
		double reco_Q2 = -(V4_el-V4_beam).Mag2();
		double Q4 = reco_Q2 * reco_Q2;
		double Mott_cross_sec = (1./Q4);

		// ---------------------------------------------------------------------------------------------------------------------

		// Sanity check, especially for radiation
		if (wght < 0 || wght > 10) { std::cout << "Something is really wrong with your weights !!!" << std::endl; }

		double WeightIncl = wght / Mott_cross_sec;

		// Securing ourselves against infinities
		if ( fabs(WeightIncl) != WeightIncl ) { continue; }

		//Calculation of Reconstructed Energy from ELectron only
		//using the same value of single nucleon separation E Ecal and Eqe
		double E_rec = (m_prot*bind_en[ftarget]+m_prot*V4_el.E())/(m_prot-V4_el.E()+V4_el.Rho()*cos(el_theta));
		double EQE_Reso = (E_rec - en_beam_Ecal[fbeam_en]) / en_beam_Ecal[fbeam_en]; 

		//Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
		double nu = -(V4_el-V4_beam).E();
		double x_bjk = reco_Q2/(2*m_prot*nu);

		// QE selection
		//if ( fabs(x_bjk - 1.) > 0.2) { continue; }

		// ---------------------------------------------------------------------------------------------------------------------

		TVector3 V3_q = (V4_beam-V4_el).Vect();
		double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

		//converting theta to degrees
		el_theta = el_theta*TMath::RadToDeg();

		//Cuts on Q2 and W, only keep events with Q2 > Q2cut and W < Wcut
		if ( reco_Q2 < Q2cut || W_var > Wcut) continue;

		// ---------------------------------------------------------------------------------------------------------------------

		TrueElectronsAboveThreshold++;

		//if (!fiducialcut->EFiducialCut(fbeam_en,V3_el) ) continue; // Electron theta & phi fiducial cuts 		

		// electron upper limit for angle
		if (el_theta > 45) { continue; }

		// ---------------------------------------------------------------------------------------------------------------------

		int ElectronSector = el_phi_mod / 60.;

		// ---------------------------------------------------------------------------------------------------------------------

		if (el_theta > MinThetaSlice && el_theta < MaxThetaSlice) {

			h1_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector[ElectronSector]->Fill(nu,WeightIncl/Q4);
			if (Interaction > -1) { h1_InteractionBreakDown_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector[Interaction][ElectronSector]->Fill(nu,WeightIncl/Q4); }

		}
 
		//Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

		//Index variables for hadrons (p and pions)
		int index_p[20]; //index for each proton
		int index_pi[20]; //index for each pion
		int ind_pi_phot[20]; //index for pions and photons
		int index_pipl[20]; //index for each pi plus
		int index_pimi[20]; //index for each pi minus

		int charge_pi[20]; //Charge for the pions and photons
		//Smeared Momentum and Energy values for GENIE (simulation) data
		double Smeared_Pp[20]; //smeared momentum values for protons
		double Smeared_Ep[20]; //smeared energy values for protons
		double Smeared_Ppi[20]; //smeared momentum values for pions
		double Smeared_Epi[20]; //smeared energy values for pions

		//Number of hadrons
		int num_p = 0;
		int num_pi = 0;
		int num_pi_phot = 0; //couting all pions and photons
		int num_phot = 0; //couting all photons
		int num_pimi = 0;
		int num_pipl = 0;
		int num_pi_phot_nonrad = 0; //counting all pions and non-radiation photons
		int num_phot_rad = 0; //counting radiation photons
		int num_phot_nonrad = 0;
		//Index and number variables for neutral particles
		int ec_index_n[20];
		int ec_num_n = 0;
		bool ec_radstat_n[20];

		//Array initialize to -1 or false
		for (int i = 0; i < 20; i++) {
			index_p[i] = -1;   index_pi[i] = -1;   index_pipl[i] = -1;   index_pimi[i] = -1;   ind_pi_phot[i] = -1;
			ec_index_n[i] = -1;   ec_radstat_n[i] = false;
			charge_pi[i] = -2; //default number should be not a possible real charge
			Smeared_Pp[i]  = 0; Smeared_Ep[i]  = 0;  //default 0 momentum and energy after smearing
			Smeared_Ppi[i] = 0; Smeared_Epi[i] = 0;  //default 0 momentum and energy after smearing
		}

		const double phot_rad_cut = 40;
		const double phot_e_phidiffcut = 30; //electron - photon phi difference cut

		// Creating vectors to store id of particles in the array
		vector <int> ProtonID; vector <int> PiPlusID; vector <int> PiMinusID; vector <int> PhotonID;
		ProtonID.clear(); PiPlusID.clear(); PiMinusID.clear();  PhotonID.clear();

		//Loop for Hadrons

		for (int i = 0; i < nf; i++) {

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			//Start of proton selection

			if (pdgf[i] == 2212  && pf[i] > 0.3) {

				//Smearing of proton
				double temp_smear_P = pf[i];
				double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_prot * m_prot );

				TVector3 V3_prot_corr(pxf[i],pyf[i],pzf[i]);
				double phi_prot = V3_prot_corr.Phi();
				//V3_prot_corr.SetPhi(phi_prot + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS

				//if (PFiducialCutExtra(StoreEnergy, V3_prot_corr)) { TrueProtonsAboveThreshold++; }

				double theta_proton = V3_prot_corr.Theta() * 180./TMath::Pi();
				if (theta_proton > 12) { TrueProtonsAboveThreshold++; }

				num_p = num_p + 1;
				index_p[num_p - 1] = i;
				ProtonID.push_back(i);
				Smeared_Pp[num_p - 1] = temp_smear_P;
				Smeared_Ep[num_p - 1] = temp_smear_E;

				//phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

			}

			// -------------------------------------------------------------------------------------------------------------------

			if (pdgf[i] == -211  && pf[i] > 0.15)  { //Pi minus

				//Smearing of pi minus
				double temp_smear_P = pf[i];
				double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

				TVector3 V3_pi_corr(pxf[i],pyf[i],pzf[i]);
				double phi_pion = V3_pi_corr.Phi();
				V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)

				if (PimiFiducialCutExtra(StoreEnergy, V3_pi_corr)) { 

					TrueChargedPionsAboveThreshold++;
					TruePiMinusAboveThreshold++;

				}

				num_pimi = num_pimi + 1;
				num_pi = num_pi + 1;
				num_pi_phot = num_pi_phot + 1;
				num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
				index_pimi[num_pi_phot - 1] = i;
				index_pi[num_pi_phot - 1] = i;
				ind_pi_phot[num_pi_phot - 1] = i;
				PiMinusID.push_back(i);
				charge_pi[num_pi_phot - 1] = -1;
				Smeared_Ppi[num_pi_phot - 1] = temp_smear_P;
				Smeared_Epi[num_pi_phot - 1] = temp_smear_E;

				phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

			}

			// -------------------------------------------------------------------------------------------------------------------

			if ( pdgf[i] == 211  && pf[i] > 0.15)  {

				//Smearing of pi plus
				double temp_smear_P = pf[i];
				double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

				TVector3 V3_pi_corr(pxf[i],pyf[i],pzf[i]);
				double phi_pion = V3_pi_corr.Phi();
				V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)

				if (PiplFiducialCutExtra(StoreEnergy, V3_pi_corr)) { 

					TrueChargedPionsAboveThreshold++;
					TruePiPlusAboveThreshold++;

				}

				num_pipl = num_pipl + 1;
				num_pi  = num_pi + 1;
				num_pi_phot = num_pi_phot + 1;
				num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
				index_pipl[num_pi_phot - 1] = i;
				index_pi[num_pi_phot - 1] = i;
				ind_pi_phot[num_pi_phot - 1] = i;
				PiPlusID.push_back(i);
				charge_pi[num_pi_phot - 1] = 1;
				Smeared_Ppi[num_pi_phot - 1] = temp_smear_P;
				Smeared_Epi[num_pi_phot - 1] = temp_smear_E;

				phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

			}

			// ---------------------------------------------------------------------------------------------------------------------------

			if (pdgf[i] == 22  && pf[i] > 0.3) {

				//Determine photon vector for the cut on radiation photon via angle with respect to the electron
				TVector3 V3_phot_angles(pxf[i],pyf[i],pzf[i]);

				//no smearing of GENIE photons
				double phi_photon = V3_phot_angles.Phi();
				V3_phot_angles.SetPhi(phi_photon + TMath::Pi()); // Vec.Phi() is between (-180,180)

				if (Phot_fidExtra(V3_phot_angles)) { 
	
					TrueGammasAboveThreshold++;

				}

				double neut_phi_mod = V3_phot_angles.Phi()*TMath::RadToDeg() + 30; //Add 30 degree
				if (neut_phi_mod < 0) neut_phi_mod = neut_phi_mod + 360;  //Neutral particle is between 0 and 360 degree

				ec_num_n = ec_num_n + 1;
				num_pi_phot = num_pi_phot + 1;
				ind_pi_phot[num_pi_phot - 1] = i;
				PhotonID.push_back(i);

				Smeared_Ppi[num_pi_phot - 1] = V3_phot_angles.Mag();
				Smeared_Epi[num_pi_phot - 1] = V3_phot_angles.Mag();

				 // within 40 degrees in theta and 30 degrees in phi. Electron phi has already added 30 degree and between 0 to 360

				 if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg() < phot_rad_cut && fabs(neut_phi_mod-el_phi_mod) < phot_e_phidiffcut ) {

					ec_radstat_n[num_pi_phot - 1] = true; // select radiation photons
					num_phot_rad = num_phot_rad + 1;

				 }

				 if(!ec_radstat_n[num_pi_phot - 1]) {

					num_phot_nonrad = num_phot_nonrad + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					charge_pi[num_pi_phot - 1] = 0;

				 }

			}

		} // end of hadron loop

		// ----------------------------------------------------------------------------------------------------------------------------

		// Truth level studies
		// Requiring true level signal 1e 1p 0pi+/- 0 gammas
		// With smearing / fiducial cuts / acceptance maps

		if (TruthLevel1p0piSignalStudy || TruthLevel0piSignalStudy) {

			if (TrueElectronsAboveThreshold != 1) { continue; }
			if (TruthLevel1p0piSignalStudy) { if (TrueProtonsAboveThreshold != 1) { continue; } }
			if (TrueChargedPionsAboveThreshold != 0) { continue; }
			if (TrueGammasAboveThreshold != 0) { continue; }

		}

		// ----------------------------------------------------------------------------------------------------------------------------

		//Skip event if there is at least one radiation photon

		if (num_phot_rad > 0) {

			continue;

		}

		// -----------------------------------------------------------------------------------------------------------------------------

		// Filling Histograms with multiplicities

		h1_Npi->Fill(num_pi);
		h1_Nprot->Fill(num_p);

		if (num_p > 0) {

			h1_Nprot_NonZeroProt->Fill(num_p);
			h1_Npi_NonZeroProt->Fill(num_pi);

		}

		h1_Nphot->Fill(ec_num_n);
		h1_Npipl->Fill(num_pipl);
		h1_Npimi->Fill(num_pimi);
		h1_Npiphot->Fill(num_pi_phot);
		h1_Npiphot_norad->Fill(num_pi_phot_nonrad);

		// ---------------------------------------------------------------------------------------------------------------------------------------

		if (num_pi_phot == 0) {

			h1_E_rec_0pi->Fill(E_rec,WeightIncl);
			h1_E_rec_0pi_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],WeightIncl);

			InclusiveEQE_BreakDown[0]->Fill(E_rec,WeightIncl);
			if (WeightIncl > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,WeightIncl); }

		}

		//------------------------------------------requiring there to be a proton -------------------------------------

		//Events with exactly one proton

		if( num_p == 1) {

			//Vector for proton without momentum smearing
			TLorentzVector V4_prot_uncorr(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
			TVector3 V3_prot_uncorr = V4_prot_uncorr.Vect();

			//Vector for proton with momentum smearing or correction (energy loss)
			TVector3 V3_prot_corr;
			TLorentzVector V4_prot_corr;

			//Fiducial cuts are done in the hadron loop
			//Vector for proton with momentum smearing
			V3_prot_corr.SetXYZ(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]]);
			V4_prot_corr.SetPxPyPzE(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],Smeared_Ep[0]);

			double phi_prot = V3_prot_corr.Phi(); //in Radians
			V3_prot_corr.SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
			phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

			//Proton kinematic variables
			double p_theta = V3_prot_corr.Theta();
			double prot_mom_corr = V3_prot_corr.Mag();

			TLorentzVector V4_prot_el_tot = V4_prot_corr + V4_el;

			double p_perp_tot = TMath::Sqrt(V4_prot_el_tot.Px()*V4_prot_el_tot.Px() + V4_prot_el_tot.Py()*V4_prot_el_tot.Py());
			double E_tot = V4_el.E() + V4_prot_corr.E() - m_prot + bind_en[ftarget];

			//---------------------------------- 1p 0pi   ----------------------------------------------

			// Main Plots

			if(num_pi_phot == 0){

				double ECalReso = (E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en];
				double EQEReso = (E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en];

				SignalEvents++;

				if (Interaction == 1) { QESignalEvents++; }
				else if (Interaction == 2) { MECSignalEvents++; }
				else if (Interaction == 3) { RESSignalEvents++; }
				else if (Interaction == 4) { DISSignalEvents++; }
				else { OtherSignalEvents++; }

				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = wght/Mott_cross_sec;

				h1_E_tot_cut2->Fill(E_tot,histoweight);
				h1_E_tot_cut2_fracfeed->Fill(ECalReso,histoweight);

				double ProtonTheta_Deg = V3_prot_corr.Theta() *180. / TMath::Pi();
				double ProtonMag = V3_prot_corr.Mag();

				h1_MissMomentum->Fill(p_perp_tot,histoweight);

				double LocalWeight = histoweight;

				ECal_BreakDown[0]->Fill(E_tot,LocalWeight);
				EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
				Pmiss_BreakDown[0]->Fill(p_perp_tot,LocalWeight);
				
				double deltaphiT = DeltaPhiTFunction(V3_el,V3_prot_corr);		
				double deltaalphaT = DeltaAlphaTFunction(V3_el,V3_prot_corr);

				double PTmiss = p_perp_tot;
				double Ecal = E_tot;
				double EcalReso = (Ecal-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en];
													     
				DeltaPhiT_BreakDown[0]->Fill(deltaphiT,LocalWeight);	     
				DeltaAlphaT_BreakDown[0]->Fill(deltaalphaT,LocalWeight);				

				ECal_BreakDown[Interaction]->Fill(E_tot,LocalWeight);
				EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
				Pmiss_BreakDown[Interaction]->Fill(p_perp_tot,LocalWeight);				
				DeltaPhiT_BreakDown[Interaction]->Fill(deltaphiT,LocalWeight);	     
				DeltaAlphaT_BreakDown[Interaction]->Fill(deltaalphaT,LocalWeight);					

				if (p_perp_tot < pperp_max[0]) { 

					ECal_LowPmiss_BreakDown[Interaction]->Fill(E_tot,LocalWeight); 
					EQE_LowPmiss_BreakDown[Interaction]->Fill(E_rec,LocalWeight); 

				}

				if (p_perp_tot > pperp_min[1] && p_perp_tot < pperp_max[1]) { 

					ECal_MidPmiss_BreakDown[Interaction]->Fill(E_tot,LocalWeight); 
					EQE_MidPmiss_BreakDown[Interaction]->Fill(E_rec,LocalWeight); 

				}

				if (p_perp_tot > pperp_min[2] && p_perp_tot < pperp_max[2]) { 

					ECal_HighPmiss_BreakDown[Interaction]->Fill(E_tot,LocalWeight); 
					EQE_HighPmiss_BreakDown[Interaction]->Fill(E_rec,LocalWeight); 

				}

				// -----------------------------------------------------------------------------------------------

				for(int i = 0; i < n_slice; i++) {

					if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){

						h1_Etot_Npi0[i]->Fill(E_tot,histoweight);
						h1_Erec_Npi0[i]->Fill(E_rec,histoweight);

					}
				}

			} //num pi=0

		} // End of 1-proton case

	} //end of event loop (jentry)

	gStyle->SetOptFit(1);
	gDirectory->Write("hist_Files", TObject::kOverwrite);

	// --------------------------------------------------------------------------------------------------------

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "# Processed Events = " << TotalCounter << std::endl;
	std::cout << std::endl << "1e1p0pi Signal # Events = " << SignalEvents << std::endl;
	std::cout << std::endl << "Passing Rate = " << int(double(SignalEvents) / double(TotalCounter)*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "QE Fractional Contribution = " << int(double(QESignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
	std::cout << std::endl << "MEC Fractional Contribution = " << int(double(MECSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
	std::cout << std::endl << "RES Fractional Contribution = " << int(double(RESSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
	std::cout << std::endl << "DIS Fractional Contribution = " << int(double(DISSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;

	std::cout << "File " << FileName << " created" << std::endl << std::endl;

} // End of program

//End Loop function
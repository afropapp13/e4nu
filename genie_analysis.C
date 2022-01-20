#define GENIE_ANALYSIS_C

#include "genie_analysis.h"
#include "Constants.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TProfile.h>
#include <TH1D.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TMath.h>
#include <exception>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH3.h>
#include <TGraph.h>

#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;

//----------------------------------------//

// Loading all the constants from Constant.h (e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565

void genie_analysis::Loop(Int_t choice) {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	int NSectors = 6;

	// All Interactions = 0, QE = 1, MEC = 2, RES = 3, DIS = 4, COH = 5
	const int NInt = 6;

	//----------------------------------------//

	// Do we want to apply fiducials & the acceptance map weights
	// Do we want a truth level study ? if so, stop ditching sectors

	bool UseAllSectors = false;
	bool ApplyFiducials = true;
	bool ApplyAccWeights = true;
	bool ApplyReso = true;
	
	bool ApplyPhiOpeningAngle = true;
	bool ApplyThetaSlice = true;

	//----------------------------------------//

	//Choice = 0 is for analysis of CLAS data while choice = 1 is for the analysis of GENIE Simulation
	if (choice != 4 && choice != 3 && choice != 2 && choice != 1 && choice != 0) {
		std::cout << "Abort with fchoice. The given value is " << choice << std::endl;
		std::exit(0);
	}

	target_name = ftarget; //std string for target name
	en_beam["1161"]=1.161;
	en_beam["2261"]=2.261;
	en_beam["4461"]=4.461;

	//----------------------------------------//

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();

	//----------------------------------------//

	double reso_e = 0.005; // smearing for the electrons

	if (!ApplyReso) {  reso_e = 0.; }

	// Resolution defined above seems to be insufficient at 1.1 GeV -> tripled it for all particles
	if(fbeam_en == "1161") { reso_e = 3*reso_e; }

	//----------------------------------------//	

	//cut for all beam energies < 2
	double Wcut = 2;

	// cut for 1.1 GeV > 0.1, for 2.2 GeV > 0.4 and 4.4 GeV > 0.8
	double Q2cut = 0; 

	//----------------------------------------//	

	TString E_acc_file;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) //1.1 GeV Configuration parameters and cuts
	{
		E_acc_file="1_161";
		Q2cut = 0.1;
	}


	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) //2.2 GeV Configuration parameters and cuts
	{
		E_acc_file="2_261";
		Q2cut = 0.4;
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) //4.4 GeV Configuration parameters and cuts
	{
		E_acc_file="4_461";
		Q2cut = 0.8;
	}

	//----------------------------------------//

	gRandom = new TRandom3();
	gRandom->SetSeed(10);

	//----------------------------------------//	

	TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);

	//----------------------------------------//	

	//Acceptance Maps

	TString WhichMap = "e2a_maps";
	TFile* file_acceptance;

	TString Target = "12C";
	if (ftarget.c_str() == "3He") { Target = "3He"; }
	if (ftarget.c_str() == "4He") { Target = "4He"; }

	if ( choice > 0 ) { // Only need acceptance maps for GENIE simulation

		file_acceptance = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+".root");
	}

	//----------------------------------------//

	//Output file definition

	TFile *file_out;
	TString FileName = ""; 

	if (choice == 0) { FileName = Form("/work/clas/claseg2/apapadop/Inclusive_data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()); }
	if (choice == 1){ FileName = Form("Inclusive_genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test_SuSav2.root",ftarget.c_str(),fbeam_en.c_str()); }
	if (choice == 2) { FileName = Form("Inclusive_genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test_G18_10a_02_11a.root",ftarget.c_str(),fbeam_en.c_str()); }
	if (choice == 3) { FileName = Form("Inclusive_genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test_SuSav2_Rad.root",ftarget.c_str(),fbeam_en.c_str()); }
	if (choice == 4) { FileName = Form("Inclusive_genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test_G18_10a_02_11a_Rad.root",ftarget.c_str(),fbeam_en.c_str()); }

	file_out = new TFile(FileName, "Recreate");

	//----------------------------------------//

	//initialize Fiducial functions for EC limits
	fiducialcut->InitEClimits();
	std::cout << " Test InitEClimits Loop " << fiducialcut->up_lim1_ec->Eval(60) << std::endl;

	//----------------------------------------//	

	//Definition and initialization of Histograms

	TH1F *h1_Omega[ThetaSlices[fbeam_en]][NSectors];
	TH1F *h1_BreakDown_Omega[ThetaSlices[fbeam_en]][NInt][NSectors];

	for (int WhichSector = 0; WhichSector < NSectors; WhichSector++) {

		for (int WhichTheta = 0; WhichTheta < ThetaSlices[fbeam_en]; WhichTheta++) {

			double ThetaStep = (MaxThetaSlice[fbeam_en] - MinThetaSlice[fbeam_en] ) / ThetaSlices[fbeam_en];
			double MinTheta = MinThetaSlice[fbeam_en] + WhichTheta*ThetaStep;
			double MaxTheta = MinThetaSlice[fbeam_en] + (WhichTheta+1)*ThetaStep;
			int AveTheta = (MaxTheta + MinTheta) / 2.;		
			TString StringAveTheta = TString(std::to_string(AveTheta));
			TString StringSector = TString(std::to_string(WhichSector));

			h1_Omega[WhichTheta][WhichSector]  = new TH1F("h1_Omega_Theta_"+StringAveTheta+"_Sector_"+StringSector,"",6000,0.,6.);

			for (int WhichInt = 1; WhichInt < NInt; WhichInt++) {

				TString StringInt = TString(std::to_string(WhichInt));

				h1_BreakDown_Omega[WhichTheta][WhichInt][WhichSector]  = new TH1F("h1_"+StringInt+"_Omega_Theta_"+StringAveTheta+"_Sector_"+StringSector,"",6000,0.,6.);

			}

		}

	}

	//----------------------------------------//

	// Counters

	int TotalCounter = 0;	
	int SignalEvents = 0;
	int QESignalEvents = 0;
	int MECSignalEvents = 0;
	int RESSignalEvents = 0;
	int DISSignalEvents = 0;

	//----------------------------------------//

	// Justification for the parameter choice
	// https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing

	TF1 *myElectronFit = new TF1("myElectronFit","[0]+[1]/x",0.,5.);

	if (en_beam[fbeam_en] == 1.161) { myElectronFit->SetParameters(17,7); }
	if (en_beam[fbeam_en] == 2.261) { myElectronFit->SetParameters(16,10.5); }
	if (en_beam[fbeam_en] == 4.461) { myElectronFit->SetParameters(13.5,15); }

	//----------------------------------------//

	// Loop over the events	

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

		//----------------------------------------//

		std::string StoreEnergy = fbeam_en;

		//----------------------------------------//

		// For GENIE samples, identify the interaction type

		int Interaction = -1;

		if (choice > 0) {

			SignalEvents++;
			if (qel) { Interaction = 1; QESignalEvents++; }
			if (mec) { Interaction = 2; MECSignalEvents++; }
			if (res) { Interaction = 3; RESSignalEvents++; }
			if (dis) { Interaction = 4; DISSignalEvents++; }

		}

		//----------------------------------------//

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
			std::cout << " EventLoop: Finished setting up fiducial cut class " << std::endl;
		}

		//----------------------------------------//

		double SmearedPe;
		double SmearedEe;
		double e_acc_ratio = 1.;	//will be 1 for CLAS data

		// Outgoing e',	Uncorr and corrected are the same read from root file.
		//V4_el and V3_el will be changed by smearing for GENIE simulation data
		TLorentzVector V4_el(pxl,pyl,pzl,El);
		TLorentzVector V4_el_uncorr(pxl,pyl,pzl,El);
		TVector3 V3_el(pxl,pyl,pzl);

		double el_momentum = V3_el.Mag();
		double el_theta = V3_el.Theta();

		//----------------------------------------//

		if (choice > 0) { //smearing, fiducials and acceptance ratio for GENIE simulation data

			//Smearing of Electron Vector from Simulation
			SmearedPe = gRandom->Gaus(pl,reso_e*pl);
			SmearedEe = sqrt( SmearedPe*SmearedPe + e_mass * e_mass );
			V3_el.SetXYZ(SmearedPe/pl * pxl,SmearedPe/pl * pyl,SmearedPe/pl * pzl);
			V4_el.SetPxPyPzE(V3_el.X(),V3_el.Y(),V3_el.Z(),SmearedEe);
			double phi_ElectronOut = V3_el.Phi(); //in Radians

			V3_el.SetPhi(phi_ElectronOut + TMath::Pi() ); // Vec.Phi() is between (-180,180), GENIE coordinate system flipped with respect to CLAS

			phi_ElectronOut += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
			el_momentum = V3_el.Mag(); //Momentum after smearing
			el_theta = V3_el.Theta(); //Angle after smearing

			// acceptance_c takes phi in radians and here unmodified by 30 degree.

			e_acc_ratio = acceptance_c(el_momentum, cos(el_theta), phi_ElectronOut, 11,file_acceptance,ApplyAccWeights);
			if ( fabs(e_acc_ratio) != e_acc_ratio ) { continue; }

		}

		//----------------------------------------//

		double theta_min = myElectronFit->Eval(el_momentum); // in deg
		if (el_theta*180./TMath::Pi() < theta_min) { continue; }

		double el_theta_deg = el_theta*180./TMath::Pi();

		int ThetaIndex = (el_theta_deg - MinThetaSlice[fbeam_en] ) / ThetaSlices[fbeam_en];

		// ----------------------------------------------------------------------------------------------------------------------

		// Explicit cuts on electron momentum

		if (fbeam_en=="1161" && el_momentum < 0.4) { continue; }
		if (fbeam_en=="2261" && el_momentum < 0.55) { continue; }
		if (fbeam_en=="4461" && el_momentum < 1.1) { continue; }

		//Definition as for data. It is also correct for GENIE simulation data since V3_el is rotated above by 180 degree in phi
		double el_phi_mod = V3_el.Phi()*TMath::RadToDeg()  + 30; //Add 30 degree for plotting and photon phi cut
		if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree

		if (ApplyPhiOpeningAngle) { if ( !(TMath::Abs(el_phi_mod - 30)  < PhiOpeningAngle || TMath::Abs(el_phi_mod - 90)  < PhiOpeningAngle || TMath::Abs(el_phi_mod - 150)  < PhiOpeningAngle || TMath::Abs(el_phi_mod - 210)  < PhiOpeningAngle || TMath::Abs(el_phi_mod - 270)  < PhiOpeningAngle || TMath::Abs(el_phi_mod - 330)  < PhiOpeningAngle ) ) { continue; } }

		//Calculated Mott Cross Section and Weights for Inclusive Histograms
		//Wght and e_acc_ratio is 1 for CLAS data

		double reco_Q2 = -(V4_el-V4_beam).Mag2();

		// ---------------------------------------------------------------------------------------------------------------------

		// Sanity check, especially for radiation
		if (wght < 0 || wght > 10) { std::cout << "Something is really wrong with your weights !!!" << std::endl; }

		double WeightIncl = wght * e_acc_ratio;

		// Securing ourselves against infinities
		if ( fabs(WeightIncl) != WeightIncl ) { continue; }

		//Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
		double nu = -(V4_el-V4_beam).E();

		// ---------------------------------------------------------------------------------------------------------------------

		TVector3 V3_q = (V4_beam-V4_el).Vect();
		double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

		//converting theta to degrees
		el_theta = el_theta*TMath::RadToDeg();

		//Cuts on Q2 and W, only keep events with Q2 > Q2cut and W < Wcut
		if ( reco_Q2 < Q2cut || W_var > Wcut) continue;

		// ---------------------------------------------------------------------------------------------------------------------

		if (choice > 0) { 			

			//Fiducial Cuts with the smeared values
			if (ApplyFiducials)  { if (!EFiducialCut(fbeam_en,V3_el) ) continue;} // Electron theta & phi fiducial cuts 
		}

		int ElectronSector = el_phi_mod / 60.;

		// ---------------------------------------------------------------------------------------------------------------------

		// apapadop: Oct 8 2020: ditching bad sectors
		// Counting sectors from 0 to 5

		if (!UseAllSectors) { 

			if ( (ElectronSector == 2 || ElectronSector == 4) &&  fbeam_en == "1161") { continue; }
			if ( (ElectronSector == 2 || ElectronSector == 3 || ElectronSector == 4) &&  fbeam_en == "2261") { continue; }

		}

		// ---------------------------------------------------------------------------------------------------------------------

		// Fully inclusive energy transfer plots & breakdown

		if (el_theta > MinThetaSlice[fbeam_en] && el_theta < MaxThetaSlice[fbeam_en]) {

			h1_Omega[ThetaIndex][ElectronSector]->Fill(nu,WeightIncl);
			if (Interaction > -1) { h1_BreakDown_Omega[ThetaIndex][Interaction][ElectronSector]->Fill(nu,WeightIncl); }

		}

		// ---------------------------------------------------------------------------------------------------------------------------------------

	} //end of event loop (jentry)

	if (choice > 0) {

		std::cout << std::endl << "QE Fractional Contribution = " << int(double(QESignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "MEC Fractional Contribution = " << int(double(MECSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "RES Fractional Contribution = " << int(double(RESSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "DIS Fractional Contribution = " << int(double(DISSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;

	}

	std::cout << "File " << FileName << " created" << std::endl << std::endl;

} // End of program

//End Loop function

// -------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------

double genie_analysis::acceptance_c(double p, double cost, double phi, int particle_id,TFile* file_acceptance, bool ApplyAccWeights) {

	if (ApplyAccWeights) {

		//Redefinition of the phi angle
		//because the acceptance maps are defined between (-30,330)

		// Check that phi is between (0,360)

		int redef = 0;

		TH3D * acc;
		TH3D * gen;

		acc = (TH3D*)file_acceptance->Get("Accepted Particles");
		gen = (TH3D*)file_acceptance->Get("Generated Particles");

		//map 330 till 360 to [-30:0] for the acceptance map histogram
		if(phi > (2*TMath::Pi() - TMath::Pi()/6.) ) { phi -= 2*TMath::Pi(); }
		//Find number of generated events

		double pbin_gen = gen->GetXaxis()->FindBin(p);
		double tbin_gen = gen->GetYaxis()->FindBin(cost);
		double phibin_gen = gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double num_gen = gen->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

		//Find number of accepted events

		double pbin_acc = acc->GetXaxis()->FindBin(p);
		double tbin_acc = acc->GetYaxis()->FindBin(cost);
		double phibin_acc = acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double num_acc = acc->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

		double acc_ratio = (double)num_acc / (double)num_gen;
		double acc_err = (double)sqrt(acc_ratio*(1-acc_ratio)) / (double)num_gen;


		return acc_ratio;

	}

	else { return 1.; }

}
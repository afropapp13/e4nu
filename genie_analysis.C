#define GENIE_ANALYSIS_C
#include "genie_analysis.h"
#include "Constants.h"
#include <TStyle.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH3.h>

#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>

using namespace std;

//----------------------------------------//

int ReturnIndex(double value, std::vector<double> vec) {

	int length = vec.size();
	int index = -1;

	for (int i = 0; i < length-1; i ++) {

		if (value >= vec.at(i) && value < vec.at(i+1)) { return i; }

	}	

	cout << "value " << value << " results in a -1 index !!! About to crash !!!" << endl;
	return index;

}

//----------------------------------------//

TString to_string_with_precision(double a_value, const int n = 2) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());
}

//----------------------------------------//

TString ConvertToString(double value) {

	TString StringValue = to_string_with_precision(value, 2);
	StringValue.ReplaceAll(".","_");
	StringValue.ReplaceAll("-","Minus");	

	return StringValue;

}

//----------------------------------------//

void STV_Tools(TVector3 MuonVector,TVector3 ProtonVector, double MuonEnergy, double ProtonEnergy, double STLV[],double P, double N, double BindE, double EE, double TrueBeamEnergy) {

//----------------------------------------//

	// 0: Pt
	// 1: DeltaAlphaT
	// 2: DeltaPhiT
	// 3: ECal
	// 4: EQE
	// 5: Q2
	// 6: Emiss
	// 7: Pmiss
	// 8: PmissMinus
	// 9: kMiss
	// 10:alpha
	// 11: PL
	// 12: Pn,proxy
	// 13: Pt,x
	// 14: Pt,y 	
	// 15: PL using PMiss projection
	// 16: Enu with QE assumption

	//----------------------------------------//

	double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938272, NeutronMass_GeV = 0.939565; // GeV
	double DeltaM2 = TMath::Power(NeutronMass_GeV,2.) - TMath::Power(ProtonMass_GeV,2.);	
	double BE = 0.04; // GeV	
			
	// STV Calculation		
			
	TVector3 MuonVectorTrans;
	MuonVectorTrans.SetXYZ(MuonVector.X(),MuonVector.Y(),0.);
	double MuonVectorTransMag = MuonVectorTrans.Mag();

	TVector3 MuonVectorLong;
	MuonVectorLong.SetXYZ(0.,0.,MuonVector.Z());
	double MuonVectorLongMag = MuonVectorLong.Mag();
	
	TLorentzVector MuonLorentzVector(MuonVector,MuonEnergy);
			
	TVector3 ProtonVectorTrans;
	ProtonVectorTrans.SetXYZ(ProtonVector.X(),ProtonVector.Y(),0.);
	double ProtonVectorTransMag = ProtonVectorTrans.Mag();

	TVector3 ProtonVectorLong;
	ProtonVectorLong.SetXYZ(0.,0.,ProtonVector.Z());
	double ProtonVectorLongMag = ProtonVectorLong.Mag();
	
	TLorentzVector ProtonLorentzVector(ProtonVector,ProtonEnergy);
	double ProtonKE = ProtonEnergy - ProtonMass_GeV;		

	TVector3 PtVector = MuonVectorTrans + ProtonVectorTrans;

	STLV[0] = PtVector.Mag();

	STLV[1] = TMath::ACos( (- MuonVectorTrans * PtVector) / ( MuonVectorTransMag * STLV[0] ) ) * 180./TMath::Pi();
	if (STLV[1] > 180.) { STLV[1] -= 180.; }
	if (STLV[1] < 0.) { STLV[1] += 180.; }


	STLV[2] = TMath::ACos( (- MuonVectorTrans * ProtonVectorTrans) / ( MuonVectorTransMag * ProtonVectorTransMag ) ) * 180./TMath::Pi();
	if (STLV[2] > 180.) { STLV[2] -= 180.; }
	if (STLV[2] < 0.) { STLV[2] += 180.; }

	//----------------------------------------//

	// Calorimetric Energy Reconstruction

	STLV[3] = MuonEnergy + ProtonKE + BE; // ECal,GeV

	// QE Energy Reconstruction

	double EQENum = 2 * (NeutronMass_GeV - BE) * MuonEnergy - (BE*BE - 2 * NeutronMass_GeV *BE + MuonMass_GeV * MuonMass_GeV + DeltaM2);
	double EQEDen = 2 * ( NeutronMass_GeV - BE - MuonEnergy + MuonVector.Mag() * MuonVector.CosTheta() );
	STLV[4] = EQENum / EQEDen; // EQE, GeV

	// Reconstructed Q2

	// TLorentzVector nuLorentzVector(0.,0.,STLV[3],STLV[3]); // neutrinos
	TLorentzVector nuLorentzVector(0.,0.,TrueBeamEnergy,TrueBeamEnergy); // electrons, known incoming energy
	TLorentzVector qLorentzVector = nuLorentzVector - MuonLorentzVector;
	STLV[5] = - qLorentzVector.Mag2(); // Q2, GeV^{2}/c^{2}
	
	//----------------------------------------//
	
	// Light Cone Variables
	
	// TLorentzVector MissLorentzVector = MuonLorentzVector + ProtonLorentzVector - nuLorentzVector; // neutrinos
	TLorentzVector MissLorentzVector = ProtonLorentzVector - qLorentzVector; // electrons
	
	STLV[6] = TMath::Abs(MissLorentzVector.E());  // Emiss
	STLV[7] = (MissLorentzVector.Vect()).Mag();   // Pmiss
	STLV[15] = MissLorentzVector.Z(); // PL using PMiss projection along z direction
	
//	fPMissMinus = fEMiss - MissLorentzVector.Z();

	// Suggestion from Jackson to avoid Ecal assumption
	STLV[8] = (MuonEnergy - MuonVector.Z()) + (ProtonEnergy - ProtonVector.Z()); // PmissMinus
	
	double fkMissNum = ( TMath::Power(STLV[0],2.) + TMath::Power(ProtonMass_GeV,2.) );
	double fkMissDen = ( STLV[8] * (2*ProtonMass_GeV - STLV[8]) );
	
	double fkMiss2 = TMath::Power(ProtonMass_GeV,2.) * fkMissNum / fkMissDen - TMath::Power(ProtonMass_GeV,2.); // Jackson's GlueX note				

	STLV[9] = sqrt(fkMiss2); // kMiss

	STLV[10] = STLV[8] / ProtonMass_GeV; // alpha

	//----------------------------------------//

	// Minerva longitudinal & total variables

	// For the calcukation of the masses
	//https://journals.aps.org/prc/pdf/10.1103/PhysRevC.95.065501

	double MA = N * NeutronMass_GeV + P * ProtonMass_GeV - BindE; // GeV
	double MAPrime = MA - NeutronMass_GeV + EE; // GeV, constant obtained from table I 

	// For the calculation of p_n, back to the Minerva PRL
	// https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.022504

	double R = MA + MuonVectorLongMag + ProtonVectorLongMag - MuonEnergy - ProtonEnergy; // Equation 8

	// Equation 7

	STLV[11] = 0.5 * R - (MAPrime * MAPrime + STLV[0] * STLV[0]) / (2 * R); // PL

	STLV[12] = TMath::Sqrt( STLV[0] * STLV[0] + STLV[11] * STLV[11] ); // Pn,proxy

	TVector3 UnitZ(0,0,1);
	STLV[13] = ( UnitZ.Cross(MuonVectorTrans) ).Dot(PtVector) / MuonVectorTransMag;// Pt,x
	STLV[14] = - (MuonVectorTrans).Dot(PtVector) / MuonVectorTransMag;// Pt,y

	STLV[16] = MuonVectorLongMag + ProtonVectorLongMag- STLV[11]; // Enu using the calorimetric assumption

	//----------------------------------------//

}

//----------------------------------------//

// Loading all the constants from Constant.h (e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565,
// H3_bind_en, He4_bind_en, C12_bind_en, B_bind_en, He3_bind_en, D2_bind_en, Fe_bind_en, Mn_bind_en

void genie_analysis::Loop(Int_t choice) {

	TH1D::SetDefaultSumw2();

	//----------------------------------------//

	// Do we want to apply fiducials & the acceptance map weights
	// Do we want a truth level study ? if so, stop ditching sectors

	bool UseAllSectors = false;
	bool ApplyFiducials = true;
	bool ApplyAccWeights = true;
	bool ApplyReso = true;
	bool TruthLevel1p0piSignalStudy = false;

	//----------------------------------------//

	//Choice = 0 is for analysis of CLAS data while choice = 1 is for the analysis of GENIE Simulation
	if (choice != 4 && choice != 3 && choice != 2 && choice != 1 && choice != 0) {
		std::cout << "This parameter value is not implemented in genie_analysis::Loop(). It should be either 0 or 1. The given value is " << choice << std::endl;
		std::exit(0);
	}

	std::map<std::string,double> bind_en;
	std::map<std::string,double> target_mass;
	std::map<std::string,double> residual_target_mass;
	std::map<std::string, double> Ecal_offset; // that might not be necessary for simulation data

	target_name = ftarget; //std string for target name
	en_beam["1161"]=1.161;
	en_beam["2261"]=2.261;
	en_beam["4461"]=4.461;

	//----------------------------------------//

	// Protons in target

	std::map<std::string,double> NProtons;	
	NProtons["3He"]  = 2.;
	NProtons["4He"]  = 2.;
	NProtons["C12"]  = 6.;
	NProtons["56Fe"]  = 26.;
	NProtons["40Ar"]  = 18.;

	// Neutrons in target

	std::map<std::string,double> NNeutrons;	
	NNeutrons["3He"]  = 1.;
	NNeutrons["4He"]  = 2.;
	NNeutrons["C12"]  = 6.;
	NNeutrons["56Fe"]  = 30.;
	NNeutrons["40Ar"]  = 22.;

	// Binding Energy

	std::map<std::string,double> BindE;	
	BindE["3He"]  = 0.02316;
	BindE["4He"]  = 0.02829;
	BindE["C12"]  = 0.09216;
	BindE["56Fe"]  = 0.4922;
	BindE["40Ar"]  = 0.34381;	

	// Excitation Energy

	std::map<std::string,double> EE;	
	EE["3He"]  = 0.0042;
	EE["4He"]  = 0.0083;
	EE["C12"]  = 0.0261;
	EE["56Fe"]  = 0.0191;

	//----------------------------------------//

	en_beam_Ecal["1161"]=1.161;
	en_beam_Ecal["2261"]=2.261;
	en_beam_Ecal["4461"]=4.461;

	en_beam_Eqe["1161"]=1.161;
	en_beam_Eqe["2261"]=2.261;
	en_beam_Eqe["4461"]=4.461;

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	//Resolutions for Smearing for GENIE simulation data
	double reso_p = 0.01; // smearing for the proton
	double reso_e = 0.005; // smearing for the electrons
	double reso_pi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)

	if (!ApplyReso) {  

		reso_p = 0.;
		reso_e = 0.;
		reso_pi = 0;

	}

	// Resolution defined above seems to be insufficient at 1.1 GeV -> tripled it for all particles
	if(fbeam_en == "1161") { reso_p = 3*reso_p; reso_e = 3*reso_e; reso_pi = 3*reso_pi; }

	double Wcut = 2; //cut for all beam energies < 2
	double Q2cut = 0; // cut for 1.1 GeV > 0.1, for 2.2 GeV > 0.4 and 4.4 GeV > 0.8

	TVector3 V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;

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

	//----------------------------------------//

	gRandom = new TRandom3();
	gRandom->SetSeed(10);

	TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);
	TLorentzVector V4_target(0,0,0,target_mass[ftarget]);

	//Acceptance Maps

	TString WhichMap = "e2a_maps";
	TFile* file_acceptance = nullptr;
	TFile* file_acceptance_p = nullptr;
	TFile* file_acceptance_pip = nullptr;
	TFile* file_acceptance_pim = nullptr;

	TString Target = "12C";
	if ( TString(ftarget.c_str()) == "3He") { Target = "3He"; }
	if ( TString(ftarget.c_str()) == "4He") { Target = "4He"; }

	if ( choice > 0 ) { // Only need acceptance maps for GENIE simulation

		file_acceptance = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+".root");
		file_acceptance_p = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_p.root");
		file_acceptance_pip = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_pip.root");
		file_acceptance_pim = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_pim.root");
	}

	//----------------------------------------//

	//Output file definition
	TString FileName = ""; 

	if (choice == 0) { FileName = Form("/work/clas/claseg2/apapadop/MultiDim_data_e2a_ep_%s_%s.root",ftarget.c_str(),fbeam_en.c_str()); }
	if (choice == 1){ FileName = Form("MultiDim_genie_e2a_ep_%s_%s_SuSav2.root",ftarget.c_str(),fbeam_en.c_str()); }
	if (choice == 2) { FileName = Form("MultiDim_genie_e2a_ep_%s_%s_G18_10a_02_11a.root",ftarget.c_str(),fbeam_en.c_str()); }
	if (choice == 3) { FileName = Form("MultiDim_genie_e2a_ep_%s_%s_SuSav2_Rad.root",ftarget.c_str(),fbeam_en.c_str()); }
	if (choice == 4) { FileName = Form("MultiDim_genie_e2a_ep_%s_%s_G18_10a_02_11a_Rad.root",ftarget.c_str(),fbeam_en.c_str()); }

	TFile* file_out = new TFile(FileName, "recreate");

	//----------------------------------------//

	fiducialcut->InitPiMinusFit(fbeam_en);

	//initialize Fiducial functions for EC limits
	fiducialcut->InitEClimits();
	std::cout << " Test InitEClimits Loop " << fiducialcut->up_lim1_ec->Eval(60) << std::endl;

	//----------------------------------------//

	// KI multi dim plots
	// 1st index: 1,2... = slices
	// 2nd index: 0 = all events, QE = 1, MEC = 2, RES = 3, DIS = 4, COH = 5

	TH1D* DeltaPT_InDeltaAlphaTPlot[int(TwoDArrayNBinsDeltaAlphaT.size()) -1][NInt];
	TH1D* DeltaAlphaT_InDeltaPTPlot[int(TwoDArrayNBinsDeltaPT.size()) -1][NInt];
	TH1D* DeltaPtx_InDeltaPtyPlot[int(TwoDArrayNBinsDeltaPty.size()) -1][NInt];
	TH1D* DeltaPty_InDeltaPtxPlot[int(TwoDArrayNBinsDeltaPtx.size()) -1][NInt];			

	// Loop over the interactions
	for (int iinte = 0; iinte < NInt; iinte++) {

		// Loop over the DeltaAlphaT slices
		for (int ideltaalphat = 0; ideltaalphat < int(TwoDArrayNBinsDeltaAlphaT.size() -1);ideltaalphat++ ) {

			TString DeltaPTTwoDInDeltaAlphaTLabel = InteractionLabels[iinte]+"DeltaPT_InDeltaAlphaT_Slice_"+TString( std::to_string(ideltaalphat) )+"Plot";			
			DeltaPT_InDeltaAlphaTPlot[ideltaalphat][iinte] = new TH1D(DeltaPTTwoDInDeltaAlphaTLabel,";#deltap_{T} [GeV/c]",TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[ideltaalphat].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[ideltaalphat][0]);			

		}

		// Loop over the DeltaPT slices
		for (int ideltapt = 0; ideltapt < int(TwoDArrayNBinsDeltaPT.size() -1);ideltapt++ ) {

			TString DeltaAlphaTTwoDInDeltaPTLabel = InteractionLabels[iinte]+"DeltaAlphaT_InDeltaPT_Slice_"+TString( std::to_string(ideltapt) )+"Plot";			
			DeltaAlphaT_InDeltaPTPlot[ideltapt][iinte] = new TH1D(DeltaAlphaTTwoDInDeltaPTLabel,";#delta#alpha_{T} [deg]",TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[ideltapt].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[ideltapt][0]);			

		}	

		// Loop over the DeltaPty slices
		for (int ideltapty = 0; ideltapty < int(TwoDArrayNBinsDeltaPty.size() -1);ideltapty++ ) {

			TString DeltaPtxTwoDInDeltaPtyLabel = InteractionLabels[iinte]+"DeltaPtx_InDeltaPty_Slice_"+TString( std::to_string(ideltapty) )+"Plot";			
			DeltaPtx_InDeltaPtyPlot[ideltapty][iinte] = new TH1D(DeltaPtxTwoDInDeltaPtyLabel,";#deltap_{T,x} [GeV/c]",TwoDArrayNBinsDeltaPtxInDeltaPtySlices[ideltapty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[ideltapty][0]);			

		}	

		// Loop over the DeltaPtx slices
		for (int ideltaptx = 0; ideltaptx < int(TwoDArrayNBinsDeltaPtx.size() -1);ideltaptx++ ) {

			TString DeltaPtyTwoDInDeltaPtxLabel = InteractionLabels[iinte]+"DeltaPty_InDeltaPtx_Slice_"+TString( std::to_string(ideltaptx) ) +"Plot";			
			DeltaPty_InDeltaPtxPlot[ideltaptx][iinte] = new TH1D(DeltaPtyTwoDInDeltaPtxLabel,";#deltap_{T,y} [GeV/c]",TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[ideltaptx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[ideltaptx][0]);			

		}				


	} // End of the loop over the interactions				

	//----------------------------------------//

	int TotalCounter = 0;
	int SignalEvents = 0;
	int QESignalEvents = 0;
	int MECSignalEvents = 0;
	int RESSignalEvents = 0;
	int DISSignalEvents = 0;
	int OtherSignalEvents = 0;

	//----------------------------------------//

	// Justification for the parameter choice
	// https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing

	TF1 *myElectronFit = new TF1("myElectronFit","[0]+[1]/x",0.,5.);

	if (en_beam[fbeam_en] == 1.161) { myElectronFit->SetParameters(17,7); }
	if (en_beam[fbeam_en] == 2.261) { myElectronFit->SetParameters(16,10.5); }
	if (en_beam[fbeam_en] == 4.461) { myElectronFit->SetParameters(13.5,15); }

	//----------------------------------------//
	//----------------------------------------//

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

			if (qel) { Interaction = 1; }
			if (mec) { Interaction = 2; }
			if (res) { Interaction = 3; }
			if (dis) { Interaction = 4; }

		}

	//----------------------------------------//

		if (jentry == 0) { //first entry to initialize TorusCurrent, Fiducials and Subtraction classes

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
//			rotation->InitSubtraction(fbeam_en, target_name, bind_en, N_tot, fiducialcut);
			rotation->InitSubtraction(StoreEnergy, target_name, bind_en, N_tot, fiducialcut);
			std::cout << " EventLoop: Finished setting up rotation initialize " << std::endl;
		}

		//Resets q vector to (0,0,0)
		rotation->ResetQVector();

		//----------------------------------------//

		// Counters for truth level studies

		int TrueElectronsAboveThreshold = 0;
		int TrueProtonsAboveThreshold = 0;
		int TrueChargedPionsAboveThreshold = 0;
		int TruePiPlusAboveThreshold = 0;
		int TruePiMinusAboveThreshold = 0;
		int TrueGammasAboveThreshold = 0;

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

//			//Fiducial Cuts with the smeared values // moved it further down after W & Q2 cuts
//			if (ApplyFiducials)  { if (!EFiducialCut(fbeam_en,V3_el) ) continue;} // Electron theta & phi fiducial cuts

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

		//----------------------------------------//

		// Explicit cuts on electron momentum

		if (fbeam_en=="1161" && el_momentum < 0.4) { continue; }
		if (fbeam_en=="2261" && el_momentum < 0.55) { continue; }
		if (fbeam_en=="4461" && el_momentum < 1.1) { continue; }

		//Definition as for data. It is also correct for GENIE simulation data since V3_el is rotated above by 180 degree in phi
		double el_phi_mod = V3_el.Phi()*TMath::RadToDeg()  + 30; //Add 30 degree for plotting and photon phi cut
		if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree

		//Calculated Mott Cross Section and Weights for Inclusive Histograms
		//Wght and e_acc_ratio is 1 for CLAS data
		//double Mott_cross_sec = ( pow(fine_struc_const,2.)*(cos(el_theta)+1))/(2*pow(El,2.)*pow((1-cos(el_theta)),2.));

		double reco_Q2 = -(V4_el-V4_beam).Mag2();
		double Q4 = reco_Q2 * reco_Q2;
		double Mott_cross_sec = (1./Q4);

		//----------------------------------------//

		// Sanity check, especially for radiation
		if (wght < 0 || wght > 10) { std::cout << "Something is really wrong with your radiation weights !!!" << std::endl; }

		double WeightIncl = wght*e_acc_ratio / Mott_cross_sec;

		// Securing ourselves against infinities
		if ( fabs(WeightIncl) != WeightIncl ) { continue; }

		//Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
		double nu = -(V4_el-V4_beam).E();
		//double x_bjk = reco_Q2/(2*m_prot*nu);

		// QE selection
		//if ( fabs(x_bjk - 1.) > 0.2) { continue; }

		//----------------------------------------//

		TVector3 V3_q = (V4_beam-V4_el).Vect();
		double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

		//converting theta to degrees
		el_theta = el_theta*TMath::RadToDeg();

		//Cuts on Q2 and W, only keep events with Q2 > Q2cut and W < Wcut
		if ( reco_Q2 < Q2cut || W_var > Wcut) continue;

		//----------------------------------------//

		// apapadop Nov 4 2020: true electron counter for truth level studies
		TrueElectronsAboveThreshold++;

		if (choice > 0) { 			

			//Fiducial Cuts with the smeared values
			if (ApplyFiducials)  { if (!EFiducialCut(fbeam_en,V3_el) ) continue;} // Electron theta & phi fiducial cuts 
		}

		//----------------------------------------//

		//Set q vector for the following rotations for the subtraction procedure
		rotation->SetQVector(V3_q);
//		rotation->PrintQVector();

		int ElectronSector = el_phi_mod / 60.;

		//----------------------------------------//

		// apapadop: Oct 8 2020: ditching bad sectors
		// Counting sectors from 0 to 5

		if (!UseAllSectors) { 

			if ( (ElectronSector == 2 || ElectronSector == 4) &&  fbeam_en == "1161") { continue; }
			if ( (ElectronSector == 2 || ElectronSector == 3 || ElectronSector == 4) &&  fbeam_en == "2261") { continue; }

		}

		//----------------------------------------//
 
		//Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

		//Index variables for hadrons (p and pions)
		int index_p[20]; //index for each proton
		//int index_pi[20]; //index for each pion
		int ind_pi_phot[20]; //index for pions and photons
		//int index_pipl[20]; //index for each pi plus
		//int index_pimi[20]; //index for each pi minus

		int charge_pi[20]; //Charge for the pions and photons
		//Smeared Momentum and Energy values for GENIE (simulation) data
		double Smeared_Pp[20]; //smeared momentum values for protons
		double Smeared_Ep[20]; //smeared energy values for protons
		double Smeared_Ppi[20]; //smeared momentum values for pions
		//double Smeared_Epi[20]; //smeared energy values for pions

		//Number of hadrons
		int num_p = 0;
		int num_pi = 0;
		int num_pi_phot = 0; //couting all pions and photons
		//int num_phot = 0; //couting all photons
		int num_pimi = 0;
		int num_pipl = 0;
		int num_pi_phot_nonrad = 0; //counting all pions and non-radiation photons
		int num_phot_rad = 0; //counting radiation photons
		int num_phot_nonrad = 0;
		//Index and number variables for neutral particles
		//int ec_index_n[20];
		int ec_num_n = 0;
		bool ec_radstat_n[20];

		//Array initialize to -1 or false
		for (int i = 0; i < 20; i++) {
			index_p[i] = -1;   
			//index_pi[i] = -1;   
			//index_pipl[i] = -1;   
			//index_pimi[i] = -1;   
			ind_pi_phot[i] = -1;
			//ec_index_n[i] = -1;   
			ec_radstat_n[i] = false;
			charge_pi[i] = -2; //default number should be not a possible real charge
			Smeared_Pp[i]  = 0; Smeared_Ep[i]  = 0;  //default 0 momentum and energy after smearing
			Smeared_Ppi[i] = 0; 
			//Smeared_Epi[i] = 0;  //default 0 momentum and energy after smearing
		}

		const double phot_rad_cut = 40;
		const double phot_e_phidiffcut = 30; //electron - photon phi difference cut

		// Creating vectors to store id of particles in the array
		vector <int> ProtonID; vector <int> PiPlusID; vector <int> PiMinusID; vector <int> PhotonID;
		ProtonID.clear(); PiPlusID.clear(); PiMinusID.clear();  PhotonID.clear();

		//----------------------------------------//

		//Loop for Hadrons
		for (int i = 0; i < nf; i++) {

			//Start of proton selection

			if (pdgf[i] == 2212  && pf[i] > 0.3) {

				if ( choice > 0 ) { //GENIE data

					//Smearing of proton
					double temp_smear_P = gRandom->Gaus(pf[i],reso_p*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_prot * m_prot );

					TVector3 V3_prot_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_prot = V3_prot_corr.Phi();
					V3_prot_corr.SetPhi(phi_prot + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS

					// apapadop Nov 4 2020: true proton counter for truth level studies above a min theta threshold (12 deg)
					if (PFiducialCutExtra(StoreEnergy, V3_prot_corr)) { TrueProtonsAboveThreshold++; }

					if (ApplyFiducials) { if (!PFiducialCut(StoreEnergy, V3_prot_corr) ) { continue; } } // Proton theta & phi fiducial cuts

					num_p = num_p + 1;
					index_p[num_p - 1] = i;
					ProtonID.push_back(i);
					Smeared_Pp[num_p - 1] = temp_smear_P;
					Smeared_Ep[num_p - 1] = temp_smear_E;

					phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				}
				else { // CLAS data does not need Fiducial Cut again

					num_p = num_p + 1;
					index_p[num_p - 1] = i;
					ProtonID.push_back(i);

					TVector3 V3_prot_corr(pxf[i+60],pyf[i+60],pzf[i+60]);

				}

			}

			//----------------------------------------//

			if (pdgf[i] == -211  && pf[i] > 0.15)  { //Pi minus

				if ( choice > 0) { //GENIE data

					//Smearing of pi minus
					double temp_smear_P = gRandom->Gaus(pf[i],reso_pi*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

					TVector3 V3_pi_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_pion = V3_pi_corr.Phi();
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)

					// apapadop Nov 4 2020: true charged pion counter for truth level studies above a min theta threshold 
					// given by a functional form A + B / P

					if (PimiFiducialCutExtra(StoreEnergy, V3_pi_corr)) { 

						TrueChargedPionsAboveThreshold++;
						TruePiMinusAboveThreshold++;

					}

					if (ApplyFiducials) { if ( !Pi_phot_fid_united(StoreEnergy, V3_pi_corr, -1) ) {  continue; } }

					num_pimi = num_pimi + 1;
					num_pi = num_pi + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					//index_pimi[num_pi_phot - 1] = i;
					//index_pi[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiMinusID.push_back(i);
					charge_pi[num_pi_phot - 1] = -1;
					Smeared_Ppi[num_pi_phot - 1] = temp_smear_P;
					//Smeared_Epi[num_pi_phot - 1] = temp_smear_E;

					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				}
				else { // CLAS data does not need Fiducial Cut again

					num_pimi = num_pimi + 1;
					num_pi = num_pi + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					//index_pimi[num_pi_phot - 1] = i;
					//index_pi[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiMinusID.push_back(i);
					charge_pi[num_pi_phot - 1] = -1;

					TVector3 V3_pi_corr(pxf[i],pyf[i],pzf[i]);

				}

			}

			//----------------------------------------//

			if ( pdgf[i] == 211  && pf[i] > 0.15)  {

				if ( choice > 0) { //GENIE data
					//Smearing of pi plus
					double temp_smear_P = gRandom->Gaus(pf[i],reso_pi*pf[i]);
					//double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

					TVector3 V3_pi_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_pion = V3_pi_corr.Phi();
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)

					// apapadop Nov 4 2020: true charged pion counter for truth level studies with min theta threshold (12 deg)

					if (PiplFiducialCutExtra(StoreEnergy, V3_pi_corr)) { 

						TrueChargedPionsAboveThreshold++;
						TruePiPlusAboveThreshold++;

					}

					if (ApplyFiducials) { if ( !Pi_phot_fid_united(StoreEnergy, V3_pi_corr, 1) )     {  continue; } }

					num_pipl = num_pipl + 1;
					num_pi  = num_pi + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					//index_pipl[num_pi_phot - 1] = i;
					//index_pi[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiPlusID.push_back(i);
					charge_pi[num_pi_phot - 1] = 1;
					Smeared_Ppi[num_pi_phot - 1] = temp_smear_P;
					//Smeared_Epi[num_pi_phot - 1] = temp_smear_E;

					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				}
				else { //CLAS data does not need Fiducial Cut again

					num_pipl = num_pipl + 1;
					num_pi  = num_pi + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					//index_pipl[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiPlusID.push_back(i);
					charge_pi[num_pi_phot - 1] = 1;

					TVector3 V3_pi_corr(pxf[i],pyf[i],pzf[i]);

				}

			}

			//----------------------------------------//

			if (pdgf[i] == 22  && pf[i] > 0.3) {

				//Determine photon vector for the cut on radiation photon via angle with respect to the electron
				TVector3 V3_phot_angles(pxf[i],pyf[i],pzf[i]);
				if (choice > 0) { //GENIE data
					//no smearing of GENIE photons
					double phi_photon = V3_phot_angles.Phi();
					V3_phot_angles.SetPhi(phi_photon + TMath::Pi()); // Vec.Phi() is between (-180,180)

					// apapadop Nov 4 2020: true photon counter for truth level studies

					if (Phot_fidExtra(V3_phot_angles)) { 
	
						TrueGammasAboveThreshold++;

					}

					if (ApplyFiducials) { if ( !Pi_phot_fid_united(fbeam_en, V3_phot_angles, 0) )  { continue;} }

				}

				double neut_phi_mod = V3_phot_angles.Phi()*TMath::RadToDeg() + 30; //Add 30 degree
				if (neut_phi_mod < 0) neut_phi_mod = neut_phi_mod + 360;  //Neutral particle is between 0 and 360 degree

				ec_num_n = ec_num_n + 1;
				num_pi_phot = num_pi_phot + 1;
				ind_pi_phot[num_pi_phot - 1] = i;
				PhotonID.push_back(i);

				Smeared_Ppi[num_pi_phot - 1] = V3_phot_angles.Mag();
				//Smeared_Epi[num_pi_phot - 1] = V3_phot_angles.Mag();

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

		//----------------------------------------//

		// Truth level studies
		// Requiring true level signal 1e 1p 0pi+/- 0 gammas
		// With smearing / fiducial cuts / acceptance maps

		if (TruthLevel1p0piSignalStudy) {

			if (TrueElectronsAboveThreshold != 1) { continue; }
			if (TruthLevel1p0piSignalStudy) { if (TrueProtonsAboveThreshold != 1) { continue; } }
			if (TrueChargedPionsAboveThreshold != 0) { continue; }
			if (TrueGammasAboveThreshold != 0) { continue; }

		}

		//----------------------------------------//

		//Skip event if there is at least one radiation photon

		if (num_phot_rad > 0) {	continue; }

		//----------------------------------------//

		// Events with exactly 2 protons

		if (num_p == 2) {

			//LorentzVectors for protons without momentum smearing or corrections
			TLorentzVector V4_prot_uncorr1(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
			TLorentzVector V4_prot_uncorr2(pxf[index_p[1]],pyf[index_p[1]],pzf[index_p[1]],TMath::Sqrt(m_prot*m_prot+pf[index_p[1]]*pf[index_p[1]]));
			//LorentzVectors for protons with momentum smearing or corrections
			TVector3 V3_prot_corr1;
			TVector3 V3_prot_corr2;

			double p_acc_ratio1 = 1; //will be 1 for CLAS data
			double p_acc_ratio2 = 1; //will be 1 for CLAS data

			if (choice == 0) { //CLAS data
				V3_prot_corr1.SetXYZ(pxf[index_p[0]+60],pyf[index_p[0]+60],pzf[index_p[0]+60]);
				V3_prot_corr2.SetXYZ(pxf[index_p[1]+60],pyf[index_p[1]+60],pzf[index_p[1]+60]);
			}

			if (choice > 0) { //GENIE data, fiducials are done in hadron loop

				V3_prot_corr1.SetXYZ(Smeared_Pp[0]/pf[index_p[0]] * pxf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pyf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pzf[index_p[0]]);
				double phi_prot1 = V3_prot_corr1.Phi();
				V3_prot_corr1.SetPhi(phi_prot1 + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS
				phi_prot1 += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				double p_theta1 = V3_prot_corr1.Theta();
				double prot_mom_corr1 = V3_prot_corr1.Mag();
				//Proton 1 weight
				p_acc_ratio1 = acceptance_c(prot_mom_corr1, cos(p_theta1), phi_prot1, 2212,file_acceptance_p,ApplyAccWeights);
				if ( fabs(p_acc_ratio1) != p_acc_ratio1 ) { continue; }

				V3_prot_corr2.SetXYZ(Smeared_Pp[1]/pf[index_p[1]] * pxf[index_p[1]],Smeared_Pp[1]/pf[index_p[1]] * pyf[index_p[1]],Smeared_Pp[1]/pf[index_p[1]] * pzf[index_p[1]]);
				double phi_prot2 = V3_prot_corr2.Phi();
 				V3_prot_corr2.SetPhi(phi_prot2 + TMath::Pi()); // Vec.Phi() is between (-180,180) // GENIE coordinate system flipped with respect to CLAS
				phi_prot2 += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				double p_theta2 = V3_prot_corr2.Theta();
				double prot_mom_corr2 = V3_prot_corr2.Mag();
				//Proton 2 weight
				p_acc_ratio2 = acceptance_c(prot_mom_corr2, cos(p_theta2), phi_prot2, 2212,file_acceptance_p,ApplyAccWeights);
				if ( fabs(p_acc_ratio2) != p_acc_ratio2 ) { continue; }

			}

			//Total proton weight
			double weight_protons = p_acc_ratio1 * p_acc_ratio2;

			TVector3 V3_2prot_uncorr[2];
			V3_2prot_uncorr[0] = V4_prot_uncorr1.Vect();
			V3_2prot_uncorr[1] = V4_prot_uncorr2.Vect();

			TVector3 V3_2prot_corr[2];
			V3_2prot_corr[0] = V3_prot_corr1;
			V3_2prot_corr[1] = V3_prot_corr2;

			//---------------------------------- 2p 0pi->  1p0pi   ----------------------------------------------

			double E_tot_2p[2]={0};
			double p_perp_tot_2p[2]={0};
			double N_prot_both = 0;
			double P_N_2p[2]={0};

			rotation->prot2_rot_func( V3_2prot_corr, V3_2prot_uncorr, V4_el, E_tot_2p, p_perp_tot_2p, P_N_2p , &N_prot_both);

			if (num_pi_phot == 0 && N_prot_both != 0){

				double histoweight = weight_protons*e_acc_ratio*wght/Mott_cross_sec; //total weight from 2p acceptance , 1e acceptance, Mott, and GENIE weight

				for (int f = 0; f < num_p; f++) { // Looping through two protons

					// -----------------------------------------------------------------------------------------------

					double LocalWeight = -P_N_2p[f]*histoweight;

					// -----------------------------------------------------------------------------------------------

					double STLV[20] = {};
					STV_Tools(V3_el,V3_2prot_corr[f],V4_el.E(),TMath::Sqrt(TMath::Power(V3_2prot_corr[f].Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
					double DeltaAlphaT = STLV[1];
					double Ptx = STLV[13];
					double Pty = STLV[14];	
					double PT = p_perp_tot_2p[f];

					int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
					int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
					int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}							

				} // End of loop over two protons

			}// End of no pion cut and N_prot_both!=0

			//---------------------------------- 2p 1pi   ----------------------------------------------
			//Const int can be placed somewhere up after if for 2 protons F.H. 05.09.19
			const int N_2prot=2;
			//Variable might/could be placed in a more local context F.H. 05.09.19
			double Ecal_2p1pi_to2p0pi[N_2prot]={0};
			double p_miss_perp_2p1pi_to2p0pi[N_2prot]={0};

			if (num_pi_phot == 1) {

				TVector3 V3_1pi_corr;
				double pion_acc_ratio = 1;

				if (choice == 0) { //CLAS data
					V3_1pi_corr.SetXYZ(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]]);
				}

				if (choice > 0) { //GENIE data
					pion_acc_ratio = 0;//reset to 0 just to be save
					V3_1pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);

					double phi_pion = V3_1pi_corr.Phi();
					V3_1pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_1pi_corr.Theta();
					double pion_mom_corr = V3_1pi_corr.Mag();

					if (charge_pi[0] == 1) { //acceptance for pi plus
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip,ApplyAccWeights);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == -1) {    //acceptance for pi minus. using electron acceptance map

						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance_pim,true);

						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == 0) {    //acceptance for neutral, setting to 1 for now F.H. 09/24/19
						pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 2proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue; }
				}

				double P_2p1pito2p0pi[2] = {0};
				double P_2p1pito1p1pi[2] = {0};
				double P_2p1pito1p0pi[2] = {0};
				double Ptot = 0;

				rotation->prot2_pi1_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_1pi_corr, charge_pi[0], V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

				double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec; 
				//Is this correct in the following loop? F.H. 09/01/19

				for (int z=0; z < N_2prot; z++){ //looping over two protons

					//---------------------------------- 2p 1pi ->2p 0pi ----------------------------------------------

					double LocalWeight = P_2p1pito2p0pi[z]*histoweight;

					// -----------------------------------------------------------------------------------------------

					double STLV[20] = {};
					STV_Tools(V3_el,V3_2prot_corr[z],V4_el.E(),TMath::Sqrt(TMath::Power(V3_2prot_corr[z].Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
					double DeltaAlphaT = STLV[1];
					double Ptx = STLV[13];
					double Pty = STLV[14];
					double PT = p_miss_perp_2p1pi_to2p0pi[z];														
			
					int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
					int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
					int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}															

					//---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

					LocalWeight = P_2p1pito1p1pi[z]*histoweight;

					// -----------------------------------------------------------------------------------------------				

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}			

					// -----------------------------------------------------------------------------------------------

					LocalWeight = -P_2p1pito1p0pi[z]*histoweight;

					// -----------------------------------------------------------------------------------------------				

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}

				} // End of filling the histograms for 2protons

			} // 1pi requirement

			//---------------------------------- 2p 2pi ----------------------------------------------

			const int N_2pi=2;
			double Ecal_2p2pi[N_2prot];
			double p_miss_perp_2p2pi[N_2prot];
			double Ptot_2p[2]={0};

			if (num_pi_phot == 2) {

				TVector3 V3_2pi_corr[N_2pi];
				double pion_acc_ratio[N_2pi] = {1};

				for (int i = 0; i < num_pi_phot; i++) {

					if (choice == 0) { //CLAS data
						V3_2pi_corr[i].SetXYZ(pxf[ind_pi_phot[i]],pyf[ind_pi_phot[i]],pzf[ind_pi_phot[i]]);
						pion_acc_ratio[i] = 1; //Acceptance 1 for CLAS data
					}

					if (choice > 0) { //GENIE data
						pion_acc_ratio[i] = 0; //reset to 0 just to be same
						V3_2pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
								Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						double phi_pion = V3_2pi_corr[i].Phi();
						V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_2pi_corr[i].Theta();
						double pion_mom_corr = V3_2pi_corr[i].Mag();

						if (charge_pi[i] == 1) { //acceptance for pi plus
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip,ApplyAccWeights);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == -1) {		//acceptance for pi minus. using electron acceptance map

							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance_pim,true);


							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == 0) {		//acceptance for photon set to 1 for now F.H. 09/24/19
							pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 2proton and 2 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue;
						}
					}

				}

				rotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi_corr,charge_pi, V4_el, Ecal_2p2pi,p_miss_perp_2p2pi,Ptot_2p);

				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
				double histoweight = weight_pions * weight_protons * e_acc_ratio * wght/Mott_cross_sec; 

				for(int z = 0; z < N_2prot; z++){ //looping over two protons

					//---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------

					double LocalWeight = Ptot_2p[z]*histoweight;

					// -----------------------------------------------------------------------------------------------

					double STLV[20] = {};
					STV_Tools(V3_el,V3_2prot_corr[z],V4_el.E(),TMath::Sqrt(TMath::Power(V3_2prot_corr[z].Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
					double DeltaAlphaT = STLV[1];	
					double Ptx = STLV[13];
					double Pty = STLV[14];	
					double PT = p_perp_tot_2p[z];													

					int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
					int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
					int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}

				} //Filling the histogram for two protons

			}//2pi requirement

		} // End of the 2-proton case

		// -------------------------------------------------------------------------------------------------------------------------------------

		// Events with exactly 3 protons

		if(num_p == 3) {

			const int N_3p = 3;
			TLorentzVector V4_p_uncorr[N_3p], V4_p_corr[N_3p],V4_prot_el[N_3p];
			TVector3 V3_prot_uncorr[N_3p],V3_prot_corr[N_3p],V3_3p_rot[N_3p];
			//double E_cal[N_3p],
			double p_miss_perp[N_3p],P_3pto1p[N_3p];
			double N_p1[N_3p]={0};
			double N_p_three=0;
			double E_cal_3pto1p[3]={0};
			double p_miss_perp_3pto1p[3]={0};
			int N_comb = 3;
			const int N_2p = 2;
			double E_cal_3pto2p[3][N_2p]={0};
			double p_miss_perp_3pto2p[3][N_2p]={0};
			double P_3pto2p[3][N_2p]={0};
			TVector3 V3_2p_rot[N_2p], V3_prot_el[N_3p][N_2p];
			double p_acc_ratio[N_3p] = {1};

			for(int i = 0; i < N_3p; i++) {

				N_p1[i] = 0;

				V4_p_uncorr[i].SetPxPyPzE(pxf[index_p[i]],pyf[index_p[i]],pzf[index_p[i]],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]]*pf[index_p[i]]));
				V3_prot_uncorr[i] = V4_p_uncorr[i].Vect();

				if (choice == 0) { //CLAS data

					V3_prot_corr[i].SetXYZ(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60]);
					V4_p_corr[i].SetPxPyPzE(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]+60]*pf[index_p[i]+60]));
					p_acc_ratio[i] = 1; //Acceptance is 1 for CLAS datafile
				}

				if (choice > 0) { //GENIE data

					p_acc_ratio[i] = 0; //Reset just to be sure
					V3_prot_corr[i].SetXYZ(Smeared_Pp[i]/pf[index_p[i]] * pxf[index_p[i]],Smeared_Pp[i]/pf[index_p[i]] * pyf[index_p[i]],
						Smeared_Pp[i]/pf[index_p[i]] * pzf[index_p[i]]);
					V4_p_corr[i].SetPxPyPzE(Smeared_Pp[i]/pf[index_p[i]] * pxf[index_p[i]],Smeared_Pp[i]/pf[index_p[i]] * pyf[index_p[i]],
						Smeared_Pp[i]/pf[index_p[i]] * pzf[index_p[i]],Smeared_Ep[i]);

					double phi_prot = V3_el.Phi(); //in Radians
					V3_prot_corr[i].SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double p_theta = V3_prot_corr[i].Theta();
					double prot_mom_corr = V3_prot_corr[i].Mag();
					//Proton acceptance weight
					p_acc_ratio[i] = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p,ApplyAccWeights);
					if ( fabs(p_acc_ratio[i]) != p_acc_ratio[i] ) { continue; }
				}

				V4_prot_el[i] = V4_p_corr[i] + V4_el;
				p_miss_perp[i] = TMath::Sqrt(V4_prot_el[i].Px()*V4_prot_el[i].Px() + V4_prot_el[i].Py()*V4_prot_el[i].Py());

			} //end loop over N_3p

			V3_prot_el[0][0]=V4_el.Vect()+V3_prot_uncorr[0];
			V3_prot_el[0][1]=V4_el.Vect()+V3_prot_uncorr[1];
			V3_prot_el[1][0]=V4_el.Vect()+V3_prot_uncorr[0];
			V3_prot_el[1][1]=V4_el.Vect()+V3_prot_uncorr[2];
			V3_prot_el[2][0]=V4_el.Vect()+V3_prot_uncorr[1];
			V3_prot_el[2][1]=V4_el.Vect()+V3_prot_uncorr[2];

			rotation->prot3_rot_func( V3_prot_corr,V3_prot_uncorr,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, E_cal_3pto1p,p_miss_perp_3pto1p,&N_p_three);

			//acceptance weight for all three protons ( = 1 for CLAS data)

			double weight_protons =	p_acc_ratio[0] * p_acc_ratio[1] * p_acc_ratio[2];

			if (num_pi_phot == 0 && N_p_three != 0){

				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = weight_protons * e_acc_ratio * wght/Mott_cross_sec; 
				//Weight for 3protons, 1 electron, GENIE weight and Mott cross section

				for(int count = 0; count < N_comb; count++) { // Loop over number of combinations

					for(int j = 0; j < N_2p; j++) { // Loop over two protons

						//-----------------------------------------  3p to 2p->1p  ------------------------------

						double LocalWeight = P_3pto2p[count][j]*histoweight;

						// -----------------------------------------------------------------------------------------------

						double STLV[20] = {};
						STV_Tools(V3_el,V3_prot_corr[j],V4_el.E(),TMath::Sqrt(TMath::Power(V3_prot_corr[j].Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
						double DeltaAlphaT = STLV[1];	
						double Ptx = STLV[13];
						double Pty = STLV[14];
						double PT = p_miss_perp_3pto2p[count][j];	
						
						int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
						int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
						int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
						int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

						if (Interaction > 0) {

							DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
							DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
							DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
							DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

						}						

					} //end loop over protons

				} //end loop over combination N_comb

				//-----------------------------------------  3p to 1p  ------------------------------------------

				for(int j = 0; j < N_3p; j++) {

					double LocalWeight = -P_3pto1p[j]*histoweight;

					// -----------------------------------------------------------------------------------------------

					double STLV[20] = {};
					STV_Tools(V3_el,V3_prot_corr[j],V4_el.E(),TMath::Sqrt(TMath::Power(V3_prot_corr[j].Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
					double DeltaAlphaT = STLV[1];
					double Ptx = STLV[13];
					double Pty = STLV[14];	
					double PT = p_miss_perp[j];													
			
					int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
					int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
					int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}																				

					// -----------------------------------------------------------------------------------------------

				} //end loop over N_3p

			} // end if num_pi_phot==0 && N_p_three!=0, no pions

			//----------------------------------3p 1pi ----------------------------------------------------------

			if (num_pi_phot == 1) {

				double P_tot_3p[N_3p]={0};
				double Ecal_3p1pi[N_3p]={0};
				double p_miss_perp_3p1pi[N_3p]={0};
				TVector3 V3_pi_corr;
				double pion_acc_ratio = 1;

				if (choice == 0) { //CLAS data
					V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]]);
					pion_acc_ratio = 1; //Acceptance is 1 for CLAS datafile
				}

				if (choice > 0){ //GENIE data

					pion_acc_ratio = 0; //Reset to 0 just to be sure
					V3_pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
							Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);

					double phi_pion = V3_pi_corr.Phi(); //in Radians
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_pi_corr.Theta();
					double pion_mom_corr = V3_pi_corr.Mag();

					if (charge_pi[0] == 1) { //acceptance for pi plus
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip,ApplyAccWeights);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == -1) {		//acceptance for pi minus. using electron acceptance map

						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance_pim,true);

						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == 0) {		//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
						 pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 3proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue; }

				}

				rotation->prot3_pi1_rot_func(V3_prot_corr,V3_prot_uncorr, V3_pi_corr, charge_pi[0] , V4_el,	Ecal_3p1pi,p_miss_perp_3p1pi, P_tot_3p);

				//for CLAS data is histoweight = 1/Mott_cross_sec
				double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec; 
				//Weight for 3protons, 1 pion, 1 electron, GENIE weight and Mott cross section

				for(int j = 0; j < N_3p; j++) { //loop over 3 protons

					double LocalWeight = P_tot_3p[j]*histoweight;

					// -----------------------------------------------------------------------------------------------

					double STLV[20] = {};
					STV_Tools(V3_el,V3_prot_corr[j],V4_el.E(),TMath::Sqrt(TMath::Power(V3_prot_corr[j].Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
					double DeltaAlphaT = STLV[1];
					double Ptx = STLV[13];
					double Pty = STLV[14];	
					double PT = p_miss_perp[j];
		
					int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
					int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
					int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}

				} //end loop over N_3p

			} // 1 pi requirement ends

		} //end of 3-proton case

		//------------------------------------------requiring there to be a proton -------------------------------------

		//Events with exactly one proton

		if( num_p == 1) {

			//Vector for proton without momentum smearing
			TLorentzVector V4_prot_uncorr(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
			TVector3 V3_prot_uncorr = V4_prot_uncorr.Vect();

			//Vector for proton with momentum smearing or correction (energy loss)
			TVector3 V3_prot_corr;
			TLorentzVector V4_prot_corr;

			double p_acc_ratio = 1; //acceptance is 1 for CLAS data

			if (choice == 0) { //CLAS data
				V3_prot_corr.SetXYZ(pxf[index_p[0]+60], pyf[index_p[0]+60], pzf[index_p[0]+60]);
				V4_prot_corr.SetPxPyPzE(pxf[index_p[0]+60], pyf[index_p[0]+60], pzf[index_p[0]+60],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]+60]*pf[index_p[0]+60]));
			}
			if (choice > 0) { //GENIE data
				p_acc_ratio = 0; //Reset just to be sure
				//Fiducial cuts are done in the hadron loop
				//Vector for proton with momentum smearing
				V3_prot_corr.SetXYZ(Smeared_Pp[0]/pf[index_p[0]] * pxf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pyf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pzf[index_p[0]]);
				V4_prot_corr.SetPxPyPzE(Smeared_Pp[0]/pf[index_p[0]] * pxf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pyf[index_p[0]],
							Smeared_Pp[0]/pf[index_p[0]] * pzf[index_p[0]],Smeared_Ep[0]);

				double phi_prot = V3_prot_corr.Phi(); //in Radians
				V3_prot_corr.SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
				phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				//Proton kinematic variables
				double p_theta = V3_prot_corr.Theta();
				double prot_mom_corr = V3_prot_corr.Mag();
				//Proton weight
				p_acc_ratio = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p,ApplyAccWeights);
				if ( fabs(p_acc_ratio) != p_acc_ratio ) { continue; }

			}

			TLorentzVector V4_prot_el_tot = V4_prot_corr + V4_el;

			double p_perp_tot = TMath::Sqrt(V4_prot_el_tot.Px()*V4_prot_el_tot.Px() + V4_prot_el_tot.Py()*V4_prot_el_tot.Py());
			//double E_tot = V4_el.E() + V4_prot_corr.E() - m_prot + bind_en[ftarget];

			//These Histograms are events with 1 electron and  1 proton and multiple pions
			//histoweight_inc is 1/Mott_cross_sec for CLAS data
			//double histoweight_inc = p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;
			//double histoweight_NoMott = p_acc_ratio * e_acc_ratio * wght;

			//---------------------------------- 1p 0pi   ----------------------------------------------

			// Main Plots

			if (num_pi_phot == 0){

				SignalEvents++;

				if (Interaction == 1) { QESignalEvents++; }
				else if (Interaction == 2) { MECSignalEvents++; }
				else if (Interaction == 3) { RESSignalEvents++; }
				else if (Interaction == 4) { DISSignalEvents++; }
				else { OtherSignalEvents++; }

				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;

				double LocalWeight = histoweight;

				// -----------------------------------------------------------------------------------------------

				double STLV[20] = {};
				STV_Tools(V3_el,V3_prot_corr,V4_el.E(),TMath::Sqrt(TMath::Power(V3_prot_corr.Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);	
				double DeltaAlphaT = STLV[1];
				double Ptx = STLV[13];
				double Pty = STLV[14];	
				double PT = p_perp_tot;									

				int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
				int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
				int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
				int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

				DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
				DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
				DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
				DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

				if (Interaction > 0) {

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

				}																																								

			} //num pi=0

			//---------------------------------- 1p 1pi   ----------------------------------------------

			if (num_pi_phot == 1) {

				double N_piphot_det;
				double N_piphot_undet;
				TVector3 V3_pi_corr;
				double pion_acc_ratio = 1;

				if (choice == 0) { //CLAS data

					pion_acc_ratio = 1; //Acceptance is 1 for CLAS data
					V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]], pyf[ind_pi_phot[0]], pzf[ind_pi_phot[0]]);
				}

				if (choice > 0) { //GENIE data

					pion_acc_ratio = 1; //Reset to 0 just to be sure
					V3_pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
							Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);
					double phi_pion = V3_pi_corr.Phi(); //in Radians
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_pi_corr.Theta();
					double pion_mom_corr = V3_pi_corr.Mag();

					if (charge_pi[0] == 1) { //acceptance for pi plus
					
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip,ApplyAccWeights);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == -1) { //acceptance for pi minus. using electron acceptance map

						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance_pim,true);

						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == 0) { //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
						pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 1 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

				}

				rotation->prot1_pi1_rot_func(V3_prot_uncorr,V3_pi_corr, charge_pi[0], &N_piphot_det,&N_piphot_undet);

				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = pion_acc_ratio * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; 
				// 1 proton, 1 Pion, 1 electron acceptance, GENIE weight and Mott

				if (N_piphot_det != 0) {

					double LocalWeight = -(N_piphot_undet/N_piphot_det)*histoweight;

					// -----------------------------------------------------------------------------------------------

					double STLV[20] = {};
					STV_Tools(V3_el,V3_prot_corr,V4_el.E(),TMath::Sqrt(TMath::Power(V3_prot_corr.Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
					double DeltaAlphaT = STLV[1];
					double Ptx = STLV[13];
					double Pty = STLV[14];	
					double PT = p_perp_tot;													
		
					int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
					int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
					int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}				

				} //end of N_piphot_det!=0

			 }//end of 1p 1pi requirement

			//---------------------------------- 1p 2pi   ----------------------------------------------

			if(num_pi_phot == 2) {

				const int N_2pi=2;
				TVector3 V3_2pi_corr[N_2pi],V3_2pi_rot[N_2pi],V3_p_rot;
				double P_1p0pi=0;
				double P_1p1pi[N_2pi]={1};

				double pion_acc_ratio[N_2pi] = {0};

				for (int i = 0; i < num_pi_phot; i++) {

					if (choice == 0) { //CLAS data
						V3_2pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
						pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
					}

					if (choice > 0) { //GENIE data

						V3_2pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						double phi_pion = V3_2pi_corr[i].Phi(); //in Radians
						V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_2pi_corr[i].Theta();
						double pion_mom_corr = V3_2pi_corr[i].Mag();

						if (charge_pi[i] == 1) { //acceptance for pi plus
							 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip,ApplyAccWeights);
							  if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == -1) {	//acceptance for pi minus. using electron acceptance map

							//pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance,ApplyAccWeights);

							// UNTIL AXEL CREATES THE CORRECT PIMINUS MAP, WE SET THE PIMINUS ACCEPTANCE TO BE 1
//							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance,false);
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance_pim,true);

							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == 0) {	//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
							 pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 1 Proton 2 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }
					}

				} //end loop over num_pi_phot

				rotation->prot1_pi2_rot_func(V3_prot_uncorr,V3_2pi_corr,charge_pi,&P_1p0pi,P_1p1pi);
				//weight_pions is 1 for CLAS data
				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = weight_pions * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; 
				//1proton, 2 Pion, 1 electron acceptance, GENIE weight and Mott

				//---------------------------------- 1p 2pi->1p1pi   ----------------------------------------------

				for (int z = 0; z < N_2pi; z++) {  //to consider 2 diff. 1pi states

					double LocalWeight = P_1p1pi[z]*histoweight;

					// -----------------------------------------------------------------------------------------------

					double STLV[20] = {};
					STV_Tools(V3_el,V3_prot_corr,V4_el.E(),TMath::Sqrt(TMath::Power(V3_prot_corr.Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);	
					double DeltaAlphaT = STLV[1];
					double Ptx = STLV[13];
					double Pty = STLV[14];	
					double PT = p_perp_tot;		

					int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
					int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
					int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

					if (Interaction > 0) {

						DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
						DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
						DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
						DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

					}																		
				
				} //end loop over N_2pi

				//---------------------------------- 1p 2pi->1p0pi   ----------------------------------------------

				double LocalWeight = -P_1p0pi*histoweight;

				// -----------------------------------------------------------------------------------------------

				double STLV[20] = {};
				STV_Tools(V3_el,V3_prot_corr,V4_el.E(),TMath::Sqrt(TMath::Power(V3_prot_corr.Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
				double DeltaAlphaT = STLV[1];
				double Ptx = STLV[13];
				double Pty = STLV[14];
				double PT = p_perp_tot;												

				int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
				int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
				int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
				int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

				DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
				DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
				DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
				DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

				if (Interaction > 0) {

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

				}		
				
			} //1p 2pi statetment ends

			//---------------------------------- 1p 3pi   ----------------------------------------------

			if(num_pi_phot == 3){

				const int N_3pi=3;
				TVector3 V3_3pi_corr[N_3pi],V3_3pi_rot[N_3pi],V3_p_rot;
				double P_1p3pi = 0;
				double pion_acc_ratio[N_3pi] = {1};

				for (int i = 0; i < num_pi_phot; i++) {

					if (choice == 0) { //CLAS data
						V3_3pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
						pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
					}

					if (choice > 0) { //GENIE data
						pion_acc_ratio[i] = 0; //Reset to 0 just to be sure
						V3_3pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Pp[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Pp[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						double phi_pion = V3_3pi_corr[i].Phi(); //in Radians
						V3_3pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_3pi_corr[i].Theta();
						double pion_mom_corr = V3_3pi_corr[i].Mag();

						if (charge_pi[i] == 1) { //acceptance for pi plus
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip,ApplyAccWeights);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == -1) {	//acceptance for pi minus. using electron acceptance map

							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance_pim,true);

							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == 0) {	//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
							pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 3 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }
					}

				} //end loop over num_pi_phot

				rotation->prot1_pi3_rot_func(V3_prot_uncorr, V3_3pi_corr, charge_pi, &P_1p3pi);
		 		//weight_pions is 1 for CLAS data
				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1] * pion_acc_ratio[2];
				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = weight_pions * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; 
				//1proton, 3 Pion, 1 electron acceptance, GENIE weight and Mott

				//---------------------------------- 1p 3pi->1p 0pi  total ?? F.H. 08/13/19 check logic here compared to 1p 2pi case ----------------------------

				double LocalWeight = P_1p3pi*histoweight;

				// -----------------------------------------------------------------------------------------------

				double STLV[20] = {};
				STV_Tools(V3_el,V3_prot_corr,V4_el.E(),TMath::Sqrt(TMath::Power(V3_prot_corr.Mag(),2.) + TMath::Power(m_prot,2.)),STLV,NProtons[target_name],NNeutrons[target_name],BindE[target_name],EE[target_name],en_beam_Ecal[fbeam_en]);
				double DeltaAlphaT = STLV[1];
				double Ptx = STLV[13];
				double Pty = STLV[14];	
				double PT = p_perp_tot;							

				int DeltaPTTwoDIndex = ReturnIndex(PT, TwoDArrayNBinsDeltaPT);
				int DeltaAlphaTTwoDIndex = ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
				int DeltaPtxTwoDIndex = ReturnIndex(Ptx, TwoDArrayNBinsDeltaPtx);
				int DeltaPtyTwoDIndex = ReturnIndex(Pty, TwoDArrayNBinsDeltaPty);					

				DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][0]->Fill(PT,LocalWeight);
				DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][0]->Fill(DeltaAlphaT,LocalWeight);
				DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][0]->Fill(Ptx,LocalWeight);		
				DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][0]->Fill(Pty,LocalWeight);																														

				if (Interaction > 0) {

					DeltaPT_InDeltaAlphaTPlot[DeltaAlphaTTwoDIndex][Interaction]->Fill(PT,LocalWeight);	
					DeltaAlphaT_InDeltaPTPlot[DeltaPTTwoDIndex][Interaction]->Fill(DeltaAlphaT,LocalWeight);											
					DeltaPtx_InDeltaPtyPlot[DeltaPtyTwoDIndex][Interaction]->Fill(Ptx,LocalWeight);		
					DeltaPty_InDeltaPtxPlot[DeltaPtxTwoDIndex][Interaction]->Fill(Pty,LocalWeight);

				}										

				// -----------------------------------------------------------------------------------------------
		
			}//end of 1p 3pi requirement

		} // End of 1-proton case

	} //end of event loop (jentry)

	gStyle->SetOptFit(1);

	file_out->cd();
	file_out->Write();
	//gDirectory->Write("hist_Files", TObject::kOverwrite);
	// skim_tree->AutoSave();

	// --------------------------------------------------------------------------------------------------------

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "# Processed Events = " << TotalCounter << std::endl;
	std::cout << std::endl << "1e1p0pi Signal # Events = " << SignalEvents << std::endl;
	std::cout << std::endl << "Passing Rate = " << int(double(SignalEvents) / double(TotalCounter)*100.) << " \%"<< std::endl << std::endl;

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
		//double acc_err = (double)sqrt(acc_ratio*(1-acc_ratio)) / (double)num_gen;

		return acc_ratio;

	}

	else { return 1.; }

}
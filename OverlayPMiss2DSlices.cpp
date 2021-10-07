#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------------------

void OverlayPMiss2DSlices(TString Id = "", TString MC = "") {

	// -----------------------------------------------------------------------------

	int font = 132;
	double size = 0.05;

	TString Label = "Data ";
	if (MC != "") { Label = MC + " "; }	

	// -----------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();	

	// -------------------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetTitleSize(0.0,"t");
	gStyle->SetTitleFont(font,"t");

	gStyle->SetTitleOffset(0.9,"xyz");
	gStyle->SetTitleSize(size,"xyz");
	gStyle->SetTitleFont(font,"xyz");
	gStyle->SetLabelSize(size,"xyz");
	gStyle->SetLabelFont(font,"xyz");
	gStyle->SetNdivisions(10,"xyz");

	gStyle->SetOptStat(0);
	gStyle->SetPadRightMargin(0.13);
	gStyle->SetPadLeftMargin(0.13);	
	gStyle->SetPadBottomMargin(0.13);

	gROOT->ForceStyle();

	TGaxis::SetMaxDigits(3);

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Energies

	std::vector<TString> EnergyTString;	std::vector<double> EnergyDouble; std::vector<TString> EnergyDoubleString;
	EnergyTString.push_back("1161"); EnergyDouble.push_back(1.161); EnergyDoubleString.push_back("1.161");
	EnergyTString.push_back("2261"); EnergyDouble.push_back(2.261); EnergyDoubleString.push_back("2.261");	
	EnergyTString.push_back("4461"); EnergyDouble.push_back(4.461); EnergyDoubleString.push_back("4.461");

	const int NEnergy = EnergyTString.size();

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Nuclei

	std::vector<TString> Nucleus; std::vector<TString> NucleusLatex;
	Nucleus.push_back("C12"); NucleusLatex.push_back("^{12}C");
	Nucleus.push_back("4He"); NucleusLatex.push_back("^{4}He");
	Nucleus.push_back("56Fe"); NucleusLatex.push_back("^{56}Fe");

	const int NNuclei = Nucleus.size();

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// PMiss Slices (e.g. 0 < PMiss < 0.1 GeV/c et al)

	std::vector<TString> PMissSlice; std::vector<TString> PMissSliceLabel; 
	PMissSlice.push_back("_0"); PMissSliceLabel.push_back("P_{Miss} < 0.1 GeV/c");
	PMissSlice.push_back("_1"); PMissSliceLabel.push_back("0.1 < P_{Miss} < 0.2 GeV/c");
	PMissSlice.push_back("_2"); PMissSliceLabel.push_back("0.2 < P_{Miss} < 0.3 GeV/c");
	PMissSlice.push_back("_3"); PMissSliceLabel.push_back("0.3 < P_{Miss} < 0.4 GeV/c");			
	PMissSlice.push_back("_4"); PMissSliceLabel.push_back("0.4 < P_{Miss} < 0.5 GeV/c");
	PMissSlice.push_back("_5"); PMissSliceLabel.push_back("0.5 < P_{Miss} < 0.6 GeV/c");
	PMissSlice.push_back("_6"); PMissSliceLabel.push_back("0.6 < P_{Miss} < 0.7 GeV/c");
	PMissSlice.push_back("_7"); PMissSliceLabel.push_back("0.7 < P_{Miss} < 0.8 GeV/c");
	PMissSlice.push_back("_8"); PMissSliceLabel.push_back("0.8 < P_{Miss} < 0.9 GeV/c");
	PMissSlice.push_back("_9"); PMissSliceLabel.push_back("P_{Miss} > 0.9 GeV/c");

	const int NPMissSlices = PMissSlice.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// TFiles to be opened & plots to be used

	TFile* Files[NEnergy][NNuclei];
	TH1D* Plots[NEnergy][NNuclei][NPMissSlices];

	// -----------------------------------------------------------------------------------------------------------------------------------------	

	// Loop over the energies
	for (int WhichEnergy = 0; WhichEnergy < NEnergy; WhichEnergy++) {

		// Loop over the nuclei
		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus++) {

			if (EnergyDouble[WhichEnergy] == 1.161 && !(Nucleus[WhichNucleus] == "C12") ) { continue; }

			// -----------------------------------------------------------------------------------------------

			// Open the relevant file
			if (MC != "") { Files[WhichEnergy][WhichNucleus] = TFile::Open("mySamples/"+Id+"LFNeutrinos_genie_e2a_ep_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+"_neutrino6_united4_radphot_test_SuSav2_Rad.root","readonly"); }
			else { Files[WhichEnergy][WhichNucleus] = TFile::Open("mySamples/"+Id+"LFNeutrinos_data_e2a_ep_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+"_neutrino6_united4_radphot_test.root","readonly"); }

			// -----------------------------------------------------------------------------------------------	

			for (int WhichSlice = 0; WhichSlice < NPMissSlices; WhichSlice++) {

				// Create the relevant canvas & Legend

				TString CanvasName= Id+MC+"PMinus_vs_Pperp_InPMissSlices_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+PMissSlice[WhichSlice];
				TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);	
				can->SetBottomMargin(0.13);

				Plots[WhichEnergy][WhichNucleus][WhichSlice] = (TH1D*)(Files[WhichEnergy][WhichNucleus]->Get("PPerp_PMinus_Slice"+PMissSlice[WhichSlice]));

				Plots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->SetTitle("P_{T} [GeV/c]");
				Plots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->CenterTitle();				

				Plots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetTitle("P^{-} [GeV/c]");				
				Plots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->CenterTitle();

				Plots[WhichEnergy][WhichNucleus][WhichSlice]->GetZaxis()->SetRangeUser(0,1.05*Plots[WhichEnergy][WhichNucleus][WhichSlice]->GetMaximum());

				Plots[WhichEnergy][WhichNucleus][WhichSlice]->SetTitle(Label + " " + NucleusLatex[WhichNucleus] + ", " + EnergyDoubleString[WhichEnergy] + " GeV," + PMissSliceLabel[WhichSlice]);
				Plots[WhichEnergy][WhichNucleus][WhichSlice]->Draw("coltz");

				// --------------------------------------------------------------------------------------------------------

				can->SaveAs("myPlots/"+CanvasName+".pdf");			
				delete can;

			} // End of the PMiss slices

		} // End of the loop over the nuclei

	} // End of the loop over the energies

	// -----------------------------------------------------------------------------------------------------------------------------------------

} // End of the program 

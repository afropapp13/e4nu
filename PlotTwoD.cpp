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

void PlotTwoD(TString Id = "") {

	// -----------------------------------------------------------------------------

	int font = 132;
	double size = 0.05;

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

	// Plots & colors to be used for plotting purposes	

	std::vector<TString> PlotName; std::vector<TString> XaxisLabel;  std::vector<TString> YaxisLabel;
	PlotName.push_back("PPerp_PMinus"); XaxisLabel.push_back("P_{#perp} [GeV/c]");  YaxisLabel.push_back("P^{-} [GeV/c]");
	PlotName.push_back("PLMinusPLFromPMiss_PMiss"); XaxisLabel.push_back("P_{Miss} [GeV/c]"); YaxisLabel.push_back("P_{L} - P_{L,Miss} [GeV/c]");
	PlotName.push_back("kMissMinusPMiss_PMiss"); XaxisLabel.push_back("P_{Miss} [GeV/c]"); YaxisLabel.push_back("k_{Miss} - P_{Miss} [GeV/c]");
	PlotName.push_back("PnProxyMinusPMiss_PMiss"); XaxisLabel.push_back("P_{Miss} [GeV/c]"); YaxisLabel.push_back("P_{n,proxy} - P_{Miss} [GeV/c]");		

	const int NPlots = PlotName.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Slices (all events, ECal reso < 15%, ECal reso > 15%)

	std::vector<TString> Slice; std::vector<TString> SliceLabel; 
	Slice.push_back("_0"); SliceLabel.push_back("");
	Slice.push_back("_1"); SliceLabel.push_back("ECal Reso < 15%");	
	Slice.push_back("_2"); SliceLabel.push_back("ECal Reso > 15%");

	const int NSlices = Slice.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// TFiles to be opened & plots to be used

	TFile* Files[NEnergy][NNuclei];
	TH2D* Plots[NEnergy][NNuclei][NSlices][NPlots];

	// -----------------------------------------------------------------------------------------------------------------------------------------	

	// Loop over the energies
	for (int WhichEnergy = 0; WhichEnergy < NEnergy; WhichEnergy++) {

		// Loop over the nuclei
		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus++) {

			if (EnergyDouble[WhichEnergy] == 1.161 && !(Nucleus[WhichNucleus] == "C12") ) { continue; }

			// -----------------------------------------------------------------------------------------------

			// Open the relevant file
			Files[WhichEnergy][WhichNucleus] = TFile::Open("mySamples/"+Id+"LFNeutrinos_data_e2a_ep_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+"_neutrino6_united4_radphot_test.root","readonly");

			// -----------------------------------------------------------------------------------------------

			for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

				// -----------------------------------------------------------------------------------------------	

				// Loop over the slices
				for (int WhichSlice = 0; WhichSlice < NSlices; WhichSlice++) {

					// --------------------------------------------------------------------------------------------------------		

					TString CanvasName2D = Id+PlotName[WhichPlot]+"_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice];
					TCanvas* can2D = new TCanvas(CanvasName2D,CanvasName2D,205,34,768,768);	
					can2D->SetBottomMargin(0.12);
					can2D->SetLeftMargin(0.12);
					can2D->SetRightMargin(0.14);				

					// 2D Pmiss vs kMiss plots

					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot] = (TH2D*)(Files[WhichEnergy][WhichNucleus]->Get(PlotName[WhichPlot]+Slice[WhichSlice]));

					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetXaxis()->CenterTitle();
					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetXaxis()->SetTitle(XaxisLabel[WhichPlot]);
					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetXaxis()->SetNdivisions(8);				

					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetYaxis()->CenterTitle();
					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetYaxis()->SetTitle(YaxisLabel[WhichPlot]);
					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetYaxis()->SetTitleOffset(1.1);				
					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetYaxis()->SetNdivisions(8);	

					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetZaxis()->SetRangeUser(0.,1.1*Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->GetMaximum());

					Plots[WhichEnergy][WhichNucleus][WhichSlice][WhichPlot]->Draw("coltz");

					can2D->SaveAs("myPlots/"+CanvasName2D+".pdf");			
					delete can2D;				
				
				} // End of the loop over the slices 

			} // End of the loop over the 2D plots	

		} // End of the loop over the nuclei

	} // End of the loop over the energies

	// -----------------------------------------------------------------------------------------------------------------------------------------

} // End of the program 

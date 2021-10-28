#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TArrayD.h>

#include <iostream>
#include <vector>

using namespace std;

#include "../myFunctions.cpp"
#include "../AfroConstants.h"

// ------------------------------------------------------------------------ //

TGraph* GraghScale(TGraph* g) {

	const int nPoints = g->GetN();

	for (int i = 0; i < nPoints; i++) {

		g->SetPointY(i,0.001 * g->GetY()[i]);

	}

	return g;

}

// ------------------------------------------------------------------------ //

void OverlayXSecs() {

	// ------------------------------------------------------------------------ //

	GlobalSettings();
	//TGaxis::SetMaxDigits(1);

	// ------------------------------------------------------------------------ //

	std::vector<TString> nucleus; 
	std::vector<TString> JustNucleus;
	std::vector<TString> E;
	std::vector<TString> LabelE;	
	std::vector<TString> Theta;	

	// ------------------------------------------------------------------------ //

	//nucleus.push_back("4He"); JustNucleus.push_back("He");
	nucleus.push_back("12C"); JustNucleus.push_back("^{12}C");
	nucleus.push_back("56Fe"); JustNucleus.push_back("^{56}Fe");

	int NNuclei = nucleus.size();			

	// ------------------------------------------------------------------------ //

	//E.push_back("1_161"); LabelE.push_back("1.161"); Theta.push_back("37.5");
	E.push_back("2_261"); LabelE.push_back("2.261"); Theta.push_back("27");	
	E.push_back("4_461"); LabelE.push_back("4.461"); Theta.push_back("21");	

	int NEnergies = E.size();	

	// ------------------------------------------------------------------------ //

	TString NameOfPlots = "h1_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_0";	
	TString OutputPlotNames = "Omega_0_AllEvents";

	// ------------------------------------------------------------------------ //

	std::vector<TString> SimLabels = {"","QE","MEC","RES","DIS"};
	const int NLabels = SimLabels.size();

	// ------------------------------------------------------------------------ //

	std::vector<TGraph*> SuSav2SimPlots; SuSav2SimPlots.resize(NLabels);
	std::vector<TGraph*> G2018SimPlots; G2018SimPlots.resize(NLabels);		

	// ------------------------------------------------------------------------ //

	// Loop over the energies

	for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

		// ------------------------------------------------------------------------ //

		// Loop over the nuclei

		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

			// ------------------------------------------------------------------------ //

			if ( nucleus[WhichNucleus] != "12C" && E[WhichEnergy] == "1_161" ) { continue; }

			// ------------------------------------------------------------------------ //

			// XSec files

			TFile* fData = TFile::Open("myXSec/XSec_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_GeV_NoxBCut.root","readonly");
			TFile* fSuSav2 = TFile::Open("myXSec/"+nucleus[WhichNucleus]+"_E_"+E[WhichEnergy]+"_SuSav2.root","readonly");	
			TFile* fG2018 = TFile::Open("myXSec/"+nucleus[WhichNucleus]+"_E_"+E[WhichEnergy]+"_G2018.root","readonly");	

			// ------------------------------------------------------------------------ //

			// XSec plots

			TH1D* hData = (TH1D*)(fData->Get("Data_"+NameOfPlots));

			for (int WhichSim = 0; WhichSim < NLabels; WhichSim++) {

				SuSav2SimPlots[WhichSim] = GraghScale( (TGraph*)(fSuSav2->Get("Simulation"+SimLabels[WhichSim])) );
				G2018SimPlots[WhichSim] = (TGraph*)(fG2018->Get("Simulation"+SimLabels[WhichSim]));				

			}

			double max = TMath::Max(hData->GetMaximum(),SuSav2SimPlots[0]->GetMaximum());													
//			double max = hData->GetMaximum();

			// ------------------------------------------------------------------------ //							

			TString CanvasName = nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots+"_NoxBCut";
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->SetBottomMargin(0.15);
			PlotCanvas->SetLeftMargin(0.18);			

			// ------------------------------------------------------------------------ //

			// Data

			hData->GetXaxis()->SetNdivisions(8);
			hData->GetXaxis()->CenterTitle();
			hData->GetXaxis()->SetTitle("Energy Transfer [GeV]");
			hData->GetXaxis()->SetLabelOffset(0.01);			

			hData->GetYaxis()->SetRangeUser(0.,1.1*max);
			TString DoubleXSecTitle = "#frac{d^{2}#sigma}{d#Omega dE} [#frac{#mub}{sr GeV "+ JustNucleus[WhichNucleus] +"}]";	
			hData->GetYaxis()->SetTitle(DoubleXSecTitle);					
			hData->GetYaxis()->SetLabelOffset(0.01);			

			hData->SetTitle(LabelE[WhichEnergy] + " GeV & #theta = " + Theta[WhichEnergy] + "^{o}");

			hData->Draw();

			for (int WhichSim = 0; WhichSim < NLabels; WhichSim++) {

				SuSav2SimPlots[WhichSim]->Draw("c same");

			}				

			// ------------------------------------------------------------------------ //

			// SuSav2

			//hSuSav2->Draw("c hist same");
			//hSuSav2QE->Draw("c hist same");

			// ------------------------------------------------------------------------ //

			// G2018

			//hG2018->SetLineStyle(kDashed);
			//hG2018->Draw("c hist same");			

			// ------------------------------------------------------------------------ //

			TString CanvasSaveName = "myPlots/"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_" +OutputPlotNames+"_OverlayXSec";
			PlotCanvas->SaveAs(CanvasSaveName+".pdf");
			PlotCanvas->SaveAs(CanvasSaveName+".eps");

			//delete PlotCanvas;
	
			// ------------------------------------------------------------------------ //			

		} // End of the loop over the nuclei

	} // End of the loop over the energies

} // End of the program
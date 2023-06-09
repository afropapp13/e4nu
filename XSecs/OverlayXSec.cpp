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

#include <iostream>
#include <vector>

using namespace std;

#include "../myFunctions.cpp"
#include "../AfroConstants.h"

// ----------------------------------------------------------------------------------------------------------------

void OverlayXSec() {

	// ------------------------------------------------------------------------

	GlobalSettings();
	//TGaxis::SetMaxDigits(3);

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut; 
	std::vector<TString> nucleus; 
	std::vector<TString> JustNucleus;
	std::vector<TString> E;
	std::vector<TString> LabelE;	
	std::vector<double>  DoubleE;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots;
	std::vector<TString> NameOfSubPlots;
	std::vector<TString> LabelOfSubPlots;
	std::vector<TString> XaxisLabel;			 
	std::vector<TString> LabelOfPlots;  	
	std::vector<TString> Yaxis;

	std::vector<int> Colors = {kBlue+2,kOrange+7,kGreen+2};	

	// ------------------------------------------------------------------------

	nucleus.push_back("4He"); JustNucleus.push_back("He");
	nucleus.push_back("12C"); JustNucleus.push_back("C");
	nucleus.push_back("56Fe"); JustNucleus.push_back("Fe");		

	// ------------------------------------------------------------------------

	E.push_back("1_161"); DoubleE.push_back(1.161); LabelE.push_back("1.159");
	E.push_back("2_261"); DoubleE.push_back(2.261); LabelE.push_back("2.257");	
	E.push_back("4_461"); DoubleE.push_back(4.461); LabelE.push_back("4.453");	

	// ------------------------------------------------------------------------

	xBCut.push_back("NoxBCut");	

	// ------------------------------------------------------------------------	

	//NameOfPlots = {"Total","Longitudinal"};
	NameOfPlots = {"Total"};			

	// ------------------------------------------------------------------------

	std::vector<TH1D*> Plots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NPlots = NameOfPlots.size();	

	TFile* File[NxBCuts][NEnergies][NNuclei];

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				// ---------------------------------------------------------------------------------------------------------------------------------------------

				if ( nucleus[WhichNucleus] != "12C" && E[WhichEnergy] == "1_161" ) { continue; }

				// ---------------------------------------------------------------------------------------------------------------------------------------------

				// File to store the extracted data/mc xsecs	

				File[WhichxBCut][WhichEnergy][WhichNucleus] = TFile::Open("myXSec/XSec_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_GeV_"+xBCut[WhichxBCut]+".root","readonly");	

				// ---------------------------------------------------------------------------------------------------------------------------------------------											

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					NameOfSubPlots.clear();
					LabelOfSubPlots.clear();					
					if (NameOfPlots[WhichPlot] == "Total") {  

						NameOfSubPlots.push_back("PMiss_0"); LabelOfSubPlots.push_back("P_{Miss}");  XaxisLabel.push_back("P_{Miss}");
						//NameOfSubPlots.push_back("kMiss_0"); LabelOfSubPlots.push_back("k_{Miss}");  XaxisLabel.push_back("P_{Miss}");
						NameOfSubPlots.push_back("PnProxy_0"); LabelOfSubPlots.push_back("P_{n,proxy}");  XaxisLabel.push_back("P_{Miss}");

					} else if (NameOfPlots[WhichPlot] == "Longitudinal") {

						NameOfSubPlots.push_back("PLFromPMiss_0"); LabelOfSubPlots.push_back("P_{L,Miss}"); XaxisLabel.push_back("P_{L,Miss}");
						NameOfSubPlots.push_back("PL_0"); LabelOfSubPlots.push_back("P_{L}"); XaxisLabel.push_back("P_{L,Miss}");

					} else { cout << "No clue what you are asking me to plot!" << endl; }

					int NSubPlots = NameOfSubPlots.size();

					std::vector<TH1D*> Plots; Plots.resize(NSubPlots);

					// ---------------------------------------------------------------------------------------------------------------------------------------------				

					TString CanvasName = nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+xBCut[WhichxBCut];
					TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
					PlotCanvas->SetBottomMargin(0.18);
					PlotCanvas->SetLeftMargin(0.2);					

					TLegend* leg = new TLegend(0.71,0.5,0.86,0.7);	
					if (NameOfPlots[WhichPlot] == "Longitudinal") { leg = new TLegend(0.21,0.45,0.36,0.7); }	

					leg->SetNColumns(1);	

					double max = -99.;
					double min = 1E12;					

					// ---------------------------------------------------------------------------------------------------------------------------------------------																

					// Loop over the plots				

					for (int WhichSubPlot = 0; WhichSubPlot < NSubPlots; WhichSubPlot ++) {

						// ---------------------------------------------------------------------------------------

						Plots[WhichSubPlot] = (TH1D*)File[WhichxBCut][WhichEnergy][WhichNucleus]->Get("Data_"+NameOfSubPlots[WhichSubPlot]);

						Plots[WhichSubPlot]->SetMarkerColor(Colors[WhichSubPlot]);
						Plots[WhichSubPlot]->SetLineColor(Colors[WhichSubPlot]);

						double max = TMath::Max(Plots[WhichSubPlot]->GetMaximum(),Plots[0]->GetMaximum());						

						Plots[WhichSubPlot]->Draw("e1x0 same");
						Plots[0]->GetYaxis()->SetRangeUser(0.,1.03*max);
						Plots[0]->GetYaxis()->SetTitleOffset(1.3);

						Plots[0]->GetXaxis()->SetTitle("(e,e'p)_{1p0#pi} " + XaxisLabel[WhichPlot] + " [GeV/c]");						

						Plots[0]->Draw("e1x0 same");						

						TString CopyString = NameOfSubPlots[WhichSubPlot];
						leg->AddEntry(Plots[WhichSubPlot],LabelOfSubPlots[WhichSubPlot],"p");

						// --------------------------------------------------------------------------------------

					} // End of the loop over the sub plots

					gPad->RedrawAxis();

					leg->SetBorderSize(0);
					leg->SetTextFont(FontStyle);
					leg->SetTextSize(0.06);					
					leg->Draw();					

					// -------------------------------------------------------------------------------------------
	
					TString ext = "";
					if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

					TString CanvasSaveName = "myPlots/"+ext+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_OverlayDataXSec_"+NameOfPlots[WhichPlot];
					PlotCanvas->SaveAs(CanvasSaveName+".pdf");
					//PlotCanvas->SaveAs(CanvasSaveName+".eps");

					delete PlotCanvas;				
	
				} // End of the loop over the plots

				File[WhichxBCut][WhichEnergy][WhichNucleus]->Close();
				
				// --------------------------------------------------------------------------------------				

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program

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

#include <iostream>
#include <vector>

using namespace std;

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"

#include "AfroConstants.h"
#include "myFunctions.cpp"

// ----------------------------------------------------------------------------------------------------------------

void OverlayPlots() {

	// ------------------------------------------------------------------------

	SetOffsetAndSize();
	TGaxis::SetMaxDigits(3);
//	TGaxis::SetExponentOffset(-0.1, 1., "y");

	int Ndivisions = 8;
	int LineWidth = 3;
	int FontStyle = 132;
	double TextSize = 0.08;
	
	// ------------------------------------------------------------------------

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> JustNucleus; std::vector<TString> LabelsOfSamples; 
	std::vector<TString> E; std::vector<double> DoubleE;
	std::vector<TString> LabelE; std::vector<TString> FSIModel; std::vector<TString> DirNames;  std::vector<int> BreakDownColors;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots; std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames; std::vector<TH1D*> BreakDownPlots;
	std::vector<int> Colors;
	std::vector<int> Style;

//	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He"); JustNucleus.push_back("He");
	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C"); JustNucleus.push_back("C");
//	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");  JustNucleus.push_back("Fe");
//	nucleus.push_back("CH2"); LabelsOfSamples.push_back("CH2"); JustNucleus.push_back("CH2");

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV"); DoubleE.push_back(1.161);
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV"); DoubleE.push_back(2.261);	
//	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");  DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");
 
//	Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kBlue); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);
	Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kBlue); Colors.push_back(kGreen+1); Colors.push_back(kOrange+1); Colors.push_back(kMagenta);

	Style.push_back(1); Style.push_back(1); Style.push_back(1); Style.push_back(1);

	BreakDownColors.push_back(kBlue); BreakDownColors.push_back(kCyan); BreakDownColors.push_back(kGreen+1); BreakDownColors.push_back(kMagenta);

	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
	FSIModel.push_back("SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger"); FSILabel.push_back("SuSav2 (Rad)");  DirNames.push_back("SuSav2_RadCorr");
	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc"); FSILabel.push_back("SuSav2 (NoRad)");  DirNames.push_back("SuSav2_NoRadCorr");	
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc"); FSILabel.push_back("G2018 (NoRad)");  DirNames.push_back("G2018_NoRadCorr");
	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("SuSav2 (NoRad Truth)");  DirNames.push_back("SuSav2_Truth");	
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("G2018 (NoRad Truth)");  DirNames.push_back("G2018_Truth");

// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	NameOfPlots.push_back("h1_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_0"); LabelOfPlots.push_back("Energy Transfer [GeV]"); OutputPlotNames.push_back("Omega_0");

// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> Plots_Clones;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	TString WhatModelsAreIncluded = "";
	for (int LoopOverFSIModels = 0 ; LoopOverFSIModels < NFSIModels ; LoopOverFSIModels ++) { WhatModelsAreIncluded += "_"+DirNames[LoopOverFSIModels]; };

	// --------------------------------------------------------------------

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					TString PlotCanvasName = nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+OutputPlotNames[WhichPlot]+"_"+xBCut[WhichxBCut];
					TCanvas* PlotCanvas = new TCanvas(PlotCanvasName, PlotCanvasName, 205,34,1024,768);

					// ----------------------------------------------------------------------------------------

					Plots.clear();

					TLegend* leg = new TLegend(0.15,0.57,0.45,0.87);
					leg->SetNColumns(1);

					double max = -99.;
					double min = 1E12;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
						CenterAxisTitle(Plots[WhichFSIModel]);

						// --------------------------------------------------------------------------------------

						// X-axis label

						Plots[WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleOffset(1.05);

						Plots[WhichFSIModel]->GetXaxis()->SetTitle(JustNucleus[WhichNucleus]+LabelOfPlots[WhichPlot]);

						// --------------------------------------------------------------------------------------

						Plots[WhichFSIModel]->GetYaxis()->SetTitle("Normalized Yield");
						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize);

						// --------------------------------------------------------------------------------------

						Plots[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.6); 
						Plots[WhichFSIModel]->SetLineWidth(LineWidth);

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) 
							{ leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "lep");}
						else { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l"); }

						// --------------------------------------------------------------------------------------

						// Scaling Factor

						UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

						double localmax = Plots[WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }
						double height = 1.05;
						Plots[0]->GetYaxis()->SetRangeUser(0.,height*max);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

						Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions);
						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);

						// ----------------------------------------------------------------------------------------------------

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(2.); 
							Plots[WhichFSIModel]->SetMarkerColor(Colors[WhichFSIModel]); 

							gStyle->SetErrorX(0); // Removing the horizontal errors
							Plots[WhichFSIModel]->Draw("e same"); 

						} else { 

							Plots[WhichFSIModel]->SetLineStyle(Style[WhichFSIModel]); 
							Plots[WhichFSIModel]->Draw("hist same"); // draw them as histos
							Plots[0]->Draw("e same"); 
							
						}

						// ----------------------------------------------------------------------------

					} // End of the loop over the FSI Models 

					leg->SetBorderSize(0);
					leg->SetTextFont(FontStyle);
					leg->SetTextSize(TextSize);
					leg->Draw(); // Just data + e.g. susav2

					//delete PlotCanvas;

					// -----------------------------------------------------------------------------------------------------------------------------------------

				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program

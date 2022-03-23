#include <TFile.h>
#include <TH1D.h>
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

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"

#include "AfroConstants.h"
#include "myFunctions.cpp"

using namespace std;

//------------------------------------//

void OverlayPlots() {

	//------------------------------------//

	SetOffsetAndSize();
	TGaxis::SetMaxDigits(3);

	int Ndivisions = 8;
	int LineWidth = 3;
	int FontStyle = 132;
	double TextSize = 0.06;
	TString xBCut = "NoxBCut";	
	
	//------------------------------------//

	std::vector<TString> nucleus; 
	std::vector<TString> JustNucleus; 
	std::vector<TString> E; 
	std::vector<double> DoubleE;
	std::vector<TString> LabelE; 
	std::vector<TString> FSIModel; 
	std::vector<TString> DirNames;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots;
	std::vector<int> Colors;
	std::vector<int> Style;	

	//------------------------------------//	

//	nucleus.push_back("4He");  JustNucleus.push_back("He");
//	nucleus.push_back("12C"); JustNucleus.push_back("C");
	nucleus.push_back("56Fe"); JustNucleus.push_back("Fe");

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV"); DoubleE.push_back(1.161);
//	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV"); DoubleE.push_back(2.261);	
	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");  DoubleE.push_back(4.461);
 
	Colors.push_back(kBlack); Colors.push_back(kBlue-7); Colors.push_back(kOrange+7); Colors.push_back(kGreen+1); Colors.push_back(kOrange+1); Colors.push_back(kMagenta);

//	Style.push_back(9); Style.push_back(9); Style.push_back(9); Style.push_back(9); // fancy dashed lines 
	Style.push_back(1); Style.push_back(1); Style.push_back(1); Style.push_back(1);

	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
//	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc"); FSILabel.push_back("SuSav2 (No Rad)");  DirNames.push_back("hA2018_Truth_RadCorr");
	FSIModel.push_back("SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger"); FSILabel.push_back("SuSav2 (Rad)");  DirNames.push_back("SuSav2_Truth_RadCorr");

	//------------------------------------//

	NameOfPlots = PlotNames;

	//NameOfPlots.push_back("DeltaPT_InDeltaAlphaT_Slice_0Plot");
	//NameOfPlots.push_back("DeltaPT_InDeltaAlphaT_Slice_1Plot");
	//NameOfPlots.push_back("DeltaPT_InDeltaAlphaT_Slice_2Plot");
	//NameOfPlots.push_back("DeltaPT_InDeltaAlphaT_Slice_3Plot");	

	//NameOfPlots.push_back("DeltaAlphaT_InDeltaPT_Slice_0Plot");
	//NameOfPlots.push_back("DeltaAlphaT_InDeltaPT_Slice_1Plot");
	//NameOfPlots.push_back("DeltaAlphaT_InDeltaPT_Slice_2Plot");	

	//NameOfPlots.push_back("DeltaPtx_InDeltaPty_Slice_0Plot");
	//NameOfPlots.push_back("DeltaPtx_InDeltaPty_Slice_1Plot");
	//NameOfPlots.push_back("DeltaPtx_InDeltaPty_Slice_2Plot");				

	//NameOfPlots.push_back("DeltaPty_InDeltaPtx_Slice_0Plot");
	//NameOfPlots.push_back("DeltaPty_InDeltaPtx_Slice_1Plot");
	//NameOfPlots.push_back("DeltaPty_InDeltaPtx_Slice_2Plot");

	//------------------------------------//	

	int NEnergies = E.size();
	int NNuclei = nucleus.size();
	int NPlots = NameOfPlots.size();
	int NFSIModels = FSIModel.size();

	TH1D* Plots[NEnergies][NNuclei][NPlots][NFSIModels];	

	//------------------------------------//	

	// Loop over the energies

	for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

		// Loop over the nuclei

		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

				//------------------------------------//				

				TString CanvasName = nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut;
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				PlotCanvas->SetBottomMargin(0.15);

				TLegend* leg = new TLegend(0.6,0.7,0.75,0.87);
				leg->SetNColumns(1);
				leg->SetBorderSize(0);
				leg->SetTextFont(FontStyle);
				leg->SetTextSize(TextSize);				

				//------------------------------------//					

				// Loop over the FSI Models

				for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

					TString PathToFiles = "myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut+"/";
					TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
					TFile* FileSample = TFile::Open(FileName);

					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel] = (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) );					
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->SetLineWidth(LineWidth);	
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->SetTitle(nucleus[WhichNucleus] + LabelE[WhichEnergy] + ", " + LatexLabel[ NameOfPlots[WhichPlot] ]);
					gStyle->SetTitleFont(FontStyle,"t");
					gStyle->SetTitleFontSize(TextSize);																			

					//------------------------------------//

					// X-axis

					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetXaxis()->CenterTitle();
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetXaxis()->SetLabelSize(TextSize);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetXaxis()->SetTitleSize(TextSize);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetXaxis()->SetTitleOffset(1.05);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions);

					// Y-axis
	
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetYaxis()->SetTitle("Normalized Yield");
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetYaxis()->CenterTitle();						
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetYaxis()->SetTitleOffset(0.8);
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);					

					if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos ) { 

						//leg->AddEntry(Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel],LatexLabel[ NameOfPlots[WhichPlot] ], "");						
						leg->AddEntry(Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel],FSILabel[WhichFSIModel], "lep");
						
					}
					else { leg->AddEntry(Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel],FSILabel[WhichFSIModel], "l"); }

					//------------------------------------//

					// Normalized yield

					UniversalE4vFunction(Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

					//------------------------------------//

					if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

						Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->SetMarkerStyle(20); 
						Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->SetMarkerSize(2.); 
						Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->SetMarkerColor(Colors[WhichFSIModel]); 

						gStyle->SetErrorX(0); // Removing the horizontal errors
						Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->Draw("e same"); 
					
					} else { 
					
						Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->SetLineStyle(Style[WhichFSIModel]); 
						Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->Draw("hist same"); // draw them as histos

						double max = TMath::Max(Plots[WhichEnergy][WhichNucleus][WhichPlot][0]->GetMaximum(),Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichFSIModel]->GetMaximum());
						Plots[WhichEnergy][WhichNucleus][WhichPlot][0]->GetYaxis()->SetRangeUser(0.,1.1*max);				
						Plots[WhichEnergy][WhichNucleus][WhichPlot][0]->Draw("e same"); 
					
					}

				} // End of the loop over the FSI Models 

				leg->Draw();				

				//delete PlotCanvas;

				//------------------------------------//

			} // End of the loop over the plots

		} // End of the loop over the nuclei

	} // End of the loop over the energies

} // End of the program
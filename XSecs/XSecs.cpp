#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TGaxis.h>
#include <TPad.h>

#include <iostream>
#include <vector>

using namespace std;

#include "../myFunctions.cpp"
//#include "../AfroConstants.h"
#include "../Constants.h"

//------------------------------//

void XSecs() {

	//------------------------------//

	GlobalSettings();

	TString xBCut = "NoxBCut";	

	//------------------------------//

	std::vector<TString> nucleus; 
	std::vector<TString> JustNucleus;
	std::vector<TString> E;
	std::vector<TString> LabelE;	
	std::vector<double> DoubleE;
	std::vector<TString> FSIModel;
	std::vector<TString> FSILabel; 
	std::vector<TString> LabelOfPlots;  	
	std::vector<TString> BreakDown;	
	std::vector<TString> OutputPlotNames;

	//------------------------------//

	nucleus.push_back("4He"); JustNucleus.push_back("He");
	nucleus.push_back("12C"); JustNucleus.push_back("C");
	nucleus.push_back("56Fe"); JustNucleus.push_back("Fe");		

	//------------------------------//

	E.push_back("1_161"); DoubleE.push_back(1.161); LabelE.push_back("1.161");
	E.push_back("2_261"); DoubleE.push_back(2.261); LabelE.push_back("2.261");	
	E.push_back("4_461"); DoubleE.push_back(4.461); LabelE.push_back("4.461");		

	//------------------------------//

	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Pinned Data");

	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("SuSav2");	
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("G2018");

	//------------------------------//

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> BreakDownPlots;

	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = PlotNames.size();

	//------------------------------//

	// Loop over the energies

	for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

		// Loop over the nuclei

		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

			//------------------------------//

			if ( nucleus[WhichNucleus] != "12C" && E[WhichEnergy] == "1_161" ) { continue; }

			//------------------------------//

			// File to store the extracted data/mc xsecs	

			TFile* f = TFile::Open("myXSec/XSec_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_GeV_"+xBCut+".root","recreate");	

			//------------------------------//

			// Loop over the plots				

			for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

				//------------------------------//

				TString CanvasName = nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+PlotNames[WhichPlot]+"_"+xBCut;
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				PlotCanvas->SetTopMargin(0.1);
				PlotCanvas->SetBottomMargin(0.18);
				PlotCanvas->SetLeftMargin(0.18); 
				PlotCanvas->SetRightMargin(0.04);				

				//------------------------------//

				Plots.clear();

				TLegend* legGenie = new TLegend(0.71,0.45,0.86,0.7);
				TLegend* legGenieBlackLine = new TLegend(0.71,0.74,0.86,0.88);
				TLegend* legGenieBreak = new TLegend(0.74,0.63,0.91,0.75);					
				TLegend* legG2018 = new TLegend(0.71,0.57,0.89,0.63);					

				legGenie->SetNColumns(1);
				legGenieBlackLine->SetNColumns(1);
				legGenieBreak->SetNColumns(2);
				legGenieBreak->SetMargin(0.35);	
				legG2018->SetNColumns(2);
				legG2018->SetMargin(0.39);									

				double max = -99.;
				double min = 1E12;

				// Loop over the FSI Models

				TH1D* DataPlot = nullptr;

				for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

					TString PathToFiles = "../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut+"/";
					TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
					TFile* FileSample = TFile::Open(FileName,"readonly");

					TH1D* ModelPlot = (TH1D*)FileSample->Get(PlotNames[WhichPlot]);
					ModelPlot->SetDirectory(0); // to decouple it from the open file directory
					Plots.push_back(ModelPlot);

					// --------------------------------------------------------------------------------------

					// Make the plot pretty

					Plots[WhichFSIModel]->SetLineColor(DataSetColors[WhichFSIModel]);
					PrettyDoubleXSecPlot(Plots[WhichFSIModel]);

					Plots[WhichFSIModel]->GetXaxis()->CenterTitle(0);
					Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(1.1);

					// -----------------------------------------------------------------------------------

					// Use the universal e4v function that 
					// scales by all the necessary scaling factors
					// divides by the bin width to obtain normalized yields / absolute cross sections
					// uses the relevant binning
					// gets the relevant x axis range
					// If data sample, also: 
					// 	apply systematics due to rotations et al
					// 	apply acceptance systematics using sector-by -sector uncertainties
					//	apply acceptance correction uncertainties	

					UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy],PlotNames[WhichPlot]);

					// ----------------------------------------------------------------------------------

					// Genie Break Down

					if (
						FSILabel[WhichFSIModel] == "SuSav2"
					) {

						legGenie->AddEntry(Plots[0],"Data", "lep"); 
						legGenieBlackLine->AddEntry(Plots[0],"Data", "lep"); 

						legGenie->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l"); 

						BreakDownPlots.clear();

						for (int j = 1; j < 5; j++) {

							TH1D* BreakDownPlot = (TH1D*)(FileSample->Get(InteractionLabels[j]+PlotNames[WhichPlot]));
							BreakDownPlot->SetDirectory(0); // to decouple it from the open file directory
							BreakDownPlots.push_back(BreakDownPlot);

							UniversalE4vFunction(BreakDownPlots[j-1],FSIModelsToLabels[FSIModel[WhichFSIModel]],\
									     		 nucleus[WhichNucleus],E[WhichEnergy],PlotNames[WhichPlot]);

							//-----------------------------------------------------------------------------------------------

							BreakDownPlots[j-1]->SetLineColor(BreakDownColors[j]);
							BreakDownPlots[j-1]->SetLineWidth(3);

							TLegendEntry* l1 = legGenie->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
							l1->SetTextColor(BreakDownColors[j]);

							TLegendEntry* l1Break = legGenieBreak->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
							l1Break->SetTextColor(BreakDownColors[j]);

							BreakDownPlots[j-1]->Draw("C hist same");

							int fraction = (int)(BreakDownPlots[j-1]->Integral() / Plots[WhichFSIModel]->Integral() * 100.);
							//cout << "Interaction " << j << " fraction = " << fraction << endl; 

							f->cd();
							BreakDownPlots[j-1]->Write(FSILabel[WhichFSIModel]+"_"+InteractionLabels[j]+PlotNames[WhichPlot]);	

						} // end of the look over the GENIE break down

					} // End of the SuSav2 if statement for breakdown

					// --------------------------------------------------------------------------------------------------

					if ( string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

						Plots[WhichFSIModel]->SetMarkerStyle(20); 
						Plots[WhichFSIModel]->SetMarkerSize(2.); 
						Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

						gStyle->SetErrorX(0); // Removing the horizontal errors

						DataPlot = Plots[WhichFSIModel];

						DataPlot = AcceptanceCorrection(Plots[WhichFSIModel],"SuSav2", nucleus[WhichNucleus],E[WhichEnergy],PlotNames[WhichPlot],xBCut);

						DataPlot->SetMarkerStyle(20); 
						DataPlot->SetMarkerSize(2.); 
						DataPlot->SetLineColor(kBlack);	
						DataPlot->SetMarkerColor(kBlack);
						max = DataPlot->GetMaximum();

						DataPlot->GetYaxis()->SetRangeUser(-0.005*max,1.1*max);	
						DataPlot->GetYaxis()->SetTitle( VarLabel[ PlotNames[WhichPlot] ] );

						DataPlot->SetTitle( nucleus[WhichNucleus] + " @" + LabelE[WhichEnergy] + " GeV, " + LatexLabel[ PlotNames[WhichPlot] ]);

						DataPlot->Draw("e same"); 

						f->cd();
						DataPlot->Write("Data_"+PlotNames[WhichPlot]);

					} else { 

						if (FSILabel[WhichFSIModel] == "G2018") { Plots[WhichFSIModel]->SetLineStyle(kDashed); }
						Plots[WhichFSIModel]->Draw("C hist same");  // draw them as lines
						if (FSILabel[WhichFSIModel] == "G2018") { BreakDownPlots[3]->Draw("C hist same"); }						

						if (FSILabel[WhichFSIModel] == "SuSav2") { legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l"); }

						max = TMath::Max(max,Plots[WhichFSIModel]->GetMaximum());
						DataPlot->GetYaxis()->SetRangeUser(-0.005*max,1.1*max);	
						DataPlot->Draw("e same"); 

						f->cd();
						Plots[WhichFSIModel]->Write(FSILabel[WhichFSIModel]+"_"+PlotNames[WhichPlot]);							

					}

					FileSample->Close();

				} // End of the loop over the FSI Models 

				gPad->RedrawAxis();

				legGenie->SetBorderSize(0);
				legGenie->SetTextFont(FontStyle);

				legGenieBlackLine->SetBorderSize(0);
				legGenieBlackLine->SetTextFont(FontStyle);

				legGenieBreak->SetBorderSize(0);
				legGenieBreak->SetTextFont(FontStyle);

				legGenie->SetTextSize(TextSize);

				legGenieBlackLine->SetNColumns(1); 
				legGenieBlackLine->SetTextSize(TextSize-0.03); 

				legGenieBreak->SetTextSize(TextSize-0.03);

				legG2018->SetBorderSize(0);
				legG2018->SetTextFont(FontStyle);
				legG2018->SetTextSize(TextSize-0.03);
				legG2018->AddEntry(Plots[2],"G2018","l");	

				legGenieBlackLine->Draw();
				legGenieBreak->Draw();
				legG2018->Draw();

				// -------------------------------------------------------------------------------------------

				TString CanvasSaveName = "myPlots/"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_" +PlotNames[WhichPlot]+"_XSec";
				PlotCanvas->SaveAs(CanvasSaveName+".pdf");
				//PlotCanvas->SaveAs(CanvasSaveName+".eps");

				//delete PlotCanvas;

				// --------------------------------------------------------------------------------------

			} // End of the loop over the plots

			f->Close();
				
			// --------------------------------------------------------------------------------------				

		} // End of the loop over the nuclei

	} // End of the loop over the energies

} // End of the program

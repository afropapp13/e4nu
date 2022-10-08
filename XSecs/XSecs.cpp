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

#include "../myFunctions.cpp"
#include "../Constants.h"

//------------------------------//

void UpdateXRange(TString Energy, int ThetaBin, TH1D* h) {

	if (Energy == "2_261" && ThetaBin == 1) { h->GetXaxis()->SetRangeUser(0.45,1.75); }
	if (Energy == "2_261" && ThetaBin == 2) { h->GetXaxis()->SetRangeUser(0.55,1.75); }
	if (Energy == "2_261" && ThetaBin == 3) { h->GetXaxis()->SetRangeUser(0.8,1.75); }
	if (Energy == "2_261" && ThetaBin == 4) { h->GetXaxis()->SetRangeUser(0.55,1.75); }			
	if (Energy == "2_261" && ThetaBin == 5) { h->GetXaxis()->SetRangeUser(0.55,1.75); }

}

// ----------------------------------------------------------------------------------------------------------------

void XSecs() {

	//------------------------------//

	std::map<TString,double> MinThetaSlice =
	{
		{ "1161", 27 }, // 36
		{ "2261", 26 }, // 25.5
		{ "4461", 20 } // 19.5
	};
	
	std::map<TString,double> MaxThetaSlice =
	{
		{ "1161", 51 }, // 39
		{ "2261", 50 }, // 28.5
		{ "4461", 44 } // 22.5
	};
	
	std::map<TString,int> ThetaSlices =
	{
		{ "1161", 6 },
		{ "2261", 6 },
		{ "4461", 6 }
	};		

	// ------------------------------------------------------------------------

	GlobalSettings();
	//TGaxis::SetMaxDigits(3);

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut; 
	std::vector<TString> nucleus; 
	std::vector<TString> JustNucleus;
	std::vector<TString> E;
	std::vector<TString> LabelE;	
	std::vector<double> DoubleE;
	std::vector<TString> FSIModel;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots; 
	std::vector<TString> BreakDown;	
	std::vector<TString> OutputPlotNames;
	std::vector<TString> Theta;	

	// ------------------------------------------------------------------------

	//nucleus.push_back("4He"); JustNucleus.push_back("He");
	nucleus.push_back("12C"); JustNucleus.push_back("^{12}C");
	//nucleus.push_back("56Fe"); JustNucleus.push_back("Fe");		

	// ------------------------------------------------------------------------

	//E.push_back("1_161"); DoubleE.push_back(1.161); LabelE.push_back("1.161");
	E.push_back("2_261"); DoubleE.push_back(2.261); LabelE.push_back("2.261");
	//E.push_back("4_461"); DoubleE.push_back(4.461); LabelE.push_back("4.461");

	// ------------------------------------------------------------------------

	xBCut.push_back("NoxBCut");	

	// ------------------------------------------------------------------------

	// Samples

	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Pinned Data");

	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("SuSav2");	
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("G2018");

	// ------------------------------------------------------------------------

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> BreakDownPlots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// TString TxtName = "myXSec/EQEFig2_XSecs.txt";
	// ofstream myTxtFile;
	// myTxtFile.open(TxtName);
	// myTxtFile << std::fixed << std::setprecision(3);
	// myTxtFile << "EQE Fig.2;" << endl << endl;
	// myTxtFile << "Bin #;Low bin edge; High bin edge;XSec [μbarn/GeV];XSec error [μbarn/GeV]" << endl << endl;

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

				TFile* f = TFile::Open("myXSec/XSec_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_GeV_"+xBCut[WhichxBCut]+".root","recreate");	

				// ---------------------------------------------------------------------------------------------------------------------------------------------							

				NameOfPlots.clear();
				OutputPlotNames.clear();
				BreakDown.clear();	

				TString CopyBeam = E[WhichEnergy];
				CopyBeam.ReplaceAll("_","");

				for (int WhichTheta = 0; WhichTheta < ThetaSlices[CopyBeam]; WhichTheta++) {

					double ThetaStep = (MaxThetaSlice[CopyBeam] - MinThetaSlice[CopyBeam] ) / ThetaSlices[CopyBeam];
					double MinTheta = MinThetaSlice[CopyBeam] + WhichTheta*ThetaStep;
					double MaxTheta = MinThetaSlice[CopyBeam] + (WhichTheta+1)*ThetaStep;
					int AveTheta = (MaxTheta + MinTheta) / 2.;		
					TString StringAveTheta = TString(std::to_string(AveTheta));	

					NameOfPlots.push_back("h1_Omega_Theta_"+StringAveTheta+"_Sector_0"); 
					OutputPlotNames.push_back("Omega_Theta_"+StringAveTheta+"_Sector_0"); 
					BreakDown.push_back("_Omega_Theta_"+StringAveTheta+"_Sector_0");
					Theta.push_back(StringAveTheta);					

				}				

				int NPlots = NameOfPlots.size();				

				// Loop over the plots				

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					TString CanvasName = nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut];
					TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
					PlotCanvas->SetBottomMargin(0.15);
					PlotCanvas->SetLeftMargin(0.17);					

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					TLegend* legGenie = new TLegend(0.18,0.45,0.33,0.7);
					TLegend* legGenieBlackLine = new TLegend(0.18,0.74,0.33,0.88);
					TLegend* legGenieBreak = new TLegend(0.21,0.63,0.38,0.75);					
					TLegend* legG2018 = new TLegend(0.18,0.57,0.36,0.63);					

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

						TString PathToFiles = "../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName,"readonly");

						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						// --------------------------------------------------------------------------------------

						// Make the plot pretty

						Plots[WhichFSIModel]->SetLineColor(DataSetColors[WhichFSIModel]);
						PrettyDoubleXSecPlot(Plots[WhichFSIModel]);

						Plots[WhichFSIModel]->GetXaxis()->SetTitleOffset(0.95);
						Plots[WhichFSIModel]->GetXaxis()->SetTitle("(e,e') Energy Transfer [GeV]");
						Plots[WhichFSIModel]->GetXaxis()->CenterTitle(0);

						Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(1.05);
						Plots[WhichFSIModel]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#Omega dE} [#frac{#mub}{sr GeV " + JustNucleus[WhichNucleus] + "}]");						

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

						UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

						// Solid angle division
						Plots[WhichFSIModel]->Scale(1./dOmega);

						// ----------------------------------------------------------------------------------

						// Genie Break Down

						if (
							FSILabel[WhichFSIModel] == "SuSav2"
//							FSILabel[WhichFSIModel] == "G2018"
						) {

							legGenie->AddEntry(Plots[0],"Data", "lep"); 
							legGenieBlackLine->AddEntry(Plots[0],"Data", "lep"); 

							legGenie->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l"); 
//							legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"GENIE (Total)", "l"); 

							BreakDownPlots.clear();

							for (int j = 1; j < 5; j++) {


								BreakDownPlots.push_back( (TH1D*)( FileSample->Get("h1_"+ToStringInt(j)+BreakDown[WhichPlot]) ) );

								UniversalE4vFunction(BreakDownPlots[j-1],FSIModelsToLabels[FSIModel[WhichFSIModel]],\
										     nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

								// Solid angle division
								BreakDownPlots[j-1]->Scale(1./dOmega);											 

								//-----------------------------------------------------------------------------------------------

								BreakDownPlots[j-1]->SetLineColor(BreakDownColors[j-1]);
								BreakDownPlots[j-1]->SetLineWidth(3);

								TLegendEntry* l1 = legGenie->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
								l1->SetTextColor(BreakDownColors[j-1]);

								TLegendEntry* l1Break = legGenieBreak->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
								l1Break->SetTextColor(BreakDownColors[j-1]);

								PlotCanvas->cd();

								//----------------------------------------//

								UpdateXRange(E[WhichEnergy],WhichPlot,BreakDownPlots[j-1]);						

								//----------------------------------------//

								BreakDownPlots[j-1]->Draw("C hist same");

								int fraction = (int)(BreakDownPlots[j-1]->Integral() / Plots[WhichFSIModel]->Integral() * 100.);
								//cout << "Interaction " << j << " fraction = " << fraction << endl; 

								f->cd();
								BreakDownPlots[j-1]->Write(FSILabel[WhichFSIModel]+"_"+BreakDown[WhichPlot]+ToStringInt(j));	


							} // end of the look over the GENIE break down

						} // End of the SuSav2 if statement for breakdown

						// --------------------------------------------------------------------------------------------------

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(2.); 
							Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

							gStyle->SetErrorX(0); // Removing the horizontal errors

							//DataPlot = Plots[WhichFSIModel];
							DataPlot = AcceptanceCorrection(Plots[WhichFSIModel],"SuSav2", nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot],xBCut[WhichxBCut]);

							DataPlot->SetMarkerStyle(20); 
							DataPlot->SetMarkerSize(2.); 
							DataPlot->SetLineColor(kBlack);	
							DataPlot->SetMarkerColor(kBlack);
							max = DataPlot->GetMaximum();
							min = DataPlot->GetMinimum();							
							DataPlot->GetYaxis()->SetRangeUser(min,1.1*max);

							//----------------------------------------//

							UpdateXRange(E[WhichEnergy],WhichPlot,DataPlot);						

							//----------------------------------------//

							DataPlot->SetTitle( LabelE[WhichEnergy] + " GeV, #theta = " + Theta[WhichPlot] + "^{o}");

							PlotCanvas->cd();
							DataPlot->Draw("e same"); 

							f->cd();
							DataPlot->Write("Data_"+NameOfPlots[WhichPlot]);

							// --------------------------------------------------------------------------------------------------				

						} else { 

							if (FSILabel[WhichFSIModel] == "G2018") { Plots[WhichFSIModel]->SetLineStyle(kDashed); }
							PlotCanvas->cd();

							//----------------------------------------//

							UpdateXRange(E[WhichEnergy],WhichPlot,Plots[WhichFSIModel]);						

							//----------------------------------------//

							Plots[WhichFSIModel]->Draw("C hist same");  // draw them as lines
							if (FSILabel[WhichFSIModel] == "G2018") { BreakDownPlots[3]->Draw("C hist same"); }						

							if (FSILabel[WhichFSIModel] == "SuSav2") { legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l"); }

							max = TMath::Max(max,Plots[WhichFSIModel]->GetMaximum());
							DataPlot->GetYaxis()->SetRangeUser(0.,1.1*max);	

							DataPlot->Draw("e same"); 

							f->cd();
							Plots[WhichFSIModel]->Write(FSILabel[WhichFSIModel]+"_"+NameOfPlots[WhichPlot]);							

						}

		                // --------------------------------------------------------------------------------------------------

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

					TString ext = "";
					if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

					TString CanvasSaveName = "myPlots/"+ext+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_" +OutputPlotNames[WhichPlot]+"_XSec";
					PlotCanvas->SaveAs(CanvasSaveName+".pdf");
					PlotCanvas->SaveAs(CanvasSaveName+".eps");

					//delete PlotCanvas;

					// --------------------------------------------------------------------------------------

				} // End of the loop over the plots

				f->Close();
				
				// --------------------------------------------------------------------------------------				

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program

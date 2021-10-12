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

#include "../../myFunctions.cpp"
#include "../../AfroConstants.h"

// ----------------------------------------------------------------------------------------------------------------

void AbsXSec_OverlayPmissFig3a_e4nuPaper() {

	// ------------------------------------------------------------------------

	GlobalSettings();

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut; 
	std::vector<TString> nucleus; 
	std::vector<TString> JustNucleus;
	std::vector<TString> E; 
	std::vector<double> DoubleE;
	std::vector<TString> LabelE; 
	std::vector<TString> FSIModel;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots; 
	std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames;

	// ------------------------------------------------------------------------

	nucleus.push_back("12C"); JustNucleus.push_back("C");

	E.push_back("2_261"); LabelE.push_back("2.257 GeV"); DoubleE.push_back(2.261);	

	xBCut.push_back("NoxBCut");

	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Pinned Data");

	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc"); FSILabel.push_back("SuSav2");
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc"); FSILabel.push_back("Genie");

//	FSIModel.push_back("SuSav2_RadCorr_LFGM"); FSILabel.push_back("SuSav2");
//	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("Genie");

	NameOfPlots.push_back("MissMomentum"); LabelOfPlots.push_back("P_{T} [GeV/c]"); OutputPlotNames.push_back("MissMomentum");

	// ------------------------------------------------------------------------

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> BreakDownPlots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				// Loop over the plots

				double max = -99.;
				double min = 1E12;
				double height = 1.1;

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					TCanvas* PlotCanvas = new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 205,34,850,1000);

					// ----------------------------------------------------------------------------------------------

					// Dimensions of TPads (pad2 will be deleted at the very end for the Ereco plots)

					double XMinPadOne = 0., XMaxPadOne = 1., YMinPadOne = 0.15, YMaxPadOne = 0.7;
					double XMinPadTwo = 0., XMaxPadTwo = 1., YMinPadTwo = YMaxPadOne+0.01, YMaxPadTwo = 0.85;

					// ----------------------------------------------------------------------------------------

					TPad* pad1 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadOne,YMinPadOne,XMaxPadOne,YMaxPadOne, 21); 
					pad1->SetFillColor(kWhite); pad1->Draw();
					TPad* pad2 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadTwo,YMinPadTwo,XMaxPadTwo,YMaxPadTwo,22); 
					pad2->SetFillColor(kWhite); pad2->Draw(); 

					pad1->SetTopMargin(0.);
					pad1->SetBottomMargin(0.19);
					pad1->SetLeftMargin(0.17);
					pad1->SetRightMargin(0.04);

					pad2->SetTopMargin(0.);
					pad2->SetBottomMargin(0.);
					pad1->cd();

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					TLegend* legData = new TLegend(0.6,0.82,0.9,0.99);					
					legData->SetNColumns(1);

					TLegend* legGenie = new TLegend(0.67,0.65,0.95,0.82);					
					legGenie->SetNColumns(2);

					TLegend* legG2018 = new TLegend(0.6,0.55,0.9,0.65);					
					legG2018->SetNColumns(1);

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						// --------------------------------------------------------------------------------------

						// Make the plot pretty

						Plots[WhichFSIModel]->SetLineColor(DataSetColors[WhichFSIModel]);
						PrettyDoubleXSecPlot(Plots[WhichFSIModel]);
						Plots[WhichFSIModel]->GetXaxis()->SetTitle(LabelOfPlots[WhichPlot]);

						// ----------------------------------------------------------------------------------

						// Use the universal e4v function that 
						// scales by all the necessary scaling factors
						// divides by the bin width to obtain normalized yields / absolute cross sections
						// uses the relevant binning
						// gets the relevant x axis range
						// If data sample: 
						//                 apply systematics due to rotations et al
						//                 apply acceptance systematics using sector-by -sector uncertainties

						UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

						// ----------------------------------------------------------------------------------

						// Genie Break Down

						if (
							FSILabel[WhichFSIModel] == "SuSav2"
						) {

							legData->AddEntry(Plots[0],"Data", "lep"); 
//							legGenie->AddEntry(Plots[WhichFSIModel],"GENIE (Total)", "l"); 

							BreakDownPlots.clear();

							for (int j = 1; j < 5; j++) {

								BreakDownPlots.push_back( (TH1D*)( FileSample->Get("Pmiss_Int_"+ToStringInt(j)) ) );

								UniversalE4vFunction(BreakDownPlots[j-1],FSIModelsToLabels[FSIModel[WhichFSIModel]],\
										     nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

								//-----------------------------------------------------------------------------------------------

								BreakDownPlots[j-1]->SetLineColor(BreakDownColors[j-1]);
								
								BreakDownPlots[j-1]->SetLineWidth(3);
								//BreakDownPlots[j-1]->SetLineStyle(Style[j-1]);

								TLegendEntry* l1 = legGenie->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
								if (j == 2) { legData->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l");  }
								//if (j == 2) { legGenie->AddEntry(UncertaintyPlots[1],"GENIE (Total)", "lf");  }
								l1->SetTextColor(BreakDownColors[j-1]);

								pad1->cd();
								BreakDownPlots[j-1]->Draw("C hist same");

							} // end of the look over the GENIE break down

						} // End of the SuSav2 / break down option

						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

						double localmax = Plots[WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }
						Plots[0]->GetYaxis()->SetRangeUser(0.,height*max);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

//						Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(5); 
//						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);

						// --------------------------------------------------------------------------------------------------

						// Drawing The Plots

						pad1->cd();

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(2.); 
							Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

							gStyle->SetErrorX(0); // Removing the horizontal errors
							Plots[WhichFSIModel]->Draw("e same"); 

						} else { 
						
							if (FSILabel[WhichFSIModel] =="Genie") { Plots[WhichFSIModel]->SetLineStyle(kDashed); }
							Plots[WhichFSIModel]->Draw("C hist same");  // draw them as lines

							Plots[0]->Draw("e same"); 

						}

					} // End of the loop over the FSI Models 

					// -----------------------------------------------------------------------------------

					// Transverse Missing Momentum Slices // Vertical Lines

					pad1->cd();
					double LowPmiss = 0.2, MidPmiss = 0.4, HighPmiss = 1.;

					TLine* line1 = new TLine(LowPmiss,0.,LowPmiss,height*max);
					line1->SetLineColor(kBlack); line1->SetLineWidth(1);
					line1->Draw(); 

					TLine* line2 = new TLine(MidPmiss,0.,MidPmiss,height*max);
					line2->SetLineColor(kBlack); line2->SetLineWidth(1);
					line2->Draw(); 

					legGenie->SetBorderSize(0);
					legGenie->SetTextFont(FontStyle);

					legData->SetBorderSize(0);
					legData->SetTextFont(FontStyle);
					legData->SetTextSize(TextSize);
					legData->Draw();

					//pad2->cd(); 
					legGenie->SetTextSize(TextSize);
					//legGenie->AddEntry(Plots[2],"G2018", "l"); 
					legGenie->Draw(); 

					legG2018->SetBorderSize(0);
					legG2018->SetTextFont(FontStyle);
					legG2018->SetTextSize(TextSize);
					legG2018->AddEntry(Plots[2],"G2018", "l"); 
					legG2018->Draw();

					// ---------------------------------------------------------------------------------------------------

					TLatex* myNucleus = new TLatex();
					myNucleus->SetTextFont(FontStyle);
					myNucleus->SetTextColor(kBlack);
					myNucleus->SetTextSize(3.5*TextSize);
					pad2->cd();
					myNucleus->DrawLatexNDC(0.43,0.05,LabelE[WhichEnergy]);

					// ---------------------------------------------------------------------------------------------------------

					TString ext = "";
					if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

					PlotCanvas->SaveAs("../../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"+E[WhichEnergy]+"/"
						+ext+"Fig3a_"+nucleus[WhichNucleus]+"_" 
						+E[WhichEnergy]+"_" +OutputPlotNames[WhichPlot]+"_SuSav2_AbsXSec.pdf");

					//delete PlotCanvas;

					// ------------------------------------------------------------------------------


				} // End of the loop over the plots

				// --------------------------------------------------------------------------------------				

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program

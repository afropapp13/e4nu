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

void AbsXSec_OverlayDeltaAlphaT_FigExtData8() {

	// ------------------------------------------------------------------------

	GlobalSettings();
	TGaxis::SetMaxDigits(4);

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut; 
	std::vector<TString> nucleus;
	std::vector<double> Energy;
	std::vector<TString> FSIModel;
	std::vector<TString> NameOfPlots; 

	// ------------------------------------------------------------------------

	nucleus.push_back("12C");
	nucleus.push_back("56Fe");

	Energy.push_back(1.161);
	Energy.push_back(2.261);	
	Energy.push_back(4.461);

	xBCut.push_back("NoxBCut");

	// ------------------------------------------------------------------------

	FSIModel.push_back("Pinned_Data_Final");
	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc");
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc");

//	FSIModel.push_back("SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger");
//	FSIModel.push_back("hA2018_Final_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger");

//	FSIModel.push_back("Pinned_Data_Final_SixSectors");
//	FSIModel.push_back("SuSav2_RadCorr_LFGM_SixSectors");
//	FSIModel.push_back("hA2018_Final_RadCorr_LFGM_SixSectors");

	NameOfPlots.push_back("DeltaAlphaT_Int_0");

	// ------------------------------------------------------------------------

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> BreakDownPlots;

	// ------------------------------------------------------------------------

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = Energy.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	// ------------------------------------------------------------------------

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		TCanvas* PlotCanvas = new TCanvas(xBCut[WhichxBCut],xBCut[WhichxBCut],205,34,1600,900);

		TLegend* legGenieBlackLine = new TLegend(0.1,0.5,0.5,1.);
		legGenieBlackLine->SetNColumns(1);
		legGenieBlackLine->SetTextFont(FontStyle); 
					
		TLegend* legGenieBreak = new TLegend(0.2,0.21,0.85,0.53);
		legGenieBreak->SetNColumns(2);
		legGenieBreak->SetTextFont(FontStyle);

		TLegend* legG2018 = new TLegend(0.1,0.0,1.,0.1);
		legG2018->SetNColumns(2);
		legG2018->SetTextFont(FontStyle);

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			// Loop over the energies

			for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

				// ---------------------------------------------------------------------------------

				// Energy from double to TString to get path to files

				TString E = ((TString)(std::to_string(Energy[WhichEnergy])));
				E.ReplaceAll(".","_");
				E.ReplaceAll("0","");

				// ---------------------------------------------------------------------------------

	 			// In order to use y-axis ticks with common scale, constraint range between (0,MaxHeight)
			
				double MaxHeight = 0.49;

				// Loop over the nuclei

				for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

					if (nucleus[WhichNucleus] == "56Fe") { MaxHeight = 2.3; }

					// ---------------------------------------------------------------------------------------------------------------

					// Create TPad for specific nucleus & energy

					TPad* pad = CreateTPad(WhichEnergy,WhichNucleus,Energy[WhichEnergy],nucleus[WhichNucleus],NameOfPlots[WhichPlot],PlotCanvas);

					// -------------------------------------------------------------------------------------------------------

					// No data on 56Fe @ 1.161 GeV

					if ( nucleus[WhichNucleus] == "56Fe" && Energy[WhichEnergy] == 1.161 ) { delete pad; continue; }

					// --------------------------------------------------------------------------------------------------------

					Plots.clear();

					double max = -99.;
					double min = 1E12;				

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						// --------------------------------------------------------------------------------------

						TString PathToFiles = "../../../myFiles/"+ E + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus] + "_" + E +"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						// --------------------------------------------------------------------------------------

						// Make the plot pretty

						PrettyDoubleXSecPlot(Plots[WhichFSIModel]);
						Plots[WhichFSIModel]->SetLineColor(DataSetColors[WhichFSIModel]);
						Plots[WhichFSIModel]->SetLineWidth(1);

						Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(1.2*TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(0.);
						Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions);

						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(0.);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelOffset(0.013);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(1.2*TextSize);
						Plots[WhichFSIModel]->GetYaxis()->SetTitle("");
						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(4);

						// --------------------------------------------------------------------------------------

						// Use the universal e4v function that 
						// scales by all the necessary scaling factors
						// divides by the bin width to obtain normalized yields / absolute cross sections
						// uses the relevant binning
						// gets the relevant x axis range
						// If data sample: 
						//                 apply systematics due to rotations et al
						//                 apply acceptance systematics using sector-by -sector uncertainties

						UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E,NameOfPlots[WhichPlot]);
						Plots[WhichFSIModel]->Scale(1000);

						// ----------------------------------------------------------------------------------

						// Genie Break Down

						if (
							FSIModelsToLabels[FSIModel[WhichFSIModel]] == "SuSav2 NoRad"
						) {

							if (Energy[WhichEnergy] == 1.161 && nucleus[WhichNucleus] == "12C") {
							
								legGenieBlackLine->AddEntry(Plots[0],"Data", "lep"); 
								legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l"); 
							}

							BreakDownPlots.clear();

							for (int j = 1; j < 5; j++) {

								BreakDownPlots.push_back( (TH1D*)( FileSample->Get("DeltaAlphaT_Int_"+ToStringInt(j)) ) );

								BreakDownPlots[j-1]->SetLineColor(BreakDownColors[j-1]);
								//BreakDownPlots[j-1]->SetLineWidth(LineWidth);

								UniversalE4vFunction(BreakDownPlots[j-1],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E,NameOfPlots[WhichPlot]);
								BreakDownPlots[j-1]->Scale(1000);

								//-----------------------------------------------------------------------------------------------

								if (Energy[WhichEnergy] == 1.161 && nucleus[WhichNucleus] == "12C") {

									TLegendEntry* l1Break = legGenieBreak->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
									l1Break->SetTextColor(BreakDownColors[j-1]);
								}

								BreakDownPlots[j-1]->Draw("C hist same");

							} // end of the look over the GENIE break down

						} // End of the SuSav2 if statement for breakdown

						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

						Plots[0]->GetYaxis()->SetRangeUser(0.,MaxHeight);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

						// -------------------------------------------------------------------------------------------------

						TLine* line = new TLine(0.95*Energy[WhichEnergy],0.,0.95*Energy[WhichEnergy],MaxHeight);
						line->SetLineColor(kBlack); 
						line->SetLineWidth(LineWidth);

						// --------------------------------------------------------------------------------------------------

						if (string(FSIModelsToLabels[FSIModel[WhichFSIModel]]).find("Data") != std::string::npos) { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(1.5); 
							Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

							gStyle->SetErrorX(0); // Removing the horizontal errors
							Plots[WhichFSIModel]->Draw("e same"); 

						} else { 

							Plots[WhichFSIModel]->SetLineStyle(Style[WhichFSIModel]);
							Plots[WhichFSIModel]->Draw("C hist same");  // "C hist same" draw them as lines // "hist same" draw them as histos
							Plots[0]->Draw("e same"); 

						}

					} // End of the loop over the FSI Models 					

					// --------------------------------------------------------------------------------------------------------

				} // End of the loop over the energies

			} // End of the loop over the nuclei

		} // End of the loop over the plots

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TLatex latexCarbon;
		latexCarbon.SetTextFont(FontStyle);
		latexCarbon.SetTextSize(0.8*TextSize);
		PlotCanvas->cd();
		latexCarbon.DrawLatexNDC(0.005,0.77,"^{12}C");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TLatex latexIron;
		latexIron.SetTextFont(FontStyle);
		latexIron.SetTextSize(0.8*TextSize);
		PlotCanvas->cd();
		latexIron.DrawLatexNDC(0.3,0.385,"^{56}Fe");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* pad1GeV = new TPad("pad1GeV","pad1GeV",0.18,0.89,0.34,0.99,21);
		pad1GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad1GeV->Draw();
		pad1GeV->cd();

		TLatex latex1GeV;
		latex1GeV.SetTextFont(FontStyle);
		latex1GeV.SetTextSize(6*TextSize);
		latex1GeV.DrawLatexNDC(0.03,0.45,"1.159 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* pad2GeV = new TPad("pad2GeV","pad2GeV",0.47,0.89,0.62,0.99,21); 
		pad2GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad2GeV->Draw();
		pad2GeV->cd();

		TLatex latex2GeV;
		latex2GeV.SetTextFont(FontStyle);
		latex2GeV.SetTextSize(6*TextSize);
		latex2GeV.DrawLatexNDC(0.06,0.45,"2.257 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* pad4GeV = new TPad("pad4GeV","pad4GeV",0.75,0.89,0.92,0.99,21); 
		pad4GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad4GeV->Draw();
		pad4GeV->cd();

		TLatex latex4GeV;
		latex4GeV.SetTextFont(FontStyle);
		latex4GeV.SetTextSize(6*TextSize);
		latex4GeV.DrawLatexNDC(0.11,0.45,"4.453 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* padPmiss = new TPad("padPmiss","padPmiss",0.405,0.01,0.765,0.14,21); 
		padPmiss->SetFillColor(kWhite); 
		PlotCanvas->cd();
		padPmiss->Draw();
		padPmiss->cd();

		TLatex latexPmiss;
		latexPmiss.SetTextFont(FontStyle);
		latexPmiss.SetTextSize(5*TextSize);
		latexPmiss.DrawLatexNDC(0.2,0.5,"(e,e'p)_{1p0#pi} #delta#alpha_{T} [deg]");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the legend

		PlotCanvas->cd();
		TPad* padLegend = new TPad("padLegend","padLegend",0.07,0.19,0.26,0.53, 21); 
		padLegend->SetFillColor(kWhite); 
		padLegend->Draw();
		padLegend->cd();

		legGenieBlackLine->SetTextSize(2.*TextSize); 
		legGenieBlackLine->SetBorderSize(0); 
		legGenieBlackLine->Draw();

		//legGenieBreak->AddEntry(Plots[2],"G2018", "l");
		legGenieBreak->SetTextSize(2.*TextSize); 
		legGenieBreak->SetBorderSize(0); 
		legGenieBreak->Draw();

		legG2018->AddEntry(Plots[2],"G2018", "l");
		legG2018->SetTextSize(2.*TextSize); 
		legG2018->SetBorderSize(0); 
		legG2018->Draw();	

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis units carbon

		PlotCanvas->cd();
		TPad* padTitle = new TPad("padTitle","padTitle",0.062,0.55,0.087,0.93, 21); 
		padTitle->SetFillColor(kWhite); 
		padTitle->Draw();
		padTitle->cd();

		TLatex latexYTitle;
		latexYTitle.SetTextFont(FontStyle);
		latexYTitle.SetTextSize(14*TextSize);
		latexYTitle.SetTextColor(kBlack);
		latexYTitle.SetTextAngle(90);
//		latexYTitle.DrawLatexNDC(0.82,0.1,DoubleXSecTitle+" x10 ^{3}");
		latexYTitle.DrawLatexNDC(0.82,0.2,DoubleXSecTitle);

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis units iron

		PlotCanvas->cd();
		TPad* padTitleFe = new TPad("padTitleFe","padTitleFe",0.37,0.15,0.395,0.555,21); 
		padTitleFe->SetFillColor(kWhite); 
		padTitleFe->Draw();
		padTitleFe->cd();

		TLatex latexYTitleFe;
		latexYTitleFe.SetTextFont(FontStyle);
		latexYTitleFe.SetTextSize(14*TextSize);
		latexYTitleFe.SetTextColor(kBlack);
		latexYTitleFe.SetTextAngle(90);
//		latexYTitleFe.DrawLatexNDC(0.82,0.05,DoubleXSecTitle+" x10 ^{3}");
		latexYTitleFe.DrawLatexNDC(0.82,0.15,DoubleXSecTitle);

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the lower X-axis to cover half zeros

		PlotCanvas->cd();
		TPad* padWhitePadTwo = new TPad("padWhitePadTwo","padWhitePadTwo",0.713,0.15,0.72,0.185,21); 
		padWhitePadTwo->SetFillColor(kWhite);
		padWhitePadTwo->Draw();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString ext = "";
		if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

		PlotCanvas->SaveAs("../../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+ext+"DeltaAlphaT_FigExtData8_SuSav2_AbsXSec.pdf");

		//delete PlotCanvas;

	} // End of the loop over the xB kinematics

} // End of the program

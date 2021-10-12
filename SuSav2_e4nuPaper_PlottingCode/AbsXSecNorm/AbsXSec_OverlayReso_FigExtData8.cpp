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

void AbsXSec_OverlayReso_FigExtData8() {

	// ------------------------------------------------------------------------

	GlobalSettings();

	// ------------------------------------------------------------------------

	int Ndivisions = 5;
	int LineWidth = 5;
	double TextSize = 0.06;

	// ------------------------------------------------------------------------

	TGaxis::SetMaxDigits(3);
	TGaxis::SetExponentOffset(-0.1, 1., "y");
	gStyle->SetLineStyleString(11,"80 60");

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> LabelsOfSamples; std::vector<TString> E;
	std::vector<TString> LabelE; std::vector<TString> JustE; std::vector<TString> FSIModel; std::vector<TString> OutputPlotNames;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots;  std::vector<TString> XLabels;

	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C");
	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");

	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV"); JustE.push_back("1.159 GeV");
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV"); JustE.push_back("2.257 GeV");
	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV"); JustE.push_back("4.453 GeV");
 
	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Data");	
//	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc"); FSILabel.push_back("SuSav2"); // We need to use the NoRad SuSav2 sample, given issue in G2018 Rad sample		
	FSIModel.push_back("SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger"); FSILabel.push_back("SuSav2");		

	xBCut.push_back("NoxBCut");
	NameOfPlots.push_back("h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed"); OutputPlotNames.push_back("EcalReso");
//	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_factor_fracfeed"); OutputPlotNames.push_back("EQEReso"); // 1p0pi
	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_noprot_frac_feed3pi"); OutputPlotNames.push_back("EQEReso"); // 0pi

	XLabels.push_back("(e,e'p)_{1p0#pi} E_{cal} Feeddown"); 
	XLabels.push_back("(e,e')_{0#pi} E_{QE} Feeddown"); 

	// -------------------------------------------------------------------------------------------------

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	const int NEnergies = E.size();
	const int NFSIModels = FSIModel.size();
	const int NPlots = NameOfPlots.size();

	TH1D* Plots[NEnergies][NFSIModels];

	// -------------------------------------------------------------------------------------------------

	vector< vector<int> > Colors{{kRed,kRed},{kGreen-3,kGreen-3},{kBlue,kBlue}};
	vector<int> LineStyle{2,11,1};
	vector<int> MarkerStyle{22,21,20};

	// -------------------------------------------------------------------------------------------------

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the nuclei

		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

				TString CanvasName = nucleus[WhichNucleus]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut];
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);

				PlotCanvas->SetBottomMargin(0.16);
				PlotCanvas->SetLeftMargin(0.11);

				TLegend* leg = new TLegend(0.12,0.69,0.4,0.835);
				if (nucleus[WhichNucleus] == "56Fe") { leg = new TLegend(0.12,0.73,0.4,0.835); }
				leg->SetNColumns(2);
				double max = -99.;

				// Loop over the energies
				
				int LowBin = -1;
				int HighBin = -1;

				if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_frac_feed3pi" || NameOfPlots[WhichPlot] == "h1_E_rec_0pi_frac_feed") 
						{ FSIModel[1] = "SuSav2_RadCorr_LFGM_Truth0pi_WithFidAcc_UpdatedSchwinger"; }
				else { FSIModel[1] = "SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger"; }				

				for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

					if (nucleus[WhichNucleus] == "56Fe" && E[WhichEnergy] == "1_161") { continue; }

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots[WhichEnergy][WhichFSIModel] =  (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) );
						Plots[WhichEnergy][WhichFSIModel]->SetLineColor(Colors[WhichEnergy][WhichFSIModel]);
						Plots[WhichEnergy][WhichFSIModel]->SetLineWidth(LineWidth-WhichEnergy);

						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetLabelSize(TextSize);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetTitleSize(TextSize);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetTitleOffset(1.05);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetNdivisions(7);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->CenterTitle();
						
//						double LowRange = -1;
//						double HighRange = -1;	
//					
//						if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed") { 
//							LowRange = -0.85; HighRange = 0.2; 
//						} else { LowRange = -0.81; HighRange = 0.07; }

//						
//						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetRangeUser(LowRange,HighRange);
//						
//						LowBin = Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->FindBin(LowRange);
//						HighBin = Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->FindBin(HighRange);

						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTickSize(0.02);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTitle("Normalized Yield");
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->CenterTitle();
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTitleOffset(1.);

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { leg->AddEntry(Plots[WhichEnergy][WhichFSIModel],"/", "ep"); }
						else { leg->AddEntry(Plots[WhichEnergy][WhichFSIModel]," "+JustE[WhichEnergy], "l"); }

						// ---------------------------------------------------------------------------------------------------

						// Absolute XSec Scaling 
						// Then use myFunctions
						// Division by bin width
						// Systematics for data

//						double ScalingFactor = 1. / Plots[WhichEnergy][WhichFSIModel]->Integral();
//						Plots[WhichEnergy][WhichFSIModel]->Scale(ScalingFactor);

						AbsoluteXSecScaling(Plots[WhichEnergy][WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy]);
						ReweightPlots(Plots[WhichEnergy][WhichFSIModel]);
						ApplyRange(Plots[WhichEnergy][WhichFSIModel],E[WhichEnergy],NameOfPlots[WhichPlot]);
						if (string(FSIModelsToLabels[FSIModel[WhichFSIModel]]).find("Data") != std::string::npos) 
							{ ApplySystUnc(Plots[WhichEnergy][WhichFSIModel], E[WhichEnergy]); }
						if (string(FSIModelsToLabels[FSIModel[WhichFSIModel]]).find("Data") != std::string::npos) 
							{ ApplySectorSystUnc(Plots[WhichEnergy][WhichFSIModel], E[WhichEnergy]); }

						// ---------------------------------------------------------------------------------------------------

						double localmax = Plots[WhichEnergy][WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }

						TString XLabel = Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[WhichEnergy][0]->GetXaxis()->SetTitle(XLabels[WhichPlot]);

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

							Plots[WhichEnergy][WhichFSIModel]->SetMarkerSize(2.); 
							Plots[WhichEnergy][WhichFSIModel]->SetMarkerColor(Colors[WhichEnergy][WhichFSIModel]);
							gStyle->SetErrorX(0);
							Plots[WhichEnergy][WhichFSIModel]->SetMarkerStyle(MarkerStyle[WhichEnergy]); 
							Plots[WhichEnergy][WhichFSIModel]->Draw("e same"); 

						}
						else { 

							double MaxHeight = 1.3;

							if (nucleus[WhichNucleus] == "12C") { Plots[0][0]->GetYaxis()->SetRangeUser(0.01,MaxHeight*max); }
							else { Plots[1][0]->GetYaxis()->SetRangeUser(0.01,MaxHeight*max); }
							Plots[WhichEnergy][WhichFSIModel]->SetLineStyle(LineStyle[WhichEnergy]);
							Plots[WhichEnergy][WhichFSIModel]->Draw("C hist same"); 

						}

					} // End of the loop over the FSI Models 
					
					// --------------------------------------------------------------------------------------				
					
					// Chi2 calculation
					
					//int NBinsX = HighBin - LowBin +1;
					//int Chi2Double = Chi2(Plots[WhichEnergy][0],Plots[WhichEnergy][1],LowBin,HighBin);
					
					//cout << endl << endl << nucleus[WhichNucleus] << "  " << E[WhichEnergy]  << " Chi2/ndof = " << Chi2Double;
					//cout << "/" << NBinsX << endl << endl;
					
					// --------------------------------------------------------------------------------------				

					// Integrated XSec Ratio Data / SuSav2

					double DataXSec = IntegratedXSec(Plots[WhichEnergy][0]);
					double SuSav2XSec = IntegratedXSec(Plots[WhichEnergy][1]);
					double DataSuSav2Ratio = DataXSec / SuSav2XSec;

					std::cout << nucleus[WhichNucleus] << "  " << E[WhichEnergy] << "  " << OutputPlotNames[WhichPlot] << "  " << DataSuSav2Ratio << std::endl;

					// --------------------------------------------------------------------------------------				

				} // End of the loop over the energies

				leg->SetBorderSize(0);
				leg->SetTextFont(FontStyle);
				leg->SetTextSize(TextSize);
				leg->Draw();

				TLatex latex;
				latex.SetTextFont(FontStyle);
				latex.SetTextSize(TextSize);
				latex.DrawLatexNDC(0.8,0.83,LabelsOfSamples[WhichNucleus]);

				TLatex latexDG;
				latexDG.SetTextFont(FontStyle);
				latexDG.SetTextSize(TextSize);
				latexDG.DrawLatexNDC(0.13,0.855,"Data/SuSav2");

				PlotCanvas->SaveAs("../../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/AbsXSec_FeedDown_"+nucleus[WhichNucleus]+"_" +OutputPlotNames[WhichPlot]+"_SuSav2.pdf");

			} // End of the loop over the plots

		} // End of the loop over the nuclei

	} // End of the loop over the xB kinematics

} // End of the program

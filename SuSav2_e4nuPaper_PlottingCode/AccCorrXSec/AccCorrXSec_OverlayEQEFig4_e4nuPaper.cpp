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

void AccCorrXSec_OverlayEQEFig4_e4nuPaper() {

	// ------------------------------------------------------------------------

	GlobalSettings();
	double height = 1.05;

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut; 
	std::vector<TString> nucleus;
	std::vector<TString> E; 
	std::vector<double> DoubleE;
	std::vector<TString> FSIModel;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots; 
	std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames; 
	std::vector<TH1D*> BreakDownPlots;
	std::vector<int> Colors;

	nucleus.push_back("12C");
	nucleus.push_back("56Fe");

	E.push_back("1_161"); DoubleE.push_back(1.161);
	E.push_back("2_261"); DoubleE.push_back(2.261);	
	E.push_back("4_461"); DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");
 
	Colors.push_back(kBlack); Colors.push_back(kBlack); Colors.push_back(kBlack); 
	Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);

	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Data");
	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc"); FSILabel.push_back("SuSav2");
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc"); FSILabel.push_back("G2018");

//	FSIModel.push_back("Pinned_Data_Final_SixSectors"); FSILabel.push_back("Pinned Data");
//	FSIModel.push_back("SuSav2_RadCorr_LFGM_SixSectors"); FSILabel.push_back("SuSav2");
//	FSIModel.push_back("hA2018_Final_RadCorr_LFGM_SixSectors"); FSILabel.push_back("G2018");

	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_noprot_3pi"); LabelOfPlots.push_back("(e,e')_{0#pi} E_{QE} [GeV]");
	OutputPlotNames.push_back("EQE0pi_slice_0");

	std::vector<TH1D*> Plots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	// ---------------------------------------------------------------------------------------------------------------------------------------------


	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		TFile* fXSecFile = TFile::Open("myXSec/EQE0pi_AccCorr_XSec_"+xBCut[WhichxBCut]+".root","recreate");

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

				// In order to use y-axis ticks with common scale, constraint range between (0,MaxHeight)
				double MaxHeight = 0.56;

				// Loop over the nuclei

				for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

					if (nucleus[WhichNucleus] == "56Fe") { MaxHeight = 2.3; }

					// ----------------------------------------------------------------------------

					// Create TPad for EQE panel

					TPad* pad = CreateTPad(WhichEnergy,WhichNucleus,DoubleE[WhichEnergy],nucleus[WhichNucleus],NameOfPlots[WhichPlot],PlotCanvas);

					// -------------------------------------------------------------------------------

					// No data on 56Fe @ 1.161 GeV

					if ( nucleus[WhichNucleus] == "56Fe" && DoubleE[WhichEnergy] == 1.161 ) { delete pad; continue; }
//					if ( nucleus[WhichNucleus] == "56Fe" && DoubleE[WhichEnergy] == 4.461 ) { delete pad; continue; }
//					if ( nucleus[WhichNucleus] == "12C" && DoubleE[WhichEnergy] == 4.461 ) { delete pad; continue; }

					// ---------------------------------------------------------------------------------------------

					Plots.clear();
//					UncertaintyPlots.clear();


					double max = -99.;
					double min = 1E12;

					// Loop over the FSI Models
					
					int LowBin = -1;
					int HighBin = -1;

					TH1D* DataPlot = nullptr;					

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						// --------------------------------------------------------------------------------------

						// Make the plot pretty

						PrettyDoubleXSecPlot(Plots[WhichFSIModel]);
						Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);

						Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(1.2*TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(0.);

						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(0.);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelOffset(0.013);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(1.2*TextSize);

						Plots[WhichFSIModel]->GetYaxis()->SetTitle("");

						// -----------------------------------------------------------------------------------

						// Use the universal e4v function that 
						// scales by all the necessary scaling factors
						// divides by the bin width to obtain normalized yields / absolute cross sections
						// uses the relevant binning
						// gets the relevant x axis range
						// If data sample: 
						//                 apply systematics due to rotations et al
						//                 apply acceptance systematics using sector-by -sector uncertainties

						UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);
						Plots[WhichFSIModel]->SetLineWidth(1);
						Plots[WhichFSIModel]->GetYaxis()->SetTitle(DoubleAccCorrXSecTitle);

						// --------------------------------------------------------------------------------------				

						// Store extracted xsec

						fXSecFile->cd();
						TString XSecName = OutputPlotNames[WhichPlot]+"_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSILabel[WhichFSIModel];
						Plots[WhichFSIModel]->Write(XSecName);

						//-----------------------------------------------------------------------------------------------

						//Larry's suggestion because ECal has a sharp peak and a low tail 
						//Thus we multiply the peak by EnhaceTail

//						if ( DoubleE[WhichEnergy] == 2.261 ) {

//							double LowE = 0.95*DoubleE[WhichEnergy];
//							int LowEBin = Plots[WhichFSIModel]->FindBin(LowE);
//							int HighEBin = Plots[WhichFSIModel]->GetNbinsX();

//							for (int i = LowEBin + 1; i <= HighEBin; i++) { 
//					
//								double content = Plots[WhichFSIModel]->GetBinContent(i);
//								double error = Plots[WhichFSIModel]->GetBinError(i);
//								double newcontent = EnhaceTail * content;
//								double newerror = EnhaceTail * error;				
//								Plots[WhichFSIModel]->SetBinContent(i,newcontent);
//								Plots[WhichFSIModel]->SetBinError(i,newerror);

//							}

//						}

						// -----------------------------------------------------------------------------------

//						// Accounting for the fact that the bin width might not be constant

//						ReweightPlots(Plots[WhichFSIModel]);

						if (DoubleE[WhichEnergy] == 1.161) { Plots[WhichFSIModel]->Scale(1./2.); }
						if (DoubleE[WhichEnergy] == 4.461) { Plots[WhichFSIModel]->Scale(5.); }

						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(5);
						Plots[WhichFSIModel]->GetYaxis()->SetTitleSize(0);

						// ----------------------------------------------------------------------------------

						// Genie Break Down

						if (
							FSILabel[WhichFSIModel] == "SuSav2"
						) {
						
							if (nucleus[WhichNucleus] == "12C" && E[WhichEnergy] == "1_161") {
								legGenieBlackLine->AddEntry(Plots[0],"Data", "lep"); 
								//legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"GENIE (Total)", "l"); 
							}

							BreakDownPlots.clear();

							for (int j = 1; j < 5; j++) {

								BreakDownPlots.push_back( (TH1D*)( FileSample->Get("InclusiveEQE_Int_"+ToStringInt(j)) ) );

								UniversalE4vFunction(BreakDownPlots[j-1],FSIModelsToLabels[FSIModel[WhichFSIModel]],\
										     nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

								//-----------------------------------------------------------------------------------------------

								BreakDownPlots[j-1]->SetLineColor(BreakDownColors[j-1]);
								
								BreakDownPlots[j-1]->SetLineWidth(1);

								if (DoubleE[WhichEnergy] == 1.161) { BreakDownPlots[j-1]->Scale(1./2.); }
								if (DoubleE[WhichEnergy] == 4.461) { BreakDownPlots[j-1]->Scale(5.); }

								//-----------------------------------------------------------------------------------------------

								if (nucleus[WhichNucleus] == "12C" && E[WhichEnergy] == "1_161") {
								
									TLegendEntry* l1Break = legGenieBreak->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
									l1Break->SetTextColor(BreakDownColors[j-1]);
								}

									BreakDownPlots[j-1]->Draw("C hist same");						
								
							} // end of the look over the GENIE break down

						} // End of the SuSav2 if statement for breakdown

						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

						Plots[0]->GetYaxis()->SetRangeUser(-0.01,MaxHeight);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

						// -------------------------------------------------------------------------------------------------

//						TLine* line = new TLine(0.95*DoubleE[WhichEnergy],0.,0.95*DoubleE[WhichEnergy],MaxHeight);
//						line->SetLineColor(kBlack); 
//						line->SetLineWidth(LineWidth);
//						if ( FSILabel[WhichFSIModel] == "SuSav2" && DoubleE[WhichEnergy] == 2.261) { line->Draw(); }

						// --------------------------------------------------------------------------------------------------

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(1.5); 
							Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

							gStyle->SetErrorX(0); // Removing the horizontal errors
							//Plots[WhichFSIModel]->Draw("e same"); 

							DataPlot = Plots[WhichFSIModel];

							DataPlot = AcceptanceCorrection(Plots[WhichFSIModel],"SuSav2", nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot],xBCut[WhichxBCut]);

							DataPlot->SetMarkerStyle(20); 
							DataPlot->SetMarkerSize(2.); 
							DataPlot->SetLineColor(kBlack);	
							DataPlot->SetMarkerColor(kBlack);
							DataPlot->GetYaxis()->SetRangeUser(-0.02,height*DataPlot->GetMaximum());	
							DataPlot->Draw("e same"); 

//							TH1D* DataPlotG2018 = AcceptanceCorrection(Plots[WhichFSIModel],"hA2018_Final", nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot],xBCut[WhichxBCut]);

//							DataPlotG2018->SetMarkerStyle(20); 
//							DataPlotG2018->SetMarkerSize(2.); 
//							DataPlotG2018->SetLineColor(kRed);	
//							DataPlotG2018->SetMarkerColor(kRed);	
//							DataPlotG2018->Draw("e same"); 

						} else { 

							TH1D* GenieTotalClone = (TH1D*)(Plots[WhichFSIModel]->Clone());
							GenieTotalClone->SetLineColor(kBlack);

							if (FSILabel[WhichFSIModel] == "G2018") { Plots[WhichFSIModel]->SetLineStyle(kDashed); }
							Plots[WhichFSIModel]->Draw("hist C same");							
													
							if (nucleus[WhichNucleus] == "12C" && E[WhichEnergy] == "1_161") {
										
								if (WhichFSIModel == 1)
									{ legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l"); }
							}
							Plots[1]->Draw("hist C same");
							//Plots[0]->Draw("e same"); 

						}

					} // End of the loop over the FSI Models 

					gPad->RedrawAxis();

					// ---------------------------------------------------------------------------------------------------------

				} // End of the loop over the energies

			} // End of the loop over the nuclei

		} // End of the loop over the plots

		// ----------------------------------------------------------------------------------

		//fXSecFile->Close();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TLatex latexCarbon;
		latexCarbon.SetTextFont(FontStyle);
		latexCarbon.SetTextSize(0.79*TextSize);
		PlotCanvas->cd();
		latexCarbon.DrawLatexNDC(0.005,0.77,"^{12}C");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TLatex latexIron;
		latexIron.SetTextFont(FontStyle);
		latexIron.SetTextSize(0.79*TextSize);
		PlotCanvas->cd();
		latexIron.DrawLatexNDC(0.28,0.365,"^{56}Fe");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* pad1GeV = new TPad("pad1GeV","pad1GeV",0.17,0.89,0.33,0.99,21);
		pad1GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad1GeV->Draw();
		pad1GeV->cd();

		TLatex latex1GeV;
		latex1GeV.SetTextFont(FontStyle);
		latex1GeV.SetTextSize(6*TextSize);
		latex1GeV.DrawLatexNDC(0.03,0.45,"1.159 GeV (x1/2)");

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
		latex4GeV.DrawLatexNDC(0.11,0.45,"4.453 GeV (x5)");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* padPmiss = new TPad("padPmiss","padPmiss",0.44,0.03,0.8,0.15,21); 
		padPmiss->SetFillColor(kWhite); 
		PlotCanvas->cd();
		padPmiss->Draw();
		padPmiss->cd();

		TLatex latexPmiss;
		latexPmiss.SetTextFont(FontStyle);
		latexPmiss.SetTextSize(6*TextSize);
		latexPmiss.DrawLatexNDC(0.1,0.5,"(e,e')_{0#pi} E_{QE} [GeV]");

		// -------------------------------------------------------------------------------------------------
	
//		TPad* padx14C12 = new TPad("padx14C12","padx14C12",0.255,0.65,0.305,0.72,21); 
//		padx14C12->SetFillColor(kWhite); 
//		PlotCanvas->cd();
//		padx14C12->Draw();
//		padx14C12->cd();

//		TLatex latexx14C12;
//		latexx14C12.SetTextFont(FontStyle);
//		latexx14C12.SetTextSize(8*TextSize);
//		latexx14C12.DrawLatexNDC(0.,0.4,"x1/2");

		// -------------------------------------------------------------------------------------------------
	
//		TPad* padx5aC12 = new TPad("padx5aC12","padx5aC12",0.855,0.7,0.905,0.77,21); 
//		padx5aC12->SetFillColor(kWhite); 
//		PlotCanvas->cd();
//		padx5aC12->Draw();
//		padx5aC12->cd();

//		TLatex latexx5aC12;
//		latexx5aC12.SetTextFont(FontStyle);
//		latexx5aC12.SetTextSize(8*TextSize);
//		latexx5aC12.DrawLatexNDC(0.,0.4,"x5");

		// -------------------------------------------------------------------------------------------------
	
//		TPad* padx5bC12 = new TPad("padx5bC12","padx5bC12",0.855,0.35,0.905,0.42,21); 
//		padx5bC12->SetFillColor(kWhite); 
//		PlotCanvas->cd();
//		padx5bC12->Draw();
//		padx5bC12->cd();

//		TLatex latexx5bC12;
//		latexx5bC12.SetTextFont(FontStyle);
//		latexx5bC12.SetTextSize(8*TextSize);
//		latexx5bC12.DrawLatexNDC(0.,0.4,"x5");

		// -------------------------------------------------------------------------------------------------
	
//		TPad* padx13C12 = new TPad("padx13C12","padx13C12",0.675,0.91,0.71,0.98,21); 
//		padx13C12->SetFillColor(kWhite); 
//		PlotCanvas->cd();
//		padx13C12->Draw();
//		padx13C12->cd();

//		TLatex latexx13C12;
//		latexx13C12.SetTextFont(FontStyle);
//		latexx13C12.SetTextSize(7.5*TextSize);
//		latexx13C12.DrawLatexNDC(0.,0.4,"x1/3");

		// -------------------------------------------------------------------------------------------------
	
//		TPad* padx13Fe56 = new TPad("padx13Fe56","padx13Fe56",0.675,0.52,0.71,0.59,21); 
//		padx13Fe56->SetFillColor(kWhite); 
//		PlotCanvas->cd();
		//padx13Fe56->Draw();
//		padx13Fe56->cd();

//		TLatex latexx13F56;
//		latexx13F56.SetTextFont(FontStyle);
//		latexx13F56.SetTextSize(7.5*TextSize);
		//latexx13F56.DrawLatexNDC(0.,0.42,"x1/3");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the legend

		PlotCanvas->cd();
		TPad* padLegend = new TPad("padLegend","padLegend",0.06,0.19,0.27,0.53, 21); 
		padLegend->SetFillColor(kWhite); 
		padLegend->Draw();
		padLegend->cd();
		
		//legGenieBreak->AddEntry(Plots[2],"G2018", "l");		

		legGenieBlackLine->SetTextSize(2.*TextSize); 
		legGenieBlackLine->SetBorderSize(0); 
		legGenieBlackLine->Draw();

		legGenieBreak->SetTextSize(2.*TextSize); 
		legGenieBreak->SetBorderSize(0); 
		legGenieBreak->Draw();

		legG2018->AddEntry(Plots[2],"G2018", "l");
		legG2018->SetTextSize(2.*TextSize); 
		legG2018->SetBorderSize(0); 
		legG2018->Draw();		

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Size for y axis label

		double YaxisSize = 7*TextSize; 

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis units carbon

		PlotCanvas->cd();
		TPad* padTitle = new TPad("padTitle","padTitle",0.048,0.58,0.098,1., 21); 
		padTitle->SetFillColor(kWhite); 
		padTitle->Draw();
		padTitle->cd();

		TLatex latexYTitle;
		latexYTitle.SetTextFont(FontStyle);
		latexYTitle.SetTextSize(0.97*YaxisSize);
		latexYTitle.SetTextColor(kBlack);
		latexYTitle.SetTextAngle(90);
		//latexYTitle.DrawLatexNDC(0.8,0.08,DoubleAccCorrXSecTitle);

		if (NameOfPlots[0] == "h_Erec_subtruct_piplpimi_noprot_3pi") 
			{ latexYTitle.DrawLatexNDC(0.6,0.25,"#frac{d#sigma}{dE_{QE}} #left[#frac{#mub}{GeV}#right]"); }
		if (NameOfPlots[0] == "epRecoEnergy_slice_0") 
			{ latexYTitle.DrawLatexNDC(0.62,0.25,"#frac{d#sigma}{dE_{cal}} #left[#frac{#mub}{GeV}#right]"); }

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis units iron

		PlotCanvas->cd();
		TPad* padTitleFe = new TPad("padTitleFe","padTitleFe",0.343,0.18,0.393,0.55,21); 
		padTitleFe->SetFillColor(kWhite); 
		padTitleFe->Draw();
		padTitleFe->cd();

		TLatex latexYTitleFe;
		latexYTitleFe.SetTextFont(FontStyle);
		latexYTitleFe.SetTextSize(0.98*YaxisSize);
		latexYTitleFe.SetTextColor(kBlack);
		latexYTitleFe.SetTextAngle(90);
		//latexYTitleFe.DrawLatexNDC(0.8,0.,DoubleAccCorrXSecTitle);

		if (NameOfPlots[0] == "h_Erec_subtruct_piplpimi_noprot_3pi") 
			{ latexYTitle.DrawLatexNDC(0.6,0.25,"#frac{d#sigma}{dE_{QE}} #left[#frac{#mub}{GeV}#right]"); }
		if (NameOfPlots[0] == "epRecoEnergy_slice_0") 
			{ latexYTitle.DrawLatexNDC(0.6245,0.25,"#frac{d#sigma}{dE_{cal}} #left[#frac{#mub}{GeV}#right]"); }

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString ext = "";
		if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

		PlotCanvas->SaveAs("../../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+ext+"EQE_Fig4_SuSav2_AccCorrXSec.pdf");

		//delete PlotCanvas;

	} // End of the loop over the xB kinematics

} // End of the program

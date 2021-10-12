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

#include "../../myFunctions.cpp"
#include "../../AfroConstants.h"

// ----------------------------------------------------------------------------------------------------------------

void AbsXSec_OverlayEQE_Fig2() {

	// ------------------------------------------------------------------------

	GlobalSettings();
	TGaxis::SetMaxDigits(3);

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut; 
	std::vector<TString> nucleus; 
	std::vector<TString> JustNucleus;
	std::vector<TString> E;
	std::vector<TString> FSIModel;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots; 
	std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames;

	nucleus.push_back("12C"); 
	JustNucleus.push_back("C");
	E.push_back("1_161");

//	nucleus.push_back("12C"); 
//	JustNucleus.push_back("C");
//	E.push_back("2_261");

//	nucleus.push_back("12C"); 
//	JustNucleus.push_back("C");
//	E.push_back("4_461");

//	nucleus.push_back("56Fe"); 
//	JustNucleus.push_back("Fe");
//	E.push_back("2_261");

//	nucleus.push_back("56Fe"); 
//	JustNucleus.push_back("Fe");
//	E.push_back("4_461");

//	nucleus.push_back("4He"); 
//	JustNucleus.push_back("He");
//	E.push_back("2_261");

//	nucleus.push_back("4He"); 
//	JustNucleus.push_back("He");
//	E.push_back("4_461");

	xBCut.push_back("NoxBCut");

	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Pinned Data");

//	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_noprot_3pi"); LabelOfPlots.push_back("(e,e')_{0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("InclusiveeRecoEnergy_slice_0");
	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");

	if (NameOfPlots[0] == "h_Erec_subtruct_piplpimi_noprot_3pi") {

		FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth0pi_WithFidAcc"); FSILabel.push_back("SuSav2");	
		FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithFidAcc"); FSILabel.push_back("G2018");

//		FSIModel.push_back("SuSav2_RadCorr_LFGM_Truth0pi_WithFidAcc_UpdatedSchwinger"); FSILabel.push_back("SuSav2");	
//		FSIModel.push_back("hA2018_Final_RadCorr_LFGM_Truth0pi_WithFidAcc_UpdatedSchwinger"); FSILabel.push_back("G2018");

	} else {

		FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc"); FSILabel.push_back("SuSav2");	
		FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc"); FSILabel.push_back("G2018");

//		FSIModel.push_back("SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger"); FSILabel.push_back("SuSav2");	
//		FSIModel.push_back("hA2018_Final_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger"); FSILabel.push_back("G2018");

	}

//	FSIModel.push_back("SuSav2_RadCorr_LFGM"); FSILabel.push_back("SuSav2");	
//	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("G2018");

	// ------------------------------------------------------------------------

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> BreakDownPlots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				// Loop over the plots				

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					TCanvas* PlotCanvas = new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 205,34,1024,768);

					// ---------------------------------------------------------------------------

					// Dimensions of TPad

					double XMinPadOne = 0., XMaxPadOne = 1., YMinPadOne = 0., YMaxPadOne = 1.;

					TPad* pad1 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadOne,YMinPadOne,XMaxPadOne,YMaxPadOne, 21); 
					pad1->SetFillColor(kWhite); pad1->Draw();
					pad1->SetTopMargin(0.1);
					pad1->SetBottomMargin(0.18);
					pad1->SetLeftMargin(0.17); 
					pad1->SetRightMargin(0.04);
					pad1->cd();

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					double LegXmin = 0.21, LegYmin = 0.45, YSpread = 0.35;

					TLegend* legGenie = new TLegend(LegXmin,LegYmin,LegXmin+0.15,LegYmin+YSpread);
					legGenie->SetNColumns(1);

					TLegend* legGenieBlackLine = new TLegend(LegXmin,0.74,LegXmin+0.15,0.88);
					legGenieBlackLine->SetNColumns(1);

					TLegend* legGenieBreak = new TLegend(0.24,0.63,0.41,0.75);					
					legGenieBreak->SetNColumns(2);
					legGenieBreak->SetMargin(0.35);

					TLegend* legG2018 = new TLegend(0.21,0.57,0.39,0.63);					
					legG2018->SetNColumns(2);
					legG2018->SetMargin(0.39);

					double max = -99.;
					double min = 1E12;

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

						Plots[WhichFSIModel]->GetXaxis()->SetTitle(JustNucleus[WhichNucleus]+LabelOfPlots[WhichPlot]);
						Plots[WhichFSIModel]->GetXaxis()->CenterTitle(0);

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

						// ----------------------------------------------------------------------------------

						// Genie Break Down

						if (
							FSILabel[WhichFSIModel] == "SuSav2"
						) {

							legGenie->AddEntry(Plots[0],"Data", "lep"); 
							legGenieBlackLine->AddEntry(Plots[0],"Data", "lep"); 

							legGenie->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l"); 
//							legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"GENIE (Total)", "l"); 

							BreakDownPlots.clear();

							for (int j = 1; j < 5; j++) {

								if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_3pi") 
									{ BreakDownPlots.push_back( (TH1D*)( FileSample->Get("InclusiveEQE_Int_"+ToStringInt(j)) ) ); }
 
								if (NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0") 
									{ BreakDownPlots.push_back( (TH1D*)( FileSample->Get("ECal_Int_"+ToStringInt(j)) ) ); }

								UniversalE4vFunction(BreakDownPlots[j-1],FSIModelsToLabels[FSIModel[WhichFSIModel]],\
										     nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

								//-----------------------------------------------------------------------------------------------

								BreakDownPlots[j-1]->SetLineColor(BreakDownColors[j-1]);
								BreakDownPlots[j-1]->SetLineWidth(3);

								TLegendEntry* l1 = legGenie->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
								l1->SetTextColor(BreakDownColors[j-1]);

								TLegendEntry* l1Break = legGenieBreak->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
								l1Break->SetTextColor(BreakDownColors[j-1]);

								BreakDownPlots[j-1]->Draw("C hist same");

							} // end of the look over the GENIE break down

						} // End of the SuSav2 if statement for breakdown

						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

						double localmax = Plots[WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }
						double height = 1.1;
						if ( xBCut[WhichxBCut] == "xBCut" ) { height = 1.1; }
						Plots[0]->GetYaxis()->SetRangeUser(0.,height*max);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						// --------------------------------------------------------------------------------------------------

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(2.); 
							Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

							gStyle->SetErrorX(0); // Removing the horizontal errors
							Plots[WhichFSIModel]->Draw("e same"); 

//TFile* f = TFile::Open("12C_1_161_GeV.root","recreate");
//Plots[WhichFSIModel]->Write();
//f->Close();

						} else { 

							if (FSILabel[WhichFSIModel] == "G2018") { Plots[WhichFSIModel]->SetLineStyle(kDashed); }
							Plots[WhichFSIModel]->Draw("C hist same");  // draw them as lines
							
							if (FSILabel[WhichFSIModel] == "SuSav2") 
								{ legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"SuSav2 (Total)", "l"); }

							Plots[0]->Draw("e same"); 

						}

		                                // --------------------------------------------------------------------------------------------------

					} // End of the loop over the FSI Models 

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
					//legGenieBreak->AddEntry(Plots[2],"G2018","l");					

					legG2018->SetBorderSize(0);
					legG2018->SetTextFont(FontStyle);
					legG2018->SetTextSize(TextSize-0.03);
					legG2018->AddEntry(Plots[2],"G2018","l");
					
//					legGenieBlackLine->Draw(); 
//					legGenieBreak->Draw();
//					legG2018->Draw();

					TLatex* myNucleus = new TLatex();
					myNucleus->SetTextFont(FontStyle);
					myNucleus->SetTextColor(kBlack);
					myNucleus->SetTextSize(TextSize-0.02);
					//if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_3pi") { myNucleus->DrawLatexNDC(0.22,0.85,JustNucleus[WhichNucleus]+"(e,e')_{0#pi}"); }
					//if (NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0") { myNucleus->DrawLatexNDC(0.22,0.85,JustNucleus[WhichNucleus]+"(e,e')_{1p0#pi}"); }

					TLatex* myEbeam = new TLatex();
					myEbeam->SetTextFont(FontStyle);
					myEbeam->SetTextColor(kAzure+4);
					myEbeam->SetTextSize(TextSize-0.02);
					if (E[WhichEnergy] == "1_161" && NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_3pi") { myEbeam->DrawLatexNDC(0.72,0.3,"E_{beam}"); }

					TLatex* myArrow = new TLatex();
					myArrow->SetTextFont(FontStyle);
					myArrow->SetTextColor(kAzure+4);
					myArrow->SetTextSize(1.2*TextSize);
					if (E[WhichEnergy] == "1_161" && NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_3pi") { myArrow->DrawLatex(1.141,0.005,"#Downarrow"); }

					// Monitor where 1.161 GeV is
					//TLine* line = new TLine(1.161,0.,1.161,2.);
					//line->Draw();

// -------------------------------------------------------------------------------------------

// Data / SuSav2 ratio

double IntXSecData = IntegratedXSec(Plots[0]);
double IntXSecSuSav2 = IntegratedXSec(Plots[1]);
//cout << "IntXSecData = " << IntXSecData << endl;
//cout << "IntXSecSuSav2 = " << IntXSecSuSav2 << endl;
double DataSuSav2Ratio = IntXSecData / IntXSecSuSav2;

std::cout << "Normalized Yield Data / SuSav2 Ratio = " << DataSuSav2Ratio << endl;

// -------------------------------------------------------------------------------------------

// Extra pad zooming in tail if Ecal plot

if (NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0") {

	// ---------------------------------------

	double PadNDCXmin = 0.35,PadNDCXmax = 0.75, PadNDCYmin = 0.3,PadNDCYmax = 0.65;
	double PadLeftMargin = 0.105, PadRightMargin = 0.05, PadBottomMargin = 0.1;

	double Xmin = 0, Xmax = 0, Ymin = 0, Ymax = 0;

	double LabelOffset = -0.075;

	// ---------------------------------------
		
	if (E[WhichEnergy] == "1_161") { 

		PadNDCXmin = 0.2,PadNDCXmax = 0.77;
		PadNDCYmin = 0.4,PadNDCYmax = 0.8;
		Xmin = 0.565; Xmax = 1.03; 
		if (nucleus[WhichNucleus] == "12C") { Ymin = 0.001; Ymax = 0.055; }

	}

		
	if (E[WhichEnergy] == "2_261") { 

		PadNDCXmin = 0.175,PadNDCXmax = 0.86;
		PadNDCYmin = 0.35,PadNDCYmax = 0.75;
		Xmin = 0.8; Xmax = 2.1; 
		if (nucleus[WhichNucleus] == "4He") { Ymin = 0.0001; Ymax = 0.009; LabelOffset = -0.03; }
		if (nucleus[WhichNucleus] == "12C") { Ymin = 0.001; Ymax = 0.027; }
		if (nucleus[WhichNucleus] == "56Fe") { Ymin = 0.001; Ymax = 0.11; }

	}

	if (E[WhichEnergy] == "4_461") { 

		PadNDCXmin = 0.15,PadNDCXmax = 0.92;
		PadNDCYmin = 0.4,PadNDCYmax = 0.8;
		Xmin = 1.6; Xmax = 4.25; 
		if (nucleus[WhichNucleus] == "4He") { Ymin = 0.0001; Ymax = 0.0035;  LabelOffset = -0.03; }
		if (nucleus[WhichNucleus] == "12C") { Ymin = 0.0001; Ymax = 0.0099; LabelOffset = -0.03; }
		if (nucleus[WhichNucleus] == "56Fe") { Ymin = 0.0001; Ymax = 0.035; }

	}

	// ---------------------------------------

	TH1D* DataPlotClone = (TH1D*)(Plots[0]->Clone("DataPlotClone"));
	TH1D* SuSav2PlotClone = (TH1D*)(Plots[1]->Clone("SuSav2PlotClone"));
	TH1D* G2018PlotClone = (TH1D*)(Plots[2]->Clone("G2018PlotClone"));

	TH1D* QEPlotClone = (TH1D*)(BreakDownPlots[0]->Clone("QEPlotClone"));
	TH1D* MECPlotClone = (TH1D*)(BreakDownPlots[1]->Clone("MECPlotClone"));
	TH1D* RESPlotClone = (TH1D*)(BreakDownPlots[2]->Clone("RESPlotClone"));
	TH1D* DISPlotClone = (TH1D*)(BreakDownPlots[3]->Clone("DISPlotClone"));

	// ---------------------------------------

	double XminLeftLine = DataPlotClone->GetBinCenter(DataPlotClone->FindBin(Xmin));
	double XminRightLine = DataPlotClone->GetBinCenter(DataPlotClone->FindBin(Xmax));

	double XmaxLeftLine = NDCtoX(PadNDCXmin+0.5*PadLeftMargin*(PadNDCXmax-PadNDCXmin));
	double XmaxRightLine = NDCtoX(PadNDCXmax-0.5*PadRightMargin*(PadNDCXmax-PadNDCXmin));

	double YmaxLeftLine = NDCtoY(PadNDCYmin+0.5*PadBottomMargin*(PadNDCYmax-PadNDCYmin));

	double YminLeftLine = DataPlotClone->GetBinContent(DataPlotClone->FindBin(Xmin));
	double YminRightLine = DataPlotClone->GetBinContent(DataPlotClone->FindBin(Xmax));

	TLine* Leftline = new TLine(XminLeftLine,YminLeftLine,XmaxLeftLine,YmaxLeftLine);
	Leftline->SetLineWidth(2);
	Leftline->SetLineColor(kBlack);
	Leftline->SetLineStyle(9);
	//Leftline->Draw();

	TLine* Rightline = new TLine(XminRightLine,YminRightLine,XmaxRightLine,YmaxLeftLine);
	Rightline->SetLineWidth(2);
	Rightline->SetLineColor(kBlack);
	Rightline->SetLineStyle(9);
	//Rightline->Draw();

	// ---------------------------------------

	TString PadName = "ZoomInPad";
	TPad* padZoomIn = new TPad(PadName,PadName,PadNDCXmin,PadNDCYmin,PadNDCXmax,PadNDCYmax,21); 
	padZoomIn->SetFillColor(kWhite); 
	padZoomIn->SetFrameLineWidth(3);
	padZoomIn->Draw();
	padZoomIn->SetTopMargin(0.1);
	padZoomIn->SetBottomMargin(PadBottomMargin);
	padZoomIn->SetLeftMargin(PadLeftMargin);
	padZoomIn->SetRightMargin(PadRightMargin);
	padZoomIn->SetFillStyle(4000); // make pad trasnparent
	padZoomIn->cd();

	// ---------------------------------------------------------------

	auto frame = PlotCanvas->DrawFrame(Xmin,Ymin,Xmax,Ymax);

	frame->GetXaxis()->SetNdivisions(6);
	frame->GetXaxis()->SetLabelSize(0.1);
	frame->GetXaxis()->SetLabelFont(FontStyle);

	frame->GetYaxis()->SetNdivisions(6);
	frame->GetYaxis()->SetLabelSize(0.1);
	frame->GetYaxis()->SetLabelFont(FontStyle);
	frame->GetYaxis()->SetTickLength(-0.02);
	frame->GetYaxis()->SetLabelOffset(LabelOffset);

	DataPlotClone->GetYaxis()->SetNdivisions(4);
	DataPlotClone->GetYaxis()->SetLabelOffset(0.1);
	DataPlotClone->Draw("e same");

	G2018PlotClone->Draw("c hist same");
	SuSav2PlotClone->Draw("c hist same");
	QEPlotClone->Draw("c hist same");
	MECPlotClone->Draw("c hist same");
	RESPlotClone->Draw("c hist same");
	DISPlotClone->Draw("c hist same");

	DataPlotClone->Draw("e same");

	gPad->RedrawAxis();

}

// -------------------------------------------------------------------------------------------

					// -------------------------------------------------------------------------------------------

					TString ext = "";
					if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

					PlotCanvas->SaveAs("../../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"
						+E[WhichEnergy]+"/"+ext+nucleus[WhichNucleus]+"_" 
						+E[WhichEnergy]+"_" +OutputPlotNames[WhichPlot]+"_SuSav2_AbsXSec.pdf");

					//delete PlotCanvas;

					// --------------------------------------------------------------------------------------


				} // End of the loop over the plots
				
				// --------------------------------------------------------------------------------------				

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program

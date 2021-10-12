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

#include "../../myFunctions.cpp"
#include "../../AfroConstants.h"

// -------------------------------------------------------------------------------------------------------------------------------------

double IntegratedXSecError(TH1D* h, int Min = -1, int Max = -1) {

	int NBinsX = h->GetXaxis()->GetNbins();

	int MinBin = 1;
	int MaxBin = NBinsX;

	if ( !(MinBin == -1 && MaxBin == -1) ) {  

		MinBin = Min;
		MaxBin = Max;

	}

	double IntegratedXSecErrorSquared = 0;

	for (int WhichXBin = 1; WhichXBin <= NBinsX; WhichXBin++) {

		if (WhichXBin >= MinBin && WhichXBin <= MaxBin) {

		double BinWidth = h->GetBinWidth(WhichXBin+1);
		double BinError = h->GetBinError(WhichXBin+1);

		IntegratedXSecErrorSquared += TMath::Power(BinError* BinWidth,2.);

		}

	}

	double IntegratedXSecError = TMath::Sqrt(IntegratedXSecErrorSquared);

	return IntegratedXSecError;

}

// ----------------------------------------------------------------------------------------------------------------

void ExtDataTable1() {

	// Range between which we want to know the fraction of reconstructed events

	gStyle->SetOptStat(0);

	// ---------------------------------------------------------------------------------------------------------

	std::vector<TString> xBCut; 
	std::vector<TString> nucleus;
	std::vector<TString> E; 
	std::vector<double> DoubleE;
	std::vector<TString> FSIModel; 
	std::vector<TString> DirNames;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots; 
	std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames;

//	nucleus.push_back("4He");
//	nucleus.push_back("12C");
	nucleus.push_back("56Fe");

//	E.push_back("1_161"); DoubleE.push_back(1.161);
//	E.push_back("2_261"); DoubleE.push_back(2.261);	
	E.push_back("4_461"); DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");

	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Pinned Data"); DirNames.push_back("Pinned Data");
	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("SuSav2");  DirNames.push_back("SuSav2_NoRadCorr");	
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Offset"); FSILabel.push_back("G2018");  DirNames.push_back("hA2018_Truth_RadCorr");

//	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("eRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");

//	NameOfPlots.push_back("h1_Ecal"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_EQE"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");

//	NameOfPlots.push_back("MissMomentum"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} P_{T} [GeV/c]"); OutputPlotNames.push_back("MissMomentum");
//	NameOfPlots.push_back("DeltaAlphaT_Int_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} #delta#alpha_{T} [GeV/c]"); OutputPlotNames.push_back("DeltaAlphaT_Int_0");
//	NameOfPlots.push_back("DeltaPhiT_Int_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} #delta#phi_{T} [GeV/c]"); OutputPlotNames.push_back("DeltaPhiT_Int_0");
	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} Feeddown"); OutputPlotNames.push_back("EcalReso");
//	NameOfPlots.push_back("h1_Ecal_Reso"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} Feeddown"); OutputPlotNames.push_back("ECalReso_FineBin");

//	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_noprot_3pi"); LabelOfPlots.push_back("(e,e')_{0#pi} E_{QE} [GeV]");  OutputPlotNames.push_back("eQE0pi_slice_0");
//	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_noprot_frac_feed3pi"); LabelOfPlots.push_back("(e,e')_{0#pi} E_{QE} Feeddown");  OutputPlotNames.push_back("eQE0piReso_slice_0");

	std::vector<TH1D*> Plots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	TString WhatModelsAreIncluded = "";
	for (int LoopOverFSIModels = 0 ; LoopOverFSIModels < NFSIModels ; LoopOverFSIModels ++) { WhatModelsAreIncluded += "_"+DirNames[LoopOverFSIModels]; };

	std::vector<int> Colors;
	Colors.push_back(kBlack); Colors.push_back(kBlue); Colors.push_back(kRed); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);


	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			double MinE = 1.12;
			double MaxE = 1.2;

			if (DoubleE[WhichEnergy] == 2.261) { MinE = 2.19; MaxE = 2.31; }
			if (DoubleE[WhichEnergy] == 4.461) { MinE = 4.35; MaxE = 4.55; }

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				cout << endl << "---------------------------------------------------------------" << endl;
				//cout << nucleus[WhichNucleus]+"_"+E[WhichEnergy] << "  Emin = " << MinE << "  Emax = " << MaxE << endl;
				cout << nucleus[WhichNucleus]+"_"+E[WhichEnergy] << endl;

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					TString PlotName = nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut];
					TCanvas* PlotCanvas = new TCanvas(PlotName,PlotName,205,34,1024,768);

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					TLegend* leg = new TLegend(0.2,0.5,0.45,0.8);
					leg->SetNColumns(1);

					cout << "---------------------------------------------------------------" << endl << endl;

					// Loop over the FSI Models

					TF1* g[NFSIModels];

					double max = -99.;

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						// ----------------------------------------------------------------------------

						if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_3pi" || NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_frac_feed3pi") { 

							if (FSIModel[WhichFSIModel] == "hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Offset") 
								{ FSIModel[WhichFSIModel] = "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc"; }

							if (FSIModel[WhichFSIModel] == "SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc") 
								{ FSIModel[WhichFSIModel] = "SuSav2_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc"; }

						}

						// ----------------------------------------------------------------------------

						TString PathToFiles = "../../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);
						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						// ----------------------------------------------------------------------------

						PlotCanvas->cd();

						Plots[WhichFSIModel]->GetXaxis()->SetTitle(LabelOfPlots[WhichPlot]);

						UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

						Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
						Plots[WhichFSIModel]->SetLineWidth(3);

						Plots[WhichFSIModel]->GetYaxis()->SetTitle(XSecEcalLabel);

						// ---------------------------------------------------------------------------------------------------------

						// Unfold Data plot

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

						Plots[WhichFSIModel] = AcceptanceCorrection(Plots[WhichFSIModel],"SuSav2", nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot],"NoxBCut");

						}

						// ---------------------------------------------------------------------------------------------------------

						double IntXSec = Plots[WhichFSIModel]->Integral("width"); 
						if (NPlots > 1) { cout << endl << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": integrated xsec = " << IntXSec << endl; }

						// ----------------------------------------------------------------------------

						double ExtraHeight = 0.2;

						double localmax = TMath::Max(max,Plots[WhichFSIModel]->GetMaximum());
						if (localmax > max) { max = localmax; }
						Plots[0]->GetYaxis()->SetRangeUser(0.,(1+ExtraHeight)*max);

						// ---------------------------------------------------------------------------------------------------------

						// Manually define range of integration

						int MinBin = Plots[WhichFSIModel]->FindBin(MinE);
						int MaxBin = Plots[WhichFSIModel]->FindBin(MaxE);

						// Fit a Gaussian around the peak

//						int MinBin = 0;
//						int MaxBin = Plots[WhichFSIModel]->GetXaxis()->GetNbins();

//						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos 
//						|| string(OutputPlotNames[WhichPlot]).find("h1_E") != std::string::npos) {

//							g[WhichFSIModel] = new TF1("g_"+ToString(WhichFSIModel),"gaus",MinE,MaxE);
//							Plots[WhichFSIModel]->Fit(g[WhichFSIModel],"R0Q","",MinE,MaxE); // 0: don't draw it , Q: quiet mode, don't print all this junk
//							g[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);

//							// Percentages with respect to the true beam energy

//							// Using +/- 2 sigma
//							double mean = g[WhichFSIModel]->GetParameter(1);
//							double sigma = g[WhichFSIModel]->GetParameter(2);

//							MinBin = Plots[WhichFSIModel]->FindBin(mean-3*sigma);
//							MaxBin = Plots[WhichFSIModel]->FindBin(mean+3*sigma);

//						}

						// ---------------------------------------------------------------------------------------------------------

						// Energy Reconstruction Percentages Calculation

						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos 
						|| string(OutputPlotNames[WhichPlot]).find("h1_E") != std::string::npos) {

							// ------------------------------------------------------------------------------------------------------------

							// Using hard coded integration limits
//							double IntXSecRange = IntegratedXSec(Plots[WhichFSIModel],MinBin,MaxBin); 

							double IntXSecRange = IntegratedXSec(Plots[WhichFSIModel],MinBin,MaxBin); 

							//cout << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": integrated xsec within peak range = " << IntXSecRange << endl;

							// ------------------------------------------------------------------------------------------------------------

							// Fraction

							int percentage = IntXSecRange / IntXSec * 100.;
							double ErrorOnIntXSecRange = IntegratedXSecError(Plots[WhichFSIModel],MinBin,MaxBin);
							double ErrorOnWholeXSecRange = IntegratedXSecError(Plots[WhichFSIModel]);
							int ErrorOnRatio = percentage * TMath::Sqrt( TMath::Power(ErrorOnIntXSecRange/IntXSecRange,2.) + TMath::Power(ErrorOnWholeXSecRange/IntXSec,2.) );

							if (NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0" && NPlots == 1) { 
		
								cout << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": fraction within range = " << percentage;
								cout << " +/- " << ErrorOnRatio << endl;

							}

							// ------------------------------------------------------------------------------------------------------------

						} // End of the if statement that we are talking about a reco energy plot

						if (string(FSILabel[WhichFSIModel]).find("Data") != std::string::npos) { 

							Plots[WhichFSIModel]->SetMarkerStyle(20);
							Plots[WhichFSIModel]->SetLineColor(kBlack);
							Plots[WhichFSIModel]->SetMarkerColor(kBlack);
							Plots[WhichFSIModel]->Draw("e same"); 
							leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel],"lep"); 

							if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos 
							|| string(OutputPlotNames[WhichPlot]).find("h1_E") != std::string::npos) {

								double sum = Plots[WhichFSIModel]->Integral(MinBin,MaxBin);
								double ErrorOnSum = IntegratedXSecError(Plots[WhichFSIModel],MinBin,MaxBin);

								if (NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0" && NPlots == 1) { 

									cout << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": sum within range = "<< round(sum,3);
									cout << " +/- " << round(ErrorOnSum,3)  << endl << endl;

								}

							}

						}
						else { 

//							Plots[WhichFSIModel]->Draw("c hist same");
							Plots[WhichFSIModel]->Draw("hist same");
							Plots[0]->Draw("e same"); 

							leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel],"l"); 

							if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos 
							|| string(OutputPlotNames[WhichPlot]).find("h1_E") != std::string::npos) {

								double sum = Plots[WhichFSIModel]->Integral(MinBin,MaxBin);
								double ErrorOnSum = IntegratedXSecError(Plots[WhichFSIModel],MinBin,MaxBin);

								if (NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0" && NPlots == 1) { 

									cout << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": sum within range = " << round(sum,3);
									cout << " +/- " << round(ErrorOnSum,3)  << endl << endl;

								}

							}

						}

						// ---------------------------------------------------------------------------------------------------------

						// Energy Reconstruction Percentages Calculation

						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos 
						|| string(OutputPlotNames[WhichPlot]).find("h1_E") != std::string::npos) {

							//g[WhichFSIModel]->Draw("same");

							TLine* line = new TLine(DoubleE[WhichEnergy],0,DoubleE[WhichEnergy],(1+ExtraHeight)*max);
							line->SetLineWidth(3);
							line->SetLineStyle(kDashed);
							line->Draw();

							TLine* linelow = new TLine( MinE,0,MinE,(1+ExtraHeight)*max);
							linelow->SetLineWidth(3);
							linelow->Draw();

							TLine* linehigh = new TLine( MaxE,0,MaxE,(1+ExtraHeight)*max);
							linehigh->SetLineWidth(3);
							linehigh->Draw();

						}

					} // End of the loop over the FSI Models 

					// -----------------------------------------------------------------------------------------------------------------------------------------

					leg->SetBorderSize(0);
					leg->SetTextFont(132);
					leg->SetTextSize(0.07);
					leg->Draw();

					if (NameOfPlots[WhichPlot] != "epRecoEnergy_slice_0" || NPlots > 1) { 

						delete PlotCanvas; 

					}

					cout << "---------------------------------------------------------------" << endl << endl;

					// -----------------------------------------------------------------------------------------------------------------------------------------


				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program

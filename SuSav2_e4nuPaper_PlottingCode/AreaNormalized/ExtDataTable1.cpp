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

//	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
	FSIModel.push_back("Pinned_Data_Final"); FSILabel.push_back("Pinned Data"); DirNames.push_back("Pinned Data");
	
//	FSIModel.push_back("SuSav2_NoRadCorr_LFGM"); FSILabel.push_back("SuSav2");  DirNames.push_back("SuSav2_NoRadCorr");
	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("SuSav2");  DirNames.push_back("SuSav2_NoRadCorr");	
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc"); FSILabel.push_back("G2018");  DirNames.push_back("hA2018_Truth_RadCorr");

//	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("eRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");

//	NameOfPlots.push_back("h1_Ecal"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_EQE"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");

	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_Ecal_SuperFine"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_EQE_SuperFine"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");

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

			double range = 0.06;
			if (DoubleE[WhichEnergy] == 2.261) { range = 0.08; }
			if (DoubleE[WhichEnergy] == 4.461) { range = 0.06; }

			double MinE = (1.-range)*DoubleE[WhichEnergy];
			double MaxE = (1.+range)*DoubleE[WhichEnergy];

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				cout << endl << "---------------------------------------------------------------" << endl;
				cout << nucleus[WhichNucleus]+"_"+E[WhichEnergy] << "  Emin = " << MinE << "  Emax = " << MaxE << endl;

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					TCanvas* PlotCanvas = new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 205,34,1024,768);

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					TLegend* leg = new TLegend(0.2,0.5,0.45,0.8);
					leg->SetNColumns(1);

					cout << "---------------------------------------------------------------" << endl;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);
						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						// ---------------------------------------------------------------------------------------------------------

						// Energy Reconstruction Percentages Calculation

						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos 
						|| string(OutputPlotNames[WhichPlot]).find("h1_E") != std::string::npos) {

							// ----------------------------------------------------------------------------

							// Percentages with respect to the true beam energy

//							double range = 0.06;
//							if (DoubleE[WhichEnergy] == 2.261) { range = 0.08; }
//							if (DoubleE[WhichEnergy] == 4.461) { range = 0.06; }

//							double MinE = (1.-range)*DoubleE[WhichEnergy];
//							double MaxE = (1.+range)*DoubleE[WhichEnergy];

							int MinBin = Plots[WhichFSIModel]->FindBin(MinE);
							int MaxBin = Plots[WhichFSIModel]->FindBin(MaxE);
							int percentage = Plots[WhichFSIModel]->Integral(MinBin,MaxBin) / Plots[WhichFSIModel]->Integral() * 100.;

							// ------------------------------------------------------------------------------------------------------------

							cout << endl << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": fraction within " << int(range*100.) 
							     << "% = " << percentage << endl;
							// ------------------------------------------------------------------------------------------------------------

							PlotCanvas->cd();

							Plots[WhichFSIModel]->GetXaxis()->SetTitle(LabelOfPlots[WhichPlot]);

						UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot]);

							Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
							Plots[WhichFSIModel]->SetLineWidth(3);

							Plots[WhichFSIModel]->GetYaxis()->SetTitle(XSecEcalLabel);

							if ( FSIModel[WhichFSIModel] == "Pinned_Data_Final" ) { 

								Plots[WhichFSIModel]->SetLineColor(kWhite);
								Plots[WhichFSIModel]->SetMarkerColor(kWhite);
								Plots[WhichFSIModel]->Draw("e same"); 

						TH1D* DataPlot = AcceptanceCorrection(Plots[WhichFSIModel],"SuSav2", nucleus[WhichNucleus],E[WhichEnergy],NameOfPlots[WhichPlot],"NoxBCut");

								DataPlot->SetMarkerStyle(20);
								DataPlot->SetMarkerColor(kBlack);
								DataPlot->SetLineColor(kBlack);
								DataPlot->Draw("e same");

								leg->AddEntry(DataPlot,FSILabel[WhichFSIModel],"lep"); 

								double sum = DataPlot->Integral(MinBin,MaxBin);
								cout << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": sum within " << int(range*100.) 
							     << "% = " << sum << endl;

							}
							else { 

								Plots[WhichFSIModel]->Draw("c hist same");

								leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel],"l"); 

								double sum = Plots[WhichFSIModel]->Integral(MinBin,MaxBin);
								cout << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": sum within " << int(range*100.) 
							     << "% = " << sum << endl;

							}

							double ExtraHeight = 0.1;

							TLine* line = new TLine(DoubleE[WhichEnergy],0,DoubleE[WhichEnergy],(1+ExtraHeight)*Plots[WhichFSIModel]->GetMaximum());
							line->SetLineWidth(3);
							line->SetLineStyle(kDashed);
							line->Draw();

							TLine* linelow = new TLine( MinE,0,MinE,(1+ExtraHeight)*Plots[WhichFSIModel]->GetMaximum());
							linelow->SetLineWidth(3);
							linelow->Draw();

							TLine* linehigh = new TLine( MaxE,0,MaxE,(1+ExtraHeight)*Plots[WhichFSIModel]->GetMaximum());
							linehigh->SetLineWidth(3);
							linehigh->Draw();

 							Plots[0]->GetYaxis()->SetRangeUser(0.,(1+ExtraHeight)*Plots[WhichFSIModel]->GetMaximum());
							Plots[0]->Draw("e same");

						} // End of the if statement that we are talking about a reco energy plot

					} // End of the loop over the FSI Models 

					// -----------------------------------------------------------------------------------------------------------------------------------------

//					delete PlotCanvas;

					leg->SetBorderSize(0);
					leg->SetTextFont(132);
					leg->SetTextSize(0.07);
					leg->Draw();

					cout << "---------------------------------------------------------------" << endl << endl;

					// -----------------------------------------------------------------------------------------------------------------------------------------


				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program

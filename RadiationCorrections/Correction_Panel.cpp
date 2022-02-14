#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

using namespace std;

#include "../myFunctions.cpp"
#include "../AfroConstants.h"

double LocalTextSize = 0.09;

// -------------------------------------------------------------------------------------------

TGraph* AccCorrAveragedFunc(std::vector<TH1D*> hVec, TString Energy, TString Var) {

	TH1D::SetDefaultSumw2();

	const int NBins = hVec[0]->GetXaxis()->GetNbins();

	// ---------------------------------------------------------------------------------------------------------

	TString XaxisTitle = hVec[0]->GetXaxis()->GetTitle();  

	// ---------------------------------------------------------------------------------------------------------

	double bincenter[NBins];
	double binentry[NBins];

	for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

		double BinCenter = hVec[0]->GetBinCenter(WhichBin);
		double BinContent = hVec[0]->GetBinContent(WhichBin);

		bincenter[WhichBin-1] = BinCenter;
		binentry[WhichBin-1] = BinContent;

	}

	// ---------------------------------------------------------------------------------------------------------

	TGraph* GraphClone = new TGraph(NBins,bincenter,binentry);
	GraphClone->GetXaxis()->SetTitle(XaxisTitle);

	// ---------------------------------------------------------------------------------------------------------

	return GraphClone;

}

// -------------------------------------------------------------------------------------------

void PrettyPlot(TGraph* h) {

	// ----------------------------------------------------------------------------------------------------------------

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(LocalTextSize);
	h->GetXaxis()->SetTitleSize(0);
	//h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);
	h->GetXaxis()->SetTickLength(0.02);

	// ----------------------------------------------------------------------------------------------------------------

	// Y-axis

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(LocalTextSize); 
	//h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetLabelSize(LocalTextSize);
	h->GetYaxis()->SetTitle("Radiation Correction");
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleOffset(0.9);
	h->GetYaxis()->SetTickLength(0.02);
	//h->GetYaxis()->SetNdivisions(Ndivisions);

}


TGraph* RMSAveragedFunc(TH1D* h, std::vector<TH1D*> hVec, TString Energy, TString Var) {

	TH1D::SetDefaultSumw2();

	double DoubleE = -99., reso = 0.;
	if (Energy == "1_161") { DoubleE = 1.161; reso = 0.07; }
	if (Energy == "2_261") { DoubleE = 2.261; reso = 0.08; }
	if (Energy == "4_461") { DoubleE = 4.461; reso = 0.06; }

	const int NBins = h->GetXaxis()->GetNbins();

	// ---------------------------------------------------------------------------------------------------------

	TH1D* RMSClone = (TH1D*)(hVec[0]->Clone("RMSClone"));
	TString XaxisTitle = RMSClone->GetXaxis()->GetTitle();  
	RMSClone->Add(hVec[1],-1);
	RMSClone->Divide(h);
	RMSClone->Scale(1./TMath::Sqrt(12.));

	double sum = 0; int nbins = 0;

	for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

		double BinCenter = RMSClone->GetBinCenter(WhichBin);
		double BinContent = TMath::Abs(RMSClone->GetBinContent(WhichBin)) * 100.;

		RMSClone->SetBinContent(WhichBin,BinContent);

		if (BinCenter > (1-reso) * DoubleE && BinCenter < (1+reso) * DoubleE ) {

			sum += BinContent; nbins++;

		}

	}

	sum = sum / double(nbins);

	// ---------------------------------------------------------------------------------------------------------

	double bincenter[NBins];
	double binentry[NBins];

	// For Ecal only
	// Use the average around the Ecal peak in the bins (1 +/- reso) Ebeam

	if (Var == "epRecoEnergy_slice_0") {

		for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

			double BinCenter = RMSClone->GetBinCenter(WhichBin);
			double BinContent = RMSClone->GetBinContent(WhichBin);

			if (BinCenter/DoubleE < 1.1) {

				bincenter[WhichBin-1] = BinCenter;

				if (BinCenter > (1-reso) * DoubleE && BinCenter < (1+reso) * DoubleE ) {

					RMSClone->SetBinContent(WhichBin,sum);
					binentry[WhichBin-1] = sum;

				} else { binentry[WhichBin-1] = BinContent; }

			}

		}

	} else {

		for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

			double BinCenter = RMSClone->GetBinCenter(WhichBin);
			double BinContent = RMSClone->GetBinContent(WhichBin);

			bincenter[WhichBin-1] = BinCenter;
			binentry[WhichBin-1] = BinContent;

		}

	}

	// ---------------------------------------------------------------------------------------------------------

	TGraph* GraphClone = new TGraph(NBins,bincenter,binentry);
	GraphClone->GetXaxis()->SetTitle(XaxisTitle);

	// ---------------------------------------------------------------------------------------------------------

	return GraphClone;

}

// -------------------------------------------------------------------------------------------

TGraph* MeanAveragedFunc(std::vector<TH1D*> hVec, TString Energy, TString Var) {

	TH1D::SetDefaultSumw2();

	const int NBins = hVec[0]->GetXaxis()->GetNbins();

	// ---------------------------------------------------------------------------------------------------------

	TString XaxisTitle = hVec[0]->GetXaxis()->GetTitle();  

	// ---------------------------------------------------------------------------------------------------------

	double bincenter[NBins];
	double binentry[NBins];

	for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

		double BinCenter = hVec[0]->GetBinCenter(WhichBin);
		double BinContent = hVec[0]->GetBinContent(WhichBin);

		bincenter[WhichBin-1] = BinCenter;
		binentry[WhichBin-1] = BinContent;

	}

	// ---------------------------------------------------------------------------------------------------------

	TGraph* GraphClone = new TGraph(NBins,bincenter,binentry);
	GraphClone->GetXaxis()->SetTitle(XaxisTitle);

	// ---------------------------------------------------------------------------------------------------------

	return GraphClone;

}

// ----------------------------------------------------------------------------------------------------------------

void Mega_CommonPanelCorrFactors() {

	// ------------------------------------------------------------------------

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	const std::vector<int> Colors{kBlack,610,410,kRed+1,kCyan+1,kBlue,kOrange+1};

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut;
	std::vector<TString> nucleus; 
	std::vector<TString> LabelsOfSamples; 
	std::vector<TString> E;
	std::vector<double> DoubleE; 
	std::vector<TString> FSIModel;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots;
	std::vector<TString> NameOfSubPlots;

	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C");
	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");
	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He");

	E.push_back("1_161"); DoubleE.push_back(1.161);
	E.push_back("2_261"); DoubleE.push_back(2.261);	
	E.push_back("4_461"); DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");

	FSIModel.push_back("SuSav2"); FSILabel.push_back("SuSav2");
//	FSIModel.push_back("hA2018_Final"); FSILabel.push_back("G2018"); // Not to be used for acceptance correction

	TString Var = "epRecoEnergy_slice_0";
//	TString Var = "h_Erec_subtruct_piplpimi_noprot_3pi";
//	TString Var = "MissMomentum";
//	TString Var = "DeltaAlphaT_Int_0";
//	TString Var = "DeltaPhiT_Int_0";

//	NameOfPlots.push_back("Reco_eRecoEnergy_slice_0");
//	NameOfPlots.push_back("TrueWithFid_eRecoEnergy_slice_0");
//	NameOfPlots.push_back("True_eRecoEnergy_slice_0");

//	NameOfPlots.push_back("BkgCorrection_eRecoEnergy_slice_0");
//	NameOfPlots.push_back("FidCorrection_eRecoEnergy_slice_0");
//	NameOfPlots.push_back("AccCorrection_eRecoEnergy_slice_0");

//	NameOfPlots.push_back("Reco_epRecoEnergy_slice_0");
//	NameOfPlots.push_back("TrueWithFid_epRecoEnergy_slice_0");
//	NameOfPlots.push_back("True_epRecoEnergy_slice_0");

//	NameOfPlots.push_back("AccCorrection_epRecoEnergy_slice_0");

//	NameOfPlots.push_back("InverseAccCorrection_"+Var);
	NameOfPlots.push_back("RadiationCorrection_"+Var);

	// ------------------------------------------------------------------------

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> SystPlots;
	std::vector<TH1D*> AccCorrPlots;
	std::vector<TH1D*> average;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	// ------------------------------------------------------------------------

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		TString CanvasName = "CommonPanel_"+xBCut[WhichxBCut];
//		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,2624,768);
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,2500,1500);

		// -------------------------------------------------------------------------------------

		TPad* pad1 = new TPad(E[0],E[0],0.02,0.1,0.36,0.44,21); 
		pad1->SetFillColor(kWhite); pad1->Draw();
		
		TPad* pad2 = new TPad(E[1],E[1],0.36,0.1,0.68,0.44,22); 
		pad2->SetFillColor(kWhite); pad2->Draw(); 
		
		TPad* pad3 = new TPad(E[2],E[2],0.68,0.1,1.,0.44,22); 
		pad3->SetFillColor(kWhite); pad3->Draw(); 	

		pad1->SetTopMargin(0.0);
		pad2->SetTopMargin(0.0);
		pad3->SetTopMargin(0.0);				
		
		pad1->SetBottomMargin(0.18);
		pad2->SetBottomMargin(0.18);
		pad3->SetBottomMargin(0.18);	

		// -------------------------------------------------------------------------------------

		TPad* pad4 = new TPad(E[0],E[0],0.02,0.44,0.36,0.72,21); 
		pad4->SetFillColor(kWhite); pad4->Draw();
		
		TPad* pad5 = new TPad(E[1],E[1],0.36,0.44,0.68,0.72,22); 
		pad5->SetFillColor(kWhite); pad5->Draw(); 
		
		TPad* pad6 = new TPad(E[2],E[2],0.68,0.44,1.,0.72,22); 
		pad6->SetFillColor(kWhite); pad6->Draw(); 	

		pad4->SetTopMargin(0.0);
		pad5->SetTopMargin(0.0);
		pad6->SetTopMargin(0.0);				
		
		pad4->SetBottomMargin(0.0);
		pad5->SetBottomMargin(0.0);
		pad6->SetBottomMargin(0.0);

		// -------------------------------------------------------------------------------------

		TPad* pad7 = new TPad(E[0],E[0],0.02,0.72,0.36,1.,21); 
		pad7->SetFillColor(kWhite); pad7->Draw();
		
		TPad* pad8 = new TPad(E[1],E[1],0.36,0.72,0.68,1.,22); 
		pad8->SetFillColor(kWhite); pad8->Draw(); 
		
		TPad* pad9 = new TPad(E[2],E[2],0.68,0.72,1.,1.,22); 
		pad9->SetFillColor(kWhite); pad9->Draw(); 	

		pad7->SetTopMargin(0.01);
		pad8->SetTopMargin(0.01);
		pad9->SetTopMargin(0.01);				
		
		pad7->SetBottomMargin(0.0);
		pad8->SetBottomMargin(0.0);
		pad9->SetBottomMargin(0.0);

		// -------------------------------------------------------------------------------------

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// ----------------------------------------------------------------------------------------------------------------------------------------------------

			// Radiation Corrections

			if (WhichEnergy == 0) { pad1->cd(); gStyle->SetTitleSize(LocalTextSize,"t"); pad1->SetRightMargin(0.); pad1->SetLeftMargin(0.15); pad1->SetTitle("");}
			if (WhichEnergy == 1)  { pad2->cd(); pad2->SetLeftMargin(0.0); pad2->SetRightMargin(0.0); }
			if (WhichEnergy == 2)  { pad3->cd(); pad3->SetLeftMargin(0.0); pad3->SetRightMargin(0.04); }

			int ColorCounter = 0;

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				if (E[WhichEnergy] == "1_161" && nucleus[WhichNucleus] != "12C") { continue; }

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {	

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					double max = -99.;
					double min = 1E12;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString FileName = "myFiles/Efficiency_"+FSIModel[WhichFSIModel]+"_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+xBCut[WhichxBCut]+".root";
						TFile* FileSample = TFile::Open(FileName);

						Plots.push_back( (TH1D*)( FileSample->Get(FSIModel[WhichFSIModel]+"_"+NameOfPlots[WhichPlot]) ) );

						Plots[WhichFSIModel]->SetLineColor(Colors[ColorCounter]);

					} // End of the loop over the FSI Models 

					std::vector<TH1D*> VectorPlots; VectorPlots.clear();

					for (int i = 0; i < NFSIModels; i++) {

						VectorPlots.push_back(Plots[i]);

					}

					TGraph* clone = MeanAveragedFunc(VectorPlots,E[WhichEnergy],Var);
					ApplyRange(clone,E[WhichEnergy],Var);
				
					// ---------------------------------------------------------------------------------------------------------------------

					PrettyPlot(clone);

					// ---------------------------------------------------------------------------------------------------------------------

					clone->GetYaxis()->SetNdivisions(8);

					clone->GetYaxis()->SetRangeUser(0.5,1.7);
					clone->GetYaxis()->SetLabelSize(LocalTextSize-0.01);
					clone->GetYaxis()->SetTitleSize(LocalTextSize-0.01);
					clone->GetYaxis()->SetTitleOffset(1.);

					// ---------------------------------------------------------------------------------------------------------------------

					clone->SetTitle("");
					clone->SetMarkerStyle(20);
					clone->SetMarkerSize(2.);
					clone->SetMarkerColor(Colors[ColorCounter]);

					if (E[WhichEnergy] == "1_161") { clone->GetXaxis()->SetRangeUser(0.5,1.23); }
					if (E[WhichEnergy] == "2_261") { clone->GetXaxis()->SetRangeUser(0.67,2.4); }
					if (E[WhichEnergy] == "4_461") { clone->GetXaxis()->SetRangeUser(1.5,4.6); }

					if (ColorCounter == 0) { clone->Draw("ap"); }
					else { 

						clone->GetXaxis()->SetLabelSize(LocalTextSize+0.01);
						clone->Draw("p"); 

					}

					// -----------------------------------------------------------------------------------------------------------------------------------------

				} // End of the loop over the plots

				ColorCounter++;

			} // End of the loop over the nuclei

			// -----------------------------------------------------------------------------------------------------------------------------------------

			TLine* line = new TLine(DoubleE[WhichEnergy],0.5,DoubleE[WhichEnergy],1.7);
			line->SetLineColor(kBlack);
			line->SetLineStyle(kDashed);

			line->Draw();

			// ----------------------------------------------------------------------------------------------------------------------------------------------------

			// Systematic Uncertainties

			average.clear();

			if (WhichEnergy == 0) { pad4->cd(); gStyle->SetTitleSize(LocalTextSize,"t"); pad4->SetRightMargin(0.); pad4->SetLeftMargin(0.15); pad4->SetTitle("");}
			if (WhichEnergy == 1)  { pad5->cd(); pad5->SetLeftMargin(0.0); pad5->SetRightMargin(0.0); }
			if (WhichEnergy == 2)  { pad6->cd(); pad6->SetLeftMargin(0.0); pad6->SetRightMargin(0.04); }

			int SystColorCounter = 0;

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				if (E[WhichEnergy] == "1_161" && nucleus[WhichNucleus] != "12C") { continue; }

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {	

					// ---------------------------------------------------------------------------------------

					SystPlots.clear();

					double max = -99.;
					double min = 1E12;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString FileName = "../SystAcceptanceCorrections/myFiles/Efficiency_SuSav2_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+xBCut[WhichxBCut]+".root";
						TFile* FileSample = TFile::Open(FileName);
						SystPlots.push_back( (TH1D*)( FileSample->Get("SuSav2_InverseAccCorrection_"+Var) ) );

						TString FileNameG2018 = "../SystAcceptanceCorrections/myFiles/Efficiency_hA2018_Final_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+xBCut[WhichxBCut]+".root";
						TFile* FileSampleG2018 = TFile::Open(FileNameG2018);
						SystPlots.push_back( (TH1D*)( FileSampleG2018->Get("hA2018_Final_InverseAccCorrection_"+Var) ) );

						SystPlots[0]->SetLineColor(Colors[SystColorCounter]);
						SystPlots[1]->SetLineColor(Colors[SystColorCounter]);

					} // End of the loop over the FSI Models 

					average.push_back( (TH1D*)(SystPlots[0]->Clone()) );
					average[SystColorCounter]->Add(SystPlots[1]);
					average[SystColorCounter]->Scale(0.5);

					std::vector<TH1D*> VectorPlots; VectorPlots.clear();

					VectorPlots.push_back(SystPlots[0]);
					VectorPlots.push_back(SystPlots[1]);

					TGraph* clone = RMSAveragedFunc(average[SystColorCounter],VectorPlots,E[WhichEnergy],Var);
					ApplyRange(clone,E[WhichEnergy],Var);
				
					// ---------------------------------------------------------------------------------------------------------------------

					PrettyPlot(clone);

					// ---------------------------------------------------------------------------------------------------------------------

					clone->GetYaxis()->SetNdivisions(8);

					clone->GetYaxis()->SetRangeUser(-0.5,17);
					clone->GetYaxis()->SetTitle("Acceptance Uncertainty [%]");

					// ---------------------------------------------------------------------------------------------------------------------

					clone->SetTitle("");
					clone->SetMarkerStyle(20);
					clone->SetMarkerSize(2.);
					clone->SetMarkerColor(Colors[SystColorCounter]);

					if (E[WhichEnergy] == "1_161") { clone->GetXaxis()->SetRangeUser(0.5,1.23); }
					if (E[WhichEnergy] == "2_261") { clone->GetXaxis()->SetRangeUser(0.67,2.4); }
					if (E[WhichEnergy] == "4_461") { clone->GetXaxis()->SetRangeUser(1.5,4.6); }

					if (SystColorCounter == 0) { clone->Draw("ap"); }
					else { 

						clone->GetXaxis()->SetLabelSize(LocalTextSize+0.01);
						clone->Draw("p"); 

					}

					// -----------------------------------------------------------------------------------------------------------------------------------------

				} // End of the loop over the plots

				SystColorCounter++;

			} // End of the loop over the nuclei

			// -----------------------------------------------------------------------------------------------------------------------------------------

			line = new TLine(DoubleE[WhichEnergy]+0.005,-0.4,DoubleE[WhichEnergy]+0.005,17.);
			if (WhichEnergy == 1)  { line = new TLine(DoubleE[WhichEnergy]+0.0001,-0.4,DoubleE[WhichEnergy]+0.0001,17.); }
			if (WhichEnergy == 2)  { line = new TLine(DoubleE[WhichEnergy]-0.02,-0.4,DoubleE[WhichEnergy]-0.02,17.); }

			line->SetLineColor(kBlack);
			line->SetLineStyle(kDashed);

			line->Draw();

			// ----------------------------------------------------------------------------------------------------------------------------------------------------

			// Acceptance Corrections

			if (WhichEnergy == 0) { pad7->cd(); gStyle->SetTitleSize(LocalTextSize,"t"); pad7->SetRightMargin(0.); pad7->SetLeftMargin(0.15); pad1->SetTitle("");}
			if (WhichEnergy == 1)  { pad8->cd(); pad8->SetLeftMargin(0.0); pad8->SetRightMargin(0.0); }
			if (WhichEnergy == 2)  { pad9->cd(); pad9->SetLeftMargin(0.0); pad9->SetRightMargin(0.04); }

			int AccCorrColorCounter = 0;

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				if (E[WhichEnergy] == "1_161" && nucleus[WhichNucleus] != "12C") { continue; }

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {	

					// ---------------------------------------------------------------------------------------

					AccCorrPlots.clear();

					double max = -99.;
					double min = 1E12;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString FileName = "../AcceptanceCorrections/myFiles/Efficiency_"+FSIModel[WhichFSIModel]+"_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+xBCut[WhichxBCut]+".root";
						TFile* FileSample = TFile::Open(FileName);

						AccCorrPlots.push_back( (TH1D*)( FileSample->Get("SuSav2_AverageInverseAccCorrection_"+Var) ) );

						AccCorrPlots[WhichFSIModel]->SetLineColor(Colors[AccCorrColorCounter]);

					} // End of the loop over the FSI Models 

					std::vector<TH1D*> VectorPlots; VectorPlots.clear();

					for (int i = 0; i < NFSIModels; i++) {

						VectorPlots.push_back(AccCorrPlots[i]);

					}

					TGraph* clone = AccCorrAveragedFunc(VectorPlots,E[WhichEnergy],Var);
					ApplyRange(clone,E[WhichEnergy],Var);
				
					// ---------------------------------------------------------------------------------------------------------------------

					PrettyPlot(clone);

					// ---------------------------------------------------------------------------------------------------------------------

					clone->GetYaxis()->SetNdivisions(8);

					clone->GetYaxis()->SetRangeUser(0.1,13.9);
					if (Var == "h_Erec_subtruct_piplpimi_noprot_3pi") { clone->GetYaxis()->SetRangeUser(0.1,12); }
					clone->GetYaxis()->SetTitle("Acceptance Correction");

					// ---------------------------------------------------------------------------------------------------------------------

					clone->SetTitle("");
					clone->SetMarkerStyle(20);
					clone->SetMarkerSize(2.);
					clone->SetMarkerColor(Colors[AccCorrColorCounter]);

					if (E[WhichEnergy] == "1_161") { clone->GetXaxis()->SetRangeUser(0.5,1.23); }
					if (E[WhichEnergy] == "2_261") { clone->GetXaxis()->SetRangeUser(0.67,2.4); }
					if (E[WhichEnergy] == "4_461") { clone->GetXaxis()->SetRangeUser(1.5,4.6); }

					if (AccCorrColorCounter == 0) { clone->Draw("ap"); }
					else { 

						clone->GetXaxis()->SetLabelSize(LocalTextSize+0.01);
						clone->Draw("p"); 

					}

					// -----------------------------------------------------------------------------------------------------------------------------------------

				} // End of the loop over the plots

				AccCorrColorCounter++;

			} // End of the loop over the nuclei

			// -----------------------------------------------------------------------------------------------------------------------------------------

			line = new TLine(DoubleE[WhichEnergy],0.1,DoubleE[WhichEnergy],13.9);
			line->SetLineColor(kBlack);
			line->SetLineStyle(kDashed);
			line->Draw();

			// -----------------------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the energies

		// ---------------------------------------

		TLatex* label = new TLatex();
		label->SetTextFont(FontStyle);
		label->SetTextColor(kBlack);
		label->SetTextSize(LocalTextSize);

		pad1->cd();
		label->DrawLatexNDC(0.8,0.9,"(g)");

		label->SetTextSize(LocalTextSize+0.01);
		pad4->cd();
		label->DrawLatexNDC(0.8,0.9,"(d)");

		pad7->cd();
		label->DrawLatexNDC(0.8,0.9,"(a)");

		TLatex *sample12C = new TLatex(); 
		sample12C->SetTextFont(FontStyle); 
		sample12C->SetTextColor(kBlack); 
		sample12C->SetTextSize(LocalTextSize+0.01);
		sample12C->DrawLatexNDC(0.45,0.6,"^{12}C");

		TLatex *E1 = new TLatex(); 
		E1->SetTextFont(FontStyle); 
		E1->SetTextColor(kBlack); 
		E1->SetTextSize(LocalTextSize+0.01);
		E1->DrawLatexNDC(0.23,0.9,"1.159 GeV"); 

		// ---------------------------------------

		label->SetTextSize(LocalTextSize);
		pad2->cd();
		label->DrawLatexNDC(0.8,0.9,"(h)");

		label->SetTextSize(LocalTextSize+0.01);
		pad5->cd();
		label->DrawLatexNDC(0.8,0.9,"(e)");

		pad8->cd();
		label->DrawLatexNDC(0.8,0.9,"(b)");

		TLatex *sample56Fe = new TLatex(); 
		sample56Fe->SetTextFont(FontStyle); 
		sample56Fe->SetTextColor(Colors[1]); 
		sample56Fe->SetTextSize(LocalTextSize+0.01);
		sample56Fe->DrawLatexNDC(0.45,0.19,"^{56}Fe");

		TLatex *sample4He = new TLatex(); 
		sample4He->SetTextFont(FontStyle); 
		sample4He->SetTextColor(Colors[2]); 
		sample4He->SetTextSize(LocalTextSize+0.01);
		sample4He->DrawLatexNDC(0.35,0.48,"^{4}He");

		TLatex *E2 = new TLatex(); 
		E2->SetTextFont(FontStyle); 
		E2->SetTextColor(kBlack); 
		E2->SetTextSize(LocalTextSize+0.01);
		E2->DrawLatexNDC(0.1,0.9,"2.257 GeV"); 

		// ---------------------------------------

		label->SetTextSize(LocalTextSize);
		pad3->cd();
		label->DrawLatexNDC(0.8,0.9,"(i)");

		label->SetTextSize(LocalTextSize+0.01);
		pad6->cd();
		label->DrawLatexNDC(0.8,0.9,"(f)");

		pad9->cd();
		label->DrawLatexNDC(0.8,0.9,"(c)");

		TLatex *E3 = new TLatex(); 
		E3->SetTextFont(FontStyle); 
		E3->SetTextColor(kBlack); 
		E3->SetTextSize(LocalTextSize+0.01);
		E3->DrawLatexNDC(0.1,0.9,"4.453 GeV"); 

		// -----------------------------------------------------------------------------------------------------------------------------------------
		
		PlotCanvas->cd();

		TPad* padEcal = new TPad("Ecal","Ecal",0.,0.03,1.,0.13,21); 
		padEcal->SetFillColor(kWhite); padEcal->Draw(); padEcal->cd();

		TLatex *XLabel = new TLatex(); 
		XLabel->SetTextFont(FontStyle); 
		XLabel->SetTextColor(kBlack); 
		XLabel->SetTextSize(5.*LocalTextSize);
		if (Var == "epRecoEnergy_slice_0") { XLabel->DrawLatexNDC(0.42,0.5,"(e,e'p)_{1p0#pi} E_{cal} [GeV]"); }
		if (Var == "h_Erec_subtruct_piplpimi_noprot_3pi") { XLabel->DrawLatexNDC(0.42,0.5,"(e,e')_{0#pi} E_{QE} [GeV]"); }

		// -----------------------------------------------------------------------------------------------------------------------------------------

		PlotCanvas->SaveAs("MegaCanvas_Corrections.pdf");
		PlotCanvas->SaveAs("MegaCanvas_Corrections.eps");

	} // End of the loop over the xB kinematics

} // End of the program

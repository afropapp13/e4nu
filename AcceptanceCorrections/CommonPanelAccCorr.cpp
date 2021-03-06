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

// -------------------------------------------------------------------------------------------

void PrettyPlot(TGraph* h) {

	// ----------------------------------------------------------------------------------------------------------------

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleSize(0);
	//h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);
	h->GetXaxis()->SetTickLength(0.02);

	// ----------------------------------------------------------------------------------------------------------------

	// Y-axis

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(TextSize); 
	//h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetLabelSize(TextSize);
	h->GetYaxis()->SetTitle("Acceptance Correction");
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleOffset(1.05);
	h->GetYaxis()->SetTickLength(0.02);
	//h->GetYaxis()->SetNdivisions(Ndivisions);

}

// -------------------------------------------------------------------------------------------

TGraph* RMSAveragedFunc(std::vector<TH1D*> hVec, TString Energy, TString Var) {

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

void CommonPanelAccCorr() {

	// ------------------------------------------------------------------------

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	const std::vector<int> Colors{kBlack,610,410,kRed+1,kCyan+1,kBlue,kOrange+1};

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut;
	std::vector<TString> nucleus; 
	std::vector<TString> LabelsOfSamples; 
	std::vector<TString> E;
	std::vector<TString> LabelE; 
	std::vector<TString> FSIModel;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots;
	std::vector<TString> NameOfSubPlots;

	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C");
	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");
	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He");

	E.push_back("1_161"); LabelE.push_back(", 1.1 GeV");
	E.push_back("2_261"); LabelE.push_back(", 2.2 GeV");	
	E.push_back("4_461"); LabelE.push_back(", 4.4 GeV");

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");

	FSIModel.push_back("SuSav2"); FSILabel.push_back("SuSav2");
//	FSIModel.push_back("hA2018_Final"); FSILabel.push_back("G2018"); // Not to be used for acceptance correction

//	TString Var = "epRecoEnergy_slice_0";
	TString Var = "h_Erec_subtruct_piplpimi_noprot_3pi";

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
	NameOfPlots.push_back("AverageInverseAccCorrection_"+Var);

	// ------------------------------------------------------------------------

	std::vector<TH1D*> Plots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	// ------------------------------------------------------------------------

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		TString CanvasName = "CommonPanel_"+xBCut[WhichxBCut];
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,2624,768);

		TPad* pad1 = new TPad(E[0],E[0],0.02,0,0.36,1.,21); 
		pad1->SetFillColor(kWhite); pad1->Draw();
		
		TPad* pad2 = new TPad(E[1],E[1],0.36,0,0.68,1,22); 
		pad2->SetFillColor(kWhite); pad2->Draw(); 
		
		TPad* pad3 = new TPad(E[2],E[2],0.68,0,1.,1,22); 
		pad3->SetFillColor(kWhite); pad3->Draw(); 					
		
		pad1->SetBottomMargin(0.18);
		pad2->SetBottomMargin(0.18);
		pad3->SetBottomMargin(0.18);	

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			if (WhichEnergy == 0) { pad1->cd(); gStyle->SetTitleSize(TextSize,"t"); pad1->SetRightMargin(0.); pad1->SetLeftMargin(0.15); pad1->SetTitle("");}
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

					TGraph* clone = RMSAveragedFunc(VectorPlots,E[WhichEnergy],Var);
					ApplyRange(clone,E[WhichEnergy],Var);
				
					// ---------------------------------------------------------------------------------------------------------------------

					PrettyPlot(clone);

					// ---------------------------------------------------------------------------------------------------------------------

					clone->GetYaxis()->SetNdivisions(8);

					clone->GetYaxis()->SetRangeUser(-0.3,16);
					if (Var == "h_Erec_subtruct_piplpimi_noprot_3pi") { clone->GetYaxis()->SetRangeUser(-0.1,12); }

					// ---------------------------------------------------------------------------------------------------------------------

					clone->SetTitle("");
					clone->SetMarkerStyle(20);
					clone->SetMarkerSize(2.);
					clone->SetMarkerColor(Colors[ColorCounter]);

					if (ColorCounter == 0) { clone->Draw("ap"); }
					else { 

						clone->GetXaxis()->SetLabelSize(TextSize+0.01);
						clone->Draw("p"); 

					}

					// -----------------------------------------------------------------------------------------------------------------------------------------

				} // End of the loop over the plots

				ColorCounter++;

			} // End of the loop over the nuclei

		} // End of the loop over the energies

		// ---------------------------------------

		pad1->cd();
		TLatex *sample12C = new TLatex(); 
		sample12C->SetTextFont(FontStyle); 
		sample12C->SetTextColor(kBlack); 
		sample12C->SetTextSize(TextSize);
		sample12C->DrawLatexNDC(0.45,0.54,"^{12}C");

		TLatex *E1 = new TLatex(); 
		E1->SetTextFont(FontStyle); 
		E1->SetTextColor(kBlack); 
		E1->SetTextSize(TextSize);
		E1->DrawLatexNDC(0.23,0.8,"1.159 GeV"); 

		TLatex* label = new TLatex();
		label->SetTextFont(FontStyle);
		label->SetTextColor(kBlack);
		label->SetTextSize(TextSize);

		label->DrawLatexNDC(0.82,0.8,"(a)");

		// ---------------------------------------

		pad2->cd();
		sample12C->SetTextSize(TextSize+0.006);
		//sample12C->DrawLatexNDC(0.25,0.24,"^{12}C");

		TLatex *sample56Fe = new TLatex(); 
		sample56Fe->SetTextFont(FontStyle); 
		sample56Fe->SetTextColor(Colors[1]); 
		sample56Fe->SetTextSize(TextSize);
		sample56Fe->DrawLatexNDC(0.45,0.34,"^{56}Fe");

		TLatex *sample4He = new TLatex(); 
		sample4He->SetTextFont(FontStyle); 
		sample4He->SetTextColor(Colors[2]); 
		sample4He->SetTextSize(TextSize);
		sample4He->DrawLatexNDC(0.35,0.48,"^{4}He");

		TLatex *E2 = new TLatex(); 
		E2->SetTextFont(FontStyle); 
		E2->SetTextColor(kBlack); 
		E2->SetTextSize(TextSize+0.006);
		E2->DrawLatexNDC(0.1,0.8,"2.257 GeV"); 

		label->SetTextSize(TextSize+0.006);
		label->DrawLatexNDC(0.82,0.8,"(b)");

		// ---------------------------------------

		pad3->cd();
//		sample12C->SetTextSize(TextSize+0.006);
//		sample12C->DrawLatexNDC(0.45,0.2,"^{12}C");

//		sample56Fe->DrawLatexNDC(0.68,0.26,"^{56}Fe");

//		sample4He->DrawLatexNDC(0.5,0.4,"^{4}He");

		TLatex *E3 = new TLatex(); 
		E3->SetTextFont(FontStyle); 
		E3->SetTextColor(kBlack); 
		E3->SetTextSize(TextSize+0.006);
		E3->DrawLatexNDC(0.1,0.8,"4.453 GeV"); 

		label->DrawLatexNDC(0.82,0.8,"(c)");

		// -----------------------------------------------------------------------------------------------------------------------------------------
		
		PlotCanvas->cd();

		TPad* pad4 = new TPad("Ecal","Ecal",0.,0.,1.,0.1,21); 
		pad4->SetFillColor(kWhite); pad4->Draw(); pad4->cd();

		TLatex *XLabel = new TLatex(); 
		XLabel->SetTextFont(FontStyle); 
		XLabel->SetTextColor(kBlack); 
		XLabel->SetTextSize(9.2*TextSize);
		if (Var == "epRecoEnergy_slice_0") { XLabel->DrawLatexNDC(0.4,0.5,"(e,e'p)_{1p0#pi} E_{cal} [GeV]"); }
		if (Var == "h_Erec_subtruct_piplpimi_noprot_3pi") { XLabel->DrawLatexNDC(0.4,0.5,"(e,e')_{0#pi} E_{QE} [GeV]"); }

		// -----------------------------------------------------------------------------------------------------------------------------------------

		if (Var == "epRecoEnergy_slice_0") { PlotCanvas->SaveAs("ECal_PanelAccCorr_"+CanvasName+".pdf"); }
		if (Var == "h_Erec_subtruct_piplpimi_noprot_3pi") { PlotCanvas->SaveAs("EQE_PanelAccCorr_"+CanvasName+".pdf"); }

		// -----------------------------------------------------------------------------------------------------------------------------------------

	} // End of the loop over the xB kinematics

} // End of the program

#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TGaxis.h>
#include <TROOT.h>

using namespace std;

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "../AfroConstants.h"
#include "../myFunctions.cpp"

//TString ToString(double num) {

//	std::ostringstream start;
//	start << num;
//	string start1 = start.str();
//	return start1;
//}

void OverlayXSec() {

	int FontStyle = 132;
	double TextSize = 0.07;
	TGaxis::SetMaxDigits(3);
//	TGaxis::SetExponentOffset(-0.1, 1., "y");
//	TGaxis::SetExponentOffset(-0.05, 0., "y");

	gStyle->SetTitleSize(TextSize-0.02,"t");
	gStyle->SetTitleFont(FontStyle,"t");

	// ---------------------------------------------------------------------------------------------

	TString CanvasName = "DataXSec_Inclusive_Validation";
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				
	PlotCanvas->SetLeftMargin(0.18);
	PlotCanvas->SetBottomMargin(0.17);

	// ---------------------------------------------------------------------------------------------

	// Saclay Data E = 1.299 GeV

	TFile* DataFile = TFile::Open("eQEupdated.root");
	TTree* qent = (TTree*)(DataFile->Get("qent"));
	int n = qent->Draw("v:xsec:xsec_err","Z==6 && A == 12 && E == 1.299 && theta == 37.5","goff");
	TGraphErrors *DataGraph = new TGraphErrors(n,qent->GetV1(),qent->GetV2(),0,qent->GetV3());

	DataGraph->SetTitle("(e,e') ^{12}C, #theta = 37.5^{o}");

	DataGraph->GetXaxis()->SetTitle("Energy Transfer [GeV]");
	DataGraph->GetXaxis()->CenterTitle();
	DataGraph->GetXaxis()->SetNdivisions(8);
	DataGraph->GetXaxis()->SetLabelFont(FontStyle);
	DataGraph->GetXaxis()->SetLabelSize(TextSize);
	DataGraph->GetXaxis()->SetTitleFont(FontStyle);
	DataGraph->GetXaxis()->SetTitleSize(TextSize);
	DataGraph->GetXaxis()->SetTitleOffset(1.1);
	DataGraph->GetXaxis()->SetLabelOffset(0.015);

	DataGraph->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#OmegadE_{e}} [#frac{nb}{sr GeV}]");
	DataGraph->GetYaxis()->CenterTitle();
	DataGraph->GetYaxis()->SetNdivisions(7);
	DataGraph->GetYaxis()->SetLabelFont(FontStyle);
	DataGraph->GetYaxis()->SetLabelSize(TextSize);
	DataGraph->GetYaxis()->SetTitleFont(FontStyle);
	DataGraph->GetYaxis()->SetTitleSize(TextSize);
	DataGraph->GetYaxis()->SetTitleOffset(1.1);

	DataGraph->GetXaxis()->SetRangeUser(0.0,0.85);
	DataGraph->GetYaxis()->SetRangeUser(0,8000); // Y axis in nb / sr / GeV

	DataGraph->SetMarkerStyle(20);
	DataGraph->Draw("AP");

	// ---------------------------------------------------------------------------------------------

	// Saclay Data E = 0.961 GeV

	int nLow = qent->Draw("v:xsec:xsec_err","Z==6 && A == 12 && E == 0.961 && theta == 37.5","goff");
	TGraphErrors *DataGraphLow = new TGraphErrors(nLow,qent->GetV1(),qent->GetV2(),0,qent->GetV3());

	DataGraphLow->SetMarkerStyle(20);
	DataGraphLow->SetMarkerColor(kRed+1);
	DataGraphLow->SetLineColor(kRed+1);
	DataGraphLow->Draw("P");

	// ---------------------------------------------------------------------------------------------

	// SuSav2 GENIE Out-Of-The-Box

	TFile* GenieBoxFile = TFile::Open("GenieOutOfTheBox_12C_DoubleDiff_E_1_161GeV_theta_37_5.root");

	TH1D* GenieBoxPlot = (TH1D*)GenieBoxFile->Get("v_distribution");

	GenieBoxPlot->SetLineColor(610);
	GenieBoxPlot->SetLineWidth(2);

	GenieBoxPlot->Draw("c hist same");

	// ---------------------------------------------------------------------------------------------

	TString PlotName = "h1_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_0";

	// ---------------------------------------------------------------------------------------------

	// SuSav2 e4v

//	TFile* GenieFile = TFile::Open("/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/1_161/SuSav2_RadCorr_LFGM_Truth_WithoutFidAcc_XSec/NoxBCut/12C_1_161_SuSav2_RadCorr_LFGM_Truth_WithoutFidAcc_XSec_Plots_FSI_em.root");

//	TH1D* GeniePlot = (TH1D*)GenieFile->Get(PlotName);
//	UniversalE4vFunction(GeniePlot,"SuSav2","12C","1_161",PlotName);

//	GeniePlot->SetLineColor(kBlue);
//	GeniePlot->SetLineWidth(2);

//	GeniePlot->Scale(1./0.00667); // sr -> solid angle for 36 < theta < 39 & 24 < phi < 36 
//	GeniePlot->Scale(1000); // Conversion from ub to nb
//	GeniePlot->Draw("c hist same");

	// ---------------------------------------------------------------------------------------------

	// Data e4v

	TFile* Datae4vFile = TFile::Open("/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/1_161/Pinned_Data_Final_XSec/NoxBCut/12C_1_161_Pinned_Data_Final_XSec_Plots_FSI_em.root");

	TH1D* Datae4vPlot = (TH1D*)Datae4vFile->Get(PlotName);
	UniversalE4vFunction(Datae4vPlot,"Pinned Data","12C","1_161",PlotName);

	Datae4vPlot->SetLineColor(kGreen+2);
	Datae4vPlot->SetMarkerColor(kGreen+2);
	Datae4vPlot->SetMarkerStyle(20);

	Datae4vPlot->Scale(1./0.0067); // sr: solid angle for 24 < phi < 36 && 36 < theta < 39 
	Datae4vPlot->Scale(1000); // Conversion from ub to nb
	Datae4vPlot->Draw("e1x0 same");


	// ---------------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------

	TLegend* leg = new TLegend(0.35,0.65,0.4,0.89);
	leg->SetNColumns(1);

	leg->AddEntry(DataGraph,"SLAC/Sealock Data E_{e} = 1.299 GeV","lep");

	TLegendEntry* lDataLow = leg->AddEntry(DataGraphLow,"SLAC/Sealock Data E_{e} = 0.961 GeV", "lep");
	lDataLow->SetTextColor(kRed+1);

//	TLegendEntry* lGenie = leg->AddEntry(GeniePlot,"e4v GENIE E_{e} = 1.161 GeV, 24^{o} < #phi_{e'} < 36^{o}", "l");
//	lGenie->SetTextColor(kBlue);

	TLegendEntry* lGenieBox = leg->AddEntry(GenieBoxPlot,"GENIE E_{e} = 1.161 GeV", "l");
	lGenieBox->SetTextColor(610);

	TLegendEntry* lDatae4v = leg->AddEntry(Datae4vPlot,"e4v CLAS Data E_{e} = 1.161 GeV, 24^{o} < #phi_{e'} < 36^{o}", "lep");
	lDatae4v->SetTextColor(kGreen+2);

	leg->SetBorderSize(0);
	leg->SetTextFont(FontStyle);
	leg->SetTextSize(TextSize-0.03);
	leg->Draw();

	// ---------------------------------------------------------------------------------------------

} // End of the program
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TProfile.h>
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

void e4vXSecRatio() {

	TH1D::SetDefaultSumw2();

	int FontStyle = 132;
	double TextSize = 0.07;
	TGaxis::SetMaxDigits(3);
//	TGaxis::SetExponentOffset(-0.1, 1., "y");
//	TGaxis::SetExponentOffset(-0.05, 0., "y");

	gStyle->SetTitleSize(TextSize-0.02,"t");
	gStyle->SetTitleFont(FontStyle,"t");
	gStyle->SetOptStat(0);	

	// ---------------------------------------------------------------------------------------------

	TString CanvasName = "e4vXSec";
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				
	PlotCanvas->SetLeftMargin(0.18);
	PlotCanvas->SetBottomMargin(0.17);

	PlotCanvas->SetGridx();
	PlotCanvas->SetGridy();

	// ---------------------------------------------------------------------------------------------

	TString PlotName = "h1_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_0";

	// ---------------------------------------------------------------------------------------------

	// Data e4v

	TFile* Datae4vFile = TFile::Open("/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/1_161/Pinned_Data_Final_XSec/NoxBCut/12C_1_161_Pinned_Data_Final_XSec_Plots_FSI_em.root");

	TH1D* Datae4vPlot = (TH1D*)Datae4vFile->Get(PlotName);
	UniversalE4vFunction(Datae4vPlot,"Pinned Data","12C","1_161",PlotName);

	Datae4vPlot->SetTitle("(e,e') ^{12}C, #theta = 37.5^{o}");

	Datae4vPlot->GetXaxis()->SetTitle("Energy Transfer [GeV]");
	Datae4vPlot->GetXaxis()->CenterTitle();
	Datae4vPlot->GetXaxis()->SetNdivisions(9);
	Datae4vPlot->GetXaxis()->SetLabelFont(FontStyle);
	Datae4vPlot->GetXaxis()->SetLabelSize(TextSize);
	Datae4vPlot->GetXaxis()->SetTitleFont(FontStyle);
	Datae4vPlot->GetXaxis()->SetTitleSize(TextSize);
	Datae4vPlot->GetXaxis()->SetTitleOffset(1.1);
	Datae4vPlot->GetXaxis()->SetLabelOffset(0.015);

	Datae4vPlot->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#OmegadE_{e}} [#frac{nb}{sr GeV}]");
	Datae4vPlot->GetYaxis()->CenterTitle();
	Datae4vPlot->GetYaxis()->SetNdivisions(7);
	Datae4vPlot->GetYaxis()->SetLabelFont(FontStyle);
	Datae4vPlot->GetYaxis()->SetLabelSize(TextSize);
	Datae4vPlot->GetYaxis()->SetTitleFont(FontStyle);
	Datae4vPlot->GetYaxis()->SetTitleSize(TextSize);
	Datae4vPlot->GetYaxis()->SetTitleOffset(1.1);

	Datae4vPlot->GetXaxis()->SetRangeUser(0.07,0.73);
	Datae4vPlot->GetYaxis()->SetRangeUser(0,4500); // Y axis in nb / sr / GeV

	Datae4vPlot->SetLineColor(kBlue-5);
	Datae4vPlot->SetMarkerColor(kBlue-5);
//	Datae4vPlot->SetMarkerStyle(20);
	Datae4vPlot->SetMarkerStyle(4);
	Datae4vPlot->SetMarkerSize(2.);
	Datae4vPlot->SetLineWidth(2);

	Datae4vPlot->Scale(1./0.0067); // sr: solid angle for 24 < phi < 36 && 36 < theta < 39 
	Datae4vPlot->Scale(1000); // Conversion from ub to nb
	Datae4vPlot->Draw("e1x0 same"); // Include

//	// ---------------------------------------------------------------------------------------------

//	// SuSav2 GENIE Out-Of-The-Box @ E = 1.161 GEV & theta = 37.5 deg

//	TFile* GenieBoxFile = TFile::Open("GenieOutOfTheBox_12C_DoubleDiff_E_1_161GeV_theta_37_5.root");

//	TH1D* GenieBoxPlot = (TH1D*)GenieBoxFile->Get("v_distribution");

//	GenieBoxPlot->SetLineColor(kBlue);
//	GenieBoxPlot->SetLineWidth(2);

//	GenieBoxPlot->Draw("e hist same");

//	TF1 *myfunc = new TF1("myfunc","pol10",0.05,0.8);
//	GenieBoxPlot->Fit(myfunc);
//	myfunc->Draw("same");

	// ---------------------------------------------------------------------------------------------

	// SuSav2 GENIE Out-Of-The-Box @ E = 1.161 GEV & theta = 37.5 deg & fine binning 

	TFile* GenieBoxFileFineBin = TFile::Open("GenieOutOfTheBox_12C_DoubleDiff_E_1_161GeV_theta_37_5_FineBin.root");

	TH1D* GenieBoxPlotFineBin = (TH1D*)GenieBoxFileFineBin->Get("h");

	GenieBoxPlotFineBin->SetLineColor(kBlue-5);
	GenieBoxPlotFineBin->SetLineWidth(2);

//	GenieBoxPlotFineBin->Rebin(); // Include
//	GenieBoxPlotFineBin->Scale(0.5); // Include
//	GenieBoxPlotFineBin->Draw("c hist same"); // Include

	GenieBoxPlotFineBin->Draw("e hist same"); // Include

	TF1 *myfunc = new TF1("myfunc","pol10",0.05,0.75);
// 10th order polynomial using the finely binned plot // Include
	GenieBoxPlotFineBin->Fit(myfunc); // Include
	myfunc->Draw("same"); // Include

	// ---------------------------------------------------------------------------------------------

	// Corrected e4v data to truth (corrected for radiation effects)

	TH1D* AccCorrDatae4vPlot = AcceptanceCorrection(Datae4vPlot,"SuSav2","12C","1_161",PlotName,"NoxBCut","_XSec");
	PlotCanvas->cd();
	AccCorrDatae4vPlot->SetLineColor(kYellow+1);
	AccCorrDatae4vPlot->SetLineWidth(2);
	AccCorrDatae4vPlot->SetMarkerColor(kYellow+1);
	AccCorrDatae4vPlot->SetMarkerStyle(20);

// Cross section corrected for radiation effects
	AccCorrDatae4vPlot->Draw("e1x0 same"); // Include
//	AccCorrDatae4vPlot->Draw("hist p same"); // Include

	// ---------------------------------------------------------------------------------------------

//	// GENIE prediction with veeeeeeeeeeery fine binning

//	TFile* GenieFile = TFile::Open("/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/1_161/SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_XSec/NoxBCut/12C_1_161_SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_XSec_Plots_FSI_em.root");
//	TH1D* GENIEFineBin = (TH1D*)(GenieFile->Get(PlotName)); 

//	GENIEFineBin->SetLineColor(610);
//	GENIEFineBin->SetLineWidth(2);

//	AbsoluteXSecScaling(GENIEFineBin,"SuSav2 NoRad","12C","1_161");
//	ApplyRebinning(GENIEFineBin,"1_161",PlotName);
//	ReweightPlots(GENIEFineBin);

//	GENIEFineBin->Scale(1./0.0067); // sr: solid angle for 24 < phi < 36 && 36 < theta < 39 
//	GENIEFineBin->Scale(1000); // Conversion from ub to nb

//	GENIEFineBin->Draw("c hist same");

	// ---------------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------

	// Data Average vs bin center

	TFile* TProfileFileData = TFile::Open("/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/1_161/Pinned_Data_Final_XSec/NoxBCut/12C_1_161_Pinned_Data_Final_XSec_Plots_FSI_em.root");
	TProfile* profData = (TProfile*)(TProfileFileData->Get("TProf_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_0")); 

	ApplyRebinningTProfile(profData,"1_161",PlotName); // make sure that we have the same binning as in data // Include

	// ---------------------------------------------------------------------------------------------

	TString CanvasNameProfData = "ProfCanvasData";
	TCanvas* PlotCanvasProfData = new TCanvas(CanvasNameProfData,CanvasNameProfData,205,34,1024,768);

	profData->GetXaxis()->SetRangeUser(0.05,0.75);
	profData->GetXaxis()->SetTitle("Energy Transfer [GeV] Bin Center");
	profData->GetYaxis()->SetRangeUser(0.05,0.75);
	profData->GetYaxis()->SetTitle("Energy Transfer [GeV] Average");

	profData->SetMarkerStyle(20);
	profData->SetMarkerSize(2.);
	profData->SetMarkerColor(kBlack);
	profData->Draw("hist p");
	TF1* fData = new TF1("fData","x",0,0.75);
	fData->SetLineColor(kRed);
	fData->Draw("same");

	PlotCanvas->cd();

	// ---------------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------

	// SuSav2 NoRad Average vs bin center

	TFile* TProfileFileSuSav2NoRad = TFile::Open("/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/1_161/SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_XSec/NoxBCut/12C_1_161_SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_XSec_Plots_FSI_em.root");
	TProfile* profSuSav2NoRad = (TProfile*)(TProfileFileSuSav2NoRad->Get("TProf_Omega_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_0")); 

	ApplyRebinningTProfile(profSuSav2NoRad,"1_161",PlotName); // make sure that we have the same binning as in data // Include

	// ---------------------------------------------------------------------------------------------

//	TString CanvasNameProfSuSav2NoRad = "ProfCanvasSuSav2NoRad";
//	TCanvas* PlotCanvasProfSuSav2NoRad = new TCanvas(CanvasNameProfSuSav2NoRad,CanvasNameProfSuSav2NoRad,205,34,1024,768);

//	profSuSav2NoRad->GetXaxis()->SetRangeUser(0.05,0.75);
//	profSuSav2NoRad->GetXaxis()->SetTitle("Energy Transfer [GeV] Bin Center");
//	profSuSav2NoRad->GetYaxis()->SetRangeUser(0.05,0.75);
//	profSuSav2NoRad->GetYaxis()->SetTitle("Energy Transfer [GeV] Average");

	profSuSav2NoRad->SetMarkerStyle(20);
	profSuSav2NoRad->SetMarkerSize(2.);
//	profSuSav2NoRad->SetMarkerColor(kBlack);
profSuSav2NoRad->SetMarkerColor(kBlue);
PlotCanvasProfData->cd();
	profSuSav2NoRad->Draw("hist p same");
//	TF1* fSuSav2NoRad = new TF1("fSuSav2NoRad","x",0,0.75);
//	fSuSav2NoRad->SetLineColor(kRed);
//	fSuSav2NoRad->Draw("same");

	PlotCanvas->cd();

	// ---------------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------

	// Data Average vs bin center for electron theta angle

	TProfile* profDataTheta = (TProfile*)(TProfileFileData->Get("TProf_Theta_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_0")); 

	for (int i = 0; i < 3; i++) { profDataTheta->Rebin(); }

	// ---------------------------------------------------------------------------------------------

	TString CanvasNameProfDataTheta = "ProfCanvasDataTheta";
	TCanvas* PlotCanvasProfDataTheta = new TCanvas(CanvasNameProfDataTheta,CanvasNameProfDataTheta,205,34,1024,768);

	profDataTheta->GetXaxis()->SetTitle("#theta [deg] Bin Center");
	profDataTheta->GetYaxis()->SetTitle("#theta [deg] Average");

	profDataTheta->SetMarkerStyle(20);
	profDataTheta->SetMarkerSize(2.);
	profDataTheta->SetMarkerColor(kBlack);
	profDataTheta->Draw("hist p");
	TF1* fDataTheta = new TF1("fDataTheta","x",36,39);
	fDataTheta->SetLineColor(kRed);
	fDataTheta->Draw("same");

	PlotCanvas->cd();

	// ---------------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------

	// SuSav2 NoRad Average vs bin center

	TProfile* profSuSav2NoRadTheta = (TProfile*)(TProfileFileSuSav2NoRad->Get("TProf_Theta_FullyInclusive_NoQ4Weight_Theta_Slice_InSector_0")); 

	for (int j = 0; j < 3; j++) { profSuSav2NoRadTheta->Rebin(); }

	// ---------------------------------------------------------------------------------------------

//	TString CanvasNameProfSuSav2NoRad = "ProfCanvasSuSav2NoRad";
//	TCanvas* PlotCanvasProfSuSav2NoRad = new TCanvas(CanvasNameProfSuSav2NoRad,CanvasNameProfSuSav2NoRad,205,34,1024,768);

//	profSuSav2NoRad->GetXaxis()->SetRangeUser(0.05,0.75);
//	profSuSav2NoRad->GetXaxis()->SetTitle("Energy Transfer [GeV] Bin Center");
//	profSuSav2NoRad->GetYaxis()->SetRangeUser(0.05,0.75);
//	profSuSav2NoRad->GetYaxis()->SetTitle("Energy Transfer [GeV] Average");

	profSuSav2NoRadTheta->SetMarkerStyle(20);
	profSuSav2NoRadTheta->SetMarkerSize(2.);
//	profSuSav2NoRadTheta->SetMarkerColor(kBlack);
profSuSav2NoRadTheta->SetMarkerColor(kBlue);
PlotCanvasProfDataTheta->cd();
	profSuSav2NoRadTheta->Draw("hist p same");

	PlotCanvas->cd();

	// ---------------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------

	// Bin center corrected xsec

	int NBinsAccCorr = AccCorrDatae4vPlot->GetXaxis()->GetNbins();
	TH1D* BinCentCorrDatae4vPlot = (TH1D*)(AccCorrDatae4vPlot->Clone("BinCentCorrDatae4vPlot"));
	TH1D* BinCenterCorrection = (TH1D*)(AccCorrDatae4vPlot->Clone("BinCenterCorrection"));

	for (int WhichBin = 0; WhichBin < NBinsAccCorr; WhichBin ++) {

//		double BinCenter = prof->GetBinCenter(WhichBin + 1);
//		double BinAverage = prof->GetBinContent(WhichBin + 1);

		double BinCenter = myfunc->Eval(profSuSav2NoRad->GetBinCenter(WhichBin + 1));
		double BinAverage = myfunc->Eval(profSuSav2NoRad->GetBinContent(WhichBin + 1));

		double CorrFactor = BinCenter / BinAverage;

		double CurrentEntry = AccCorrDatae4vPlot->GetBinContent(WhichBin + 1);
		double CurrentError = AccCorrDatae4vPlot->GetBinError(WhichBin + 1);
		double NewEntry = CurrentEntry / CorrFactor;
		double NewError = CurrentError / CorrFactor;

		BinCentCorrDatae4vPlot->SetBinContent(WhichBin+1,NewEntry);
		BinCentCorrDatae4vPlot->SetBinError(WhichBin+1,NewError);

		BinCenterCorrection->SetBinContent(WhichBin+1,1./CorrFactor);
		BinCenterCorrection->SetBinError(WhichBin+1,0.);

	}

	BinCentCorrDatae4vPlot->SetLineColor(kOrange+7);
	BinCentCorrDatae4vPlot->SetLineWidth(2);
	BinCentCorrDatae4vPlot->SetMarkerColor(kOrange+7);

// Final cross section, corrected for radiation & bin centering effects // Include
//	BinCentCorrDatae4vPlot->Draw("e1x0 same");
	BinCentCorrDatae4vPlot->Draw("hist p same");

	// ---------------------------------------------------------------------------------------------

	TString CanvasNameBinCenterCorrection = "BinCenterCorrectionCanvas";
	TCanvas* PlotCanvasBinCenterCorrection = new TCanvas(CanvasNameBinCenterCorrection,CanvasNameBinCenterCorrection,205,34,1024,768);

	// Correction due to bin centering

	BinCenterCorrection->GetXaxis()->SetRangeUser(0.05,0.75);	
	BinCenterCorrection->Draw("e1x0");

	PlotCanvas->cd();	

	// ---------------------------------------------------------------------------------------------
//	// ---------------------------------------------------------------------------------------------

//	TLegend* leg = new TLegend(0.35,0.65,0.4,0.89);
//	leg->SetNColumns(1);

////	TLegendEntry* lGenie = leg->AddEntry(GeniePlot,"e4v GENIE E_{e} = 1.161 GeV, 24^{o} < #phi_{e'} < 36^{o}", "l");
////	lGenie->SetTextColor(kBlue);

//	TLegendEntry* lGenieBox = leg->AddEntry(GenieBoxPlot,"GENIE E_{e} = 1.161 GeV", "l");
//	lGenieBox->SetTextColor(kBlue-5);

//	TLegendEntry* lDatae4v = leg->AddEntry(Datae4vPlot,"e4v CLAS Data E_{e} = 1.161 GeV, 24^{o} < #phi_{e'} < 36^{o}", "lep");
//	lDatae4v->SetTextColor(kBlue-5);

//	leg->SetBorderSize(0);
//	leg->SetTextFont(FontStyle);
//	leg->SetTextSize(TextSize-0.03);
//	leg->Draw();

	// ---------------------------------------------------------------------------------------------

	TLegend* leg = new TLegend(0.3,0.75,0.85,0.89);
	leg->SetNColumns(3);

	TLegendEntry* lDatae4v = leg->AddEntry(Datae4vPlot,"1.161 GeV", "");
	lDatae4v->SetTextColor(kBlue-5);

	leg->SetBorderSize(0);
	leg->SetTextFont(FontStyle);
	leg->SetTextSize(TextSize-0.01);
//	leg->Draw();

	// ---------------------------------------------------------------------------------------------

	TLatex* Latex1161 = new TLatex();
	Latex1161->SetTextFont(FontStyle);
	Latex1161->SetTextSize(TextSize-0.01);
	Latex1161->SetTextColor(kBlue-5);
	Latex1161->SetTextAlign(12);  // centered
	Latex1161->DrawLatexNDC(0.7,0.35,"#splitline{e4#nu 36^{o}-39^{o}}{1.161 GeV}");

	// ---------------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------

	// Ratio Data / GENIE

	// ---------------------------------------------------------------------------------------------

	TString RatioCanvasName = "e4vXSecRatio";
	TCanvas* RatioPlotCanvas = new TCanvas(RatioCanvasName,RatioCanvasName,205,34,1024,768);
				
	RatioPlotCanvas->SetLeftMargin(0.18);
	RatioPlotCanvas->SetBottomMargin(0.17);

	RatioPlotCanvas->SetGridx();
	RatioPlotCanvas->SetGridy();

	int NBins = Datae4vPlot->GetXaxis()->GetNbins();
	TH1D* DataClone = (TH1D*)(Datae4vPlot->Clone("DataClone"));

	for (int WhichBin = 0; WhichBin < NBins; WhichBin++) {

		double DataBinCenter = Datae4vPlot->GetBinCenter(WhichBin+1);
		double DataBinEntry = Datae4vPlot->GetBinContent(WhichBin+1);
		double DataBinError = Datae4vPlot->GetBinError(WhichBin+1);

		int GENIEBin = GenieBoxPlotFineBin->FindBin(DataBinCenter);
		double GENIEBinEntry = GenieBoxPlotFineBin->GetBinContent(GENIEBin);
		double GENIEBinError = GenieBoxPlotFineBin->GetBinError(GENIEBin);

		double NewEntry = 0.;
		double NewError = 0.;

		if (GENIEBinEntry > 0 && DataBinEntry > 0) { 

			NewEntry = DataBinEntry / GENIEBinEntry;
			NewError = NewEntry * TMath::Sqrt( TMath::Power(DataBinError/DataBinEntry,2.) + TMath::Power(GENIEBinError/GENIEBinEntry,2.) ); 

		}

		DataClone->SetBinContent(WhichBin+1,NewEntry);
		DataClone->SetBinError(WhichBin+1,NewError);

	}

	DataClone->GetXaxis()->SetNdivisions(9);
	DataClone->GetYaxis()->SetNdivisions(8);

	DataClone->GetYaxis()->SetRangeUser(0.75,1.55);
	DataClone->GetYaxis()->SetTitle("Data / GENIE");
	DataClone->Draw("e1x0");

//	delete RatioPlotCanvas;

	// ---------------------------------------------------------------------------------------------

//	TString CorrRatioCanvasName = "Corre4vXSecRatio";
//	TCanvas* CorrRatioPlotCanvas = new TCanvas(CorrRatioCanvasName,CorrRatioCanvasName,205,34,1024,768);
//				
//	CorrRatioPlotCanvas->SetLeftMargin(0.18);
//	CorrRatioPlotCanvas->SetBottomMargin(0.17);

//	CorrRatioPlotCanvas->SetGridx();
//	CorrRatioPlotCanvas->SetGridy();

	int CorrNBins = AccCorrDatae4vPlot->GetXaxis()->GetNbins();
	TH1D* CorrDataClone = (TH1D*)(AccCorrDatae4vPlot->Clone("CorrDataClone"));

	for (int WhichBin = 0; WhichBin < CorrNBins; WhichBin++) {

		double DataBinCenter = AccCorrDatae4vPlot->GetBinCenter(WhichBin+1);
		double DataBinEntry = AccCorrDatae4vPlot->GetBinContent(WhichBin+1);
		double DataBinError = AccCorrDatae4vPlot->GetBinError(WhichBin+1);

		int GENIEBin = GenieBoxPlotFineBin->FindBin(DataBinCenter);
		double GENIEBinEntry = GenieBoxPlotFineBin->GetBinContent(GENIEBin);
		double GENIEBinError = GenieBoxPlotFineBin->GetBinError(GENIEBin);

		double NewEntry = 0.;
		double NewError = 0.;

		if (GENIEBinEntry > 0 && DataBinEntry > 0) { 

			NewEntry = DataBinEntry / GENIEBinEntry;
			NewError = NewEntry * TMath::Sqrt( TMath::Power(DataBinError/DataBinEntry,2.) + TMath::Power(GENIEBinError/GENIEBinEntry,2.) ); 

		}

		CorrDataClone->SetBinContent(WhichBin+1,NewEntry);
		CorrDataClone->SetBinError(WhichBin+1,NewError);

	}

//	CorrDataClone->GetXaxis()->SetNdivisions(9);
//	CorrDataClone->GetYaxis()->SetNdivisions(8);

//	CorrDataClone->GetYaxis()->SetRangeUser(0.75,1.45);
//	CorrDataClone->GetYaxis()->SetTitle("Data / GENIE");

// Cross section Data / MC ratio, corrected for radiation // Include
	CorrDataClone->Draw("e1x0 same");
//	CorrDataClone->Draw("hist p same");

	//delete CorrRatioPlotCanvas;

	// ---------------------------------------------------------------------------------------------

//	TString BinCentCorrRatioCanvasName = "BinCentCorre4vXSecRatio";
//	TCanvas* BinCentCorrRatioPlotCanvas = new TCanvas(BinCentCorrRatioCanvasName,BinCentCorrRatioCanvasName,205,34,1024,768);
//				
//	BinCentCorrRatioPlotCanvas->SetLeftMargin(0.18);
//	BinCentCorrRatioPlotCanvas->SetBottomMargin(0.17);

//	BinCentCorrRatioPlotCanvas->SetGridx();
//	BinCentCorrRatioPlotCanvas->SetGridy();

	int BinCentCorrNBins = BinCentCorrDatae4vPlot->GetXaxis()->GetNbins();
	TH1D* BinCentCorrDataClone = (TH1D*)(BinCentCorrDatae4vPlot->Clone("BinCentCorrDataClone"));

	for (int WhichBin = 0; WhichBin < BinCentCorrNBins; WhichBin++) {

		double DataBinCenter = BinCentCorrDatae4vPlot->GetBinCenter(WhichBin+1);
		double DataBinEntry = BinCentCorrDatae4vPlot->GetBinContent(WhichBin+1);
		double DataBinError = BinCentCorrDatae4vPlot->GetBinError(WhichBin+1);

		int GENIEBin = GenieBoxPlotFineBin->FindBin(DataBinCenter);
		double GENIEBinEntry = GenieBoxPlotFineBin->GetBinContent(GENIEBin);
		double GENIEBinError = GenieBoxPlotFineBin->GetBinError(GENIEBin);

		double NewEntry = 0.;
		double NewError = 0.;

		if (GENIEBinEntry > 0 && DataBinEntry > 0) { 

			NewEntry = DataBinEntry / GENIEBinEntry;
			NewError = NewEntry * TMath::Sqrt( TMath::Power(DataBinError/DataBinEntry,2.) + TMath::Power(GENIEBinError/GENIEBinEntry,2.) ); 

		}

		BinCentCorrDataClone->SetBinContent(WhichBin+1,NewEntry);
		BinCentCorrDataClone->SetBinError(WhichBin+1,NewError);

	}

//	CorrDataClone->GetXaxis()->SetNdivisions(9);
//	CorrDataClone->GetYaxis()->SetNdivisions(8);

//	CorrDataClone->GetYaxis()->SetRangeUser(0.75,1.45);
//	CorrDataClone->GetYaxis()->SetTitle("Data / GENIE");

// Final cross section Data / MC ratio, corrected for radiation & bin centering effects // Include
//	BinCentCorrDataClone->Draw("e1x0 same");
	BinCentCorrDataClone->Draw("p hist same");

	//delete CorrRatioPlotCanvas;

	// ---------------------------------------------------------------------------------------------

} // End of the program
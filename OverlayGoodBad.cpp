#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------------------

void PrettyPlot(int font, double size,TH1D* h, int Color) {

	//h->GetXaxis()->SetRangeUser(0.01,0.79);

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetTitleFont(font);
	h->GetXaxis()->SetTitleSize(size);	
	h->GetXaxis()->SetLabelFont(font);
	h->GetXaxis()->SetLabelSize(size);
	h->GetXaxis()->SetTitleOffset(1.1);	
	h->GetXaxis()->SetTitle("p_{n} [GeV/c]");		

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleFont(font);
	h->GetYaxis()->SetTitleSize(size);	
	h->GetYaxis()->SetLabelFont(font);
	h->GetYaxis()->SetLabelSize(size);
	h->GetYaxis()->SetTitleOffset(1.1);		
	h->GetYaxis()->SetTitle("# Events");
	h->GetYaxis()->SetNdivisions(8);		

	h->SetLineColor(Color);
	h->SetMarkerColor(Color);
	h->SetMarkerSize(1.5);
	h->SetMarkerStyle(20);
	h->Rebin();	

	//h->Scale(1./h->GetMaximum());
	//h->GetYaxis()->SetRangeUser(0.01,1.1);
	h->GetYaxis()->SetRangeUser(0.,1.25*h->GetMaximum());		

}

// -------------------------------------------------------------------------------------------------------------------------------------------------

void OverlayGoodBad(TString Id = "") {

	// -----------------------------------------------------------------------------

	int font = 132;
	double size = 0.05;

	// -----------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();	

	// -------------------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetTitleSize(0.0,"t");
	gStyle->SetTitleFont(font,"t");

	gStyle->SetTitleOffset(0.9,"xyz");
	gStyle->SetTitleSize(size,"xyz");
	gStyle->SetTitleFont(font,"xyz");
	gStyle->SetLabelSize(size,"xyz");
	gStyle->SetLabelFont(font,"xyz");
	gStyle->SetNdivisions(10,"xyz");

	gStyle->SetOptStat(0);
	gStyle->SetPadRightMargin(0.13);
	gStyle->SetPadLeftMargin(0.13);	
	gStyle->SetPadBottomMargin(0.13);

	gROOT->ForceStyle();

	TGaxis::SetMaxDigits(3);

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Energies

	std::vector<TString> EnergyTString;	std::vector<double> EnergyDouble; std::vector<TString> EnergyDoubleString;
	EnergyTString.push_back("1161"); EnergyDouble.push_back(1.161); EnergyDoubleString.push_back("1.161");
	EnergyTString.push_back("2261"); EnergyDouble.push_back(2.261); EnergyDoubleString.push_back("2.261");	
	EnergyTString.push_back("4461"); EnergyDouble.push_back(4.461); EnergyDoubleString.push_back("4.461");

	const int NEnergy = EnergyTString.size();

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Nuclei

	std::vector<TString> Nucleus; std::vector<TString> NucleusLatex;
	Nucleus.push_back("C12"); NucleusLatex.push_back("^{12}C");
	Nucleus.push_back("4He"); NucleusLatex.push_back("^{4}He");
	Nucleus.push_back("56Fe"); NucleusLatex.push_back("^{56}Fe");

	const int NNuclei = Nucleus.size();

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Plots & colors to be used for plotting purposes	

	std::vector<TString> PlotName; std::vector<int> Color; std::vector<TString> PlotLatex;
	PlotName.push_back("PMiss_GoodBad"); Color.push_back(kBlue+2); PlotLatex.push_back("P_{Miss}");
	PlotName.push_back("kMiss_GoodBad"); Color.push_back(kOrange+7); PlotLatex.push_back("k_{Miss}");
	PlotName.push_back("PnProxy_GoodBad"); Color.push_back(kGreen+2);PlotLatex.push_back("P_{n,proxy}");		

	const int NPlots = PlotName.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Plots & colors to be used for PL plotting purposes	

	std::vector<TString> PLPlotName; std::vector<int> PLColor; std::vector<TString> PLPlotLatex;
	PLPlotName.push_back("PLFromPMiss_GoodBad"); PLColor.push_back(kOrange+7); PLPlotLatex.push_back("P_{Miss} projection");
	PLPlotName.push_back("PL_GoodBad"); PLColor.push_back(kBlue+2); PLPlotLatex.push_back("P_{L}");

	const int NPLPlots = PLPlotName.size();		

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Slices (ECal reso < 10 %  and > 10%)

	std::vector<TString> Slice; std::vector<TString> SliceLabel; 
	Slice.push_back("_0"); SliceLabel.push_back("E_{Cal}^{Reso} < 10 %");
	Slice.push_back("_1"); SliceLabel.push_back("E_{Cal}^{Reso} > 10 %");

	const int NSlices = Slice.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// TFiles to be opened & plots to be used

	TFile* Files[NEnergy][NNuclei];
	TH1D* Plots[NEnergy][NNuclei][NPlots][NSlices];
	TH2D* PMisskMissPlots[NEnergy][NNuclei][NSlices];
	TH2D* PLFromMissPLPlots[NEnergy][NNuclei][NSlices];	
	TH1D* PMissPlots[NEnergy][NNuclei][NPlots][NSlices];			
	TH1D* PLPlots[NEnergy][NNuclei][NPlots][NSlices];

	// -----------------------------------------------------------------------------------------------------------------------------------------	

	// Loop over the energies
	for (int WhichEnergy = 0; WhichEnergy < NEnergy; WhichEnergy++) {

		// Loop over the nuclei
		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus++) {

			if (EnergyDouble[WhichEnergy] == 1.161 && !(Nucleus[WhichNucleus] == "C12") ) { continue; }

			// -----------------------------------------------------------------------------------------------

			// Open the relevant file
			Files[WhichEnergy][WhichNucleus] = TFile::Open("mySamples/"+Id+"LFNeutrinos_data_e2a_ep_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+"_neutrino6_united4_radphot_test.root","readonly");

			// -----------------------------------------------------------------------------------------------	

			// PT slices (all events, PT < 0.3 GeV/c, PT > 0.3 GeV/c)

			for (int WhichSlice = 0; WhichSlice < NSlices; WhichSlice++) {

				// Create the relevant canvas & Legend

				TString CanvasName= Id+"LFG_Pn_GoodBad_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice];
				TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);	
				can->SetBottomMargin(0.42);

				TLegend* leg = new TLegend(0.7,0.7,0.85,0.89);
				leg->SetTextFont(font);
				leg->SetTextSize(size);
				leg->SetBorderSize(0);	

				TLegend* PLleg = new TLegend(0.2,0.7,0.35,0.85);
				PLleg->SetTextFont(font);
				PLleg->SetTextSize(size);
				PLleg->SetBorderSize(0);				

				TLatex *text = new TLatex();
				text->SetTextFont(font);
				text->SetTextSize(size);

				TLatex *textRegion = new TLatex();
				textRegion->SetTextFont(font);
				textRegion->SetTextSize(size);								

				// -----------------------------------------------------------------------------------------------	

				// Loop over the plots of interest

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice] = (TH1D*)(Files[WhichEnergy][WhichNucleus]->Get(PlotName[WhichPlot]+Slice[WhichSlice]));
					PrettyPlot(font,size,Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],Color[WhichPlot]);
					leg->AddEntry(Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],PlotLatex[WhichPlot],"p");
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("e1x0 same");

					double max = TMath::Max(Plots[WhichEnergy][WhichNucleus][0][WhichSlice]->GetMaximum(),Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->GetMaximum());
					Plots[WhichEnergy][WhichNucleus][0][WhichSlice]->GetYaxis()->SetRangeUser(0.,1.1*max);
					Plots[WhichEnergy][WhichNucleus][0][WhichSlice]->Draw("e1x0 same");					

				}	

				leg->Draw();
				text->DrawLatexNDC(0.4,0.92,NucleusLatex[WhichNucleus] + ", " + EnergyDoubleString[WhichEnergy] + " GeV");
				textRegion->DrawLatexNDC(0.49,0.85,SliceLabel[WhichSlice]);				

				// --------------------------------------------------------------------------------------------------------

				// Residual plots with respect to PMiss

				TPad* pad = new TPad("pad","pad",0.,0.,1.,0.3,21);
				pad->SetFillColor(kWhite);
				pad->SetTopMargin(0.);
				pad->Draw();
				pad->cd();
				pad->SetGridx();
				pad->SetGridy();						

				TH1D* kMissClone = (TH1D*)(Plots[WhichEnergy][WhichNucleus][1][WhichSlice]->Clone());
				kMissClone->Add(Plots[WhichEnergy][WhichNucleus][0][WhichSlice],-1);
				kMissClone->Divide(Plots[WhichEnergy][WhichNucleus][0][WhichSlice]);

				kMissClone->GetXaxis()->SetTitleSize(0.);
				kMissClone->GetXaxis()->SetLabelSize(0.);
				kMissClone->GetXaxis()->SetTickSize(0.1);				

				kMissClone->GetYaxis()->SetTitle("Residual");
				kMissClone->GetYaxis()->SetTitleSize(0.17);
				kMissClone->GetYaxis()->SetLabelSize(0.17);			
				kMissClone->GetYaxis()->SetTitleOffset(0.3);
				kMissClone->GetYaxis()->SetNdivisions(8);

				TH1D* PnProxyClone = (TH1D*)(Plots[WhichEnergy][WhichNucleus][2][WhichSlice]->Clone());	
				PnProxyClone->Add(Plots[WhichEnergy][WhichNucleus][0][WhichSlice],-1);
				PnProxyClone->Divide(Plots[WhichEnergy][WhichNucleus][0][WhichSlice]);			

				double Min = TMath::Min(kMissClone->GetMinimum(),PnProxyClone->GetMinimum());
				double Max = TMath::Max(kMissClone->GetMaximum(),PnProxyClone->GetMaximum());;
				kMissClone->GetYaxis()->SetRangeUser( 1.1*TMath::Max(-0.5,Min),1.1*TMath::Min(0.5,Max));

				kMissClone->GetYaxis()->SetTickSize(0.02);			

				kMissClone->Draw("e1x0 same");					 
				PnProxyClone->Draw("e1x0 same");

				// --------------------------------------------------------------------------------------------------------

				can->SaveAs("myPlots/"+CanvasName+".pdf");			
				delete can;

				// --------------------------------------------------------------------------------------------------------
				// --------------------------------------------------------------------------------------------------------		

				TString CanvasName2D = Id+"PMiss_kMiss_GoodBad_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice];
				TCanvas* can2D = new TCanvas(CanvasName2D,CanvasName2D,205,34,768,768);	
				can2D->SetBottomMargin(0.12);
				can2D->SetRightMargin(0.12);								

				// 2D Pmiss vs kMiss plots

				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice] = (TH2D*)(Files[WhichEnergy][WhichNucleus]->Get("PMiss_kMiss_GoodBad"+Slice[WhichSlice]));

				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->CenterTitle();
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->SetTitle("P_{Miss} [GeV/c]");
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->SetNdivisions(8);				

				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->CenterTitle();
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetTitle("k_{Miss} [GeV/c]");
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetTitleOffset(1.1);				
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetNdivisions(8);								

				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetZaxis()->SetRangeUser(0.,1.1*PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetMaximum());
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->Draw("coltz");
				can2D->SaveAs("myPlots/"+CanvasName2D+".pdf");			
				delete can2D;	

				// --------------------------------------------------------------------------------------------------------
				// --------------------------------------------------------------------------------------------------------		

				TString CanvasName2DPL = Id+"PLFromPMiss_PL_GoodBad_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice];
				TCanvas* can2DPL = new TCanvas(CanvasName2DPL,CanvasName2DPL,205,34,768,768);	
				can2DPL->SetBottomMargin(0.12);
				can2DPL->SetRightMargin(0.13);								

				// 2D Pmiss vs kMiss plots

				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice] = (TH2D*)(Files[WhichEnergy][WhichNucleus]->Get("PLFromPMiss_PL_GoodBad"+Slice[WhichSlice]));

				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->CenterTitle();
				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->SetTitle("P_{L,Miss} [GeV/c]");
				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->SetNdivisions(8);				

				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->CenterTitle();
				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetTitle("P_{L} [GeV/c]");
				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetTitleOffset(1.1);				
				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetNdivisions(8);								

				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetZaxis()->SetRangeUser(0.,1.1*PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetMaximum());
				PLFromMissPLPlots[WhichEnergy][WhichNucleus][WhichSlice]->Draw("coltz");
				can2D->SaveAs("myPlots/"+CanvasName2DPL+".pdf");			
				delete can2D;							

				// -----------------------------------------------------------------------------------------------	

				// Loop over the PL plots of interest

				TString PLCanvasName= Id+"PL_GoodBad_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice];
				TCanvas* PLcan = new TCanvas(PLCanvasName,PLCanvasName,205,34,1024,768);	
				PLcan->SetBottomMargin(0.42);

				for (int WhichPlot = 0; WhichPlot < NPLPlots; WhichPlot++) {

					PLPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice] = (TH1D*)(Files[WhichEnergy][WhichNucleus]->Get(PLPlotName[WhichPlot]+Slice[WhichSlice]));
					PrettyPlot(font,size,PLPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],PLColor[WhichPlot]);
					PLleg->AddEntry(PLPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],PLPlotLatex[WhichPlot],"p");
					PLPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->GetXaxis()->SetTitle("p_{L} [GeV/c]");
					PLPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("e1x0 same");

					double max = TMath::Max(PLPlots[WhichEnergy][WhichNucleus][0][WhichSlice]->GetMaximum(),PLPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->GetMaximum());
					PLPlots[WhichEnergy][WhichNucleus][0][WhichSlice]->GetYaxis()->SetRangeUser(0.,1.1*max);
					PLPlots[WhichEnergy][WhichNucleus][0][WhichSlice]->Draw("e1x0 same");					

				}	

				PLleg->Draw();
				text->DrawLatexNDC(0.4,0.92,NucleusLatex[WhichNucleus] + ", " + EnergyDoubleString[WhichEnergy] + " GeV");
				textRegion->DrawLatexNDC(0.39,0.85,SliceLabel[WhichSlice]);				

				// --------------------------------------------------------------------------------------------------------

				// Residual plots with respect to PMiss

				TPad* PLpad = new TPad("PLpad","PLpad",0.,0.,1.,0.3,21);
				PLpad->SetFillColor(kWhite);
				PLpad->SetTopMargin(0.);
				PLpad->Draw();
				PLpad->cd();
				PLpad->SetGridx();
				PLpad->SetGridy();						

				TH1D* PLFromPMissClone = (TH1D*)(PLPlots[WhichEnergy][WhichNucleus][1][WhichSlice]->Clone());
				PLFromPMissClone->Add(PLPlots[WhichEnergy][WhichNucleus][0][WhichSlice],-1);
				PLFromPMissClone->Divide(PLPlots[WhichEnergy][WhichNucleus][0][WhichSlice]);

				PLFromPMissClone->GetXaxis()->SetTitleSize(0.);
				PLFromPMissClone->GetXaxis()->SetLabelSize(0.);
				PLFromPMissClone->GetXaxis()->SetTickSize(0.1);				

				PLFromPMissClone->GetYaxis()->SetTitle("Residual");
				PLFromPMissClone->GetYaxis()->SetTitleSize(0.17);
				PLFromPMissClone->GetYaxis()->SetLabelSize(0.17);			
				PLFromPMissClone->GetYaxis()->SetTitleOffset(0.3);
				PLFromPMissClone->GetYaxis()->SetNdivisions(8);		

				double PLMin = PLFromPMissClone->GetMinimum();
				double PLMax = PLFromPMissClone->GetMaximum();
				PLFromPMissClone->GetYaxis()->SetRangeUser(1.1*PLMin,1.1*PLMax);

				PLFromPMissClone->GetYaxis()->SetTickSize(0.02);			

				PLFromPMissClone->Draw("e1x0 same");					 
				
				// --------------------------------------------------------------------------------------------------------

				PLcan->SaveAs("myPlots/"+PLCanvasName+".pdf");			
				delete PLcan;
			
			} // End of the loop over the slices 

			// -----------------------------------------------------------------------------------------------			

		} // End of the loop over the nuclei

	} // End of the loop over the energies

	// -----------------------------------------------------------------------------------------------------------------------------------------

} // End of the program 

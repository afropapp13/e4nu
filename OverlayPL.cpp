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

void OverlayPL(TString Id = "") {

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
	PlotName.push_back("PLFromPMiss"); Color.push_back(kBlue+2); PlotLatex.push_back("P_{L,Miss}");
	PlotName.push_back("PL"); Color.push_back(kOrange+7); PlotLatex.push_back("P_{L}");

	const int NPlots = PlotName.size();		

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Slices (all events, PT < 0.3 GeV/c, PT > 0.3 GeV/c)

	std::vector<TString> Slice; std::vector<TString> SliceLabel; 
	Slice.push_back("_0"); SliceLabel.push_back("");
	Slice.push_back("_1"); SliceLabel.push_back("P_{T} < 0.3 GeV/c");	
	Slice.push_back("_2"); SliceLabel.push_back("P_{T} > 0.3 GeV/c");

	const int NSlices = Slice.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// PMiss Slices (e.g. 0 < PMiss < 0.1 GeV/c et al)

	std::vector<TString> PMissSlice; std::vector<TString> PMissSliceLabel; 
	PMissSlice.push_back("_0"); PMissSliceLabel.push_back("P_{L,Miss} < 0.1 GeV");
	PMissSlice.push_back("_1"); PMissSliceLabel.push_back("0.1 < P_{L,Miss} < 0.2 GeV");
	PMissSlice.push_back("_2"); PMissSliceLabel.push_back("0.2 < P_{L,Miss} < 0.3 GeV");
	PMissSlice.push_back("_3"); PMissSliceLabel.push_back("0.3 < P_{L,Miss} < 0.4 GeV");			
	PMissSlice.push_back("_4"); PMissSliceLabel.push_back("0.4 < P_{L,Miss} < 0.5 GeV");
	PMissSlice.push_back("_5"); PMissSliceLabel.push_back("0.5 < P_{L,Miss} < 0.6 GeV");
	PMissSlice.push_back("_6"); PMissSliceLabel.push_back("0.6 < P_{L,Miss} < 0.7 GeV");
	PMissSlice.push_back("_7"); PMissSliceLabel.push_back("0.7 < P_{L,Miss} < 0.8 GeV");
	PMissSlice.push_back("_8"); PMissSliceLabel.push_back("0.8 < P_{L,Miss} < 0.9 GeV");
	PMissSlice.push_back("_9"); PMissSliceLabel.push_back("0.9 < P_{L,Miss} < 1 GeV");

	const int NPMissSlices = PMissSlice.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// TFiles to be opened & plots to be used

	TFile* Files[NEnergy][NNuclei];
	TH1D* Plots[NEnergy][NNuclei][NPlots][NSlices];
	TH2D* PMisskMissPlots[NEnergy][NNuclei][NSlices];
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

				// --------------------------------------------------------------------------------------------------------		

				TString CanvasName2D = Id+"PLFromMiss_PL_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice];
				TCanvas* can2D = new TCanvas(CanvasName2D,CanvasName2D,205,34,768,768);	
				can2D->SetBottomMargin(0.42);				

				// 2D PLFromPMiss vs PL plots

				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice] = (TH2D*)(Files[WhichEnergy][WhichNucleus]->Get("PLFromPMiss_PL"+Slice[WhichSlice]));

				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->CenterTitle();
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->SetTitle("P_{L,Miss} [GeV/c]");
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetXaxis()->SetNdivisions(8);				

				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->CenterTitle();
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetTitle("P_{L} [GeV/c]");
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetTitleOffset(1.1);				
				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->GetYaxis()->SetNdivisions(8);								


				PMisskMissPlots[WhichEnergy][WhichNucleus][WhichSlice]->Draw("coltz");
				can2D->SaveAs("myPlots/"+CanvasName2D+".pdf");			
				delete can2D;								 
			
			} // End of the loop over the slices 

			// -----------------------------------------------------------------------------------------------			

			// PLFromPMiss Slices (10 slices)

			for (int WhichSlice = 0; WhichSlice < NPMissSlices; WhichSlice++) {

				// Create the relevant canvas & Legend

				TString CanvasName= Id+"PL_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+"_Slice"+PMissSlice[WhichSlice];
				TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);	
				can->SetBottomMargin(0.42);

				TLegend* leg = new TLegend(0.7,0.7,0.85,0.89);
				leg->SetTextFont(font);
				leg->SetTextSize(size);
				leg->SetBorderSize(0);	

				TLatex *text = new TLatex();
				text->SetTextFont(font);
				text->SetTextSize(size);

				TLatex *textRegion = new TLatex();
				textRegion->SetTextFont(font);
				textRegion->SetTextSize(size);								

				// -----------------------------------------------------------------------------------------------	

				// Loop over the plots of interest

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

					PMissPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice] = (TH1D*)(Files[WhichEnergy][WhichNucleus]->Get(PlotName[WhichPlot]+"_Slice"+PMissSlice[WhichSlice]));
					PrettyPlot(font,size,PMissPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],Color[WhichPlot]);
					leg->AddEntry(PMissPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],PlotLatex[WhichPlot],"p");
					PMissPlots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("e1x0 same");

				}	

				leg->Draw();
				text->DrawLatexNDC(0.4,0.92,NucleusLatex[WhichNucleus] + ", " + EnergyDoubleString[WhichEnergy] + " GeV");
				textRegion->DrawLatexNDC(0.42,0.85,PMissSliceLabel[WhichSlice]);				

				// --------------------------------------------------------------------------------------------------------

				// Residual plots with respect to PMiss

				TPad* pad = new TPad("pad","pad",0.,0.,1.,0.3,21);
				pad->SetFillColor(kWhite);
				pad->SetTopMargin(0.);
				pad->Draw();
				pad->cd();
				pad->SetGridx();
				pad->SetGridy();						

				TH1D* kMissClone = (TH1D*)(PMissPlots[WhichEnergy][WhichNucleus][1][WhichSlice]->Clone());
				kMissClone->Add(PMissPlots[WhichEnergy][WhichNucleus][0][WhichSlice],-1);
				kMissClone->Divide(PMissPlots[WhichEnergy][WhichNucleus][0][WhichSlice]);

				kMissClone->GetXaxis()->SetTitleSize(0.);
				kMissClone->GetXaxis()->SetLabelSize(0.);
				kMissClone->GetXaxis()->SetTickSize(0.1);				

				kMissClone->GetYaxis()->SetTitle("Residual");
				kMissClone->GetYaxis()->SetTitleSize(0.17);
				kMissClone->GetYaxis()->SetLabelSize(0.17);			
				kMissClone->GetYaxis()->SetTitleOffset(0.3);
				kMissClone->GetYaxis()->SetNdivisions(8);		

				double Min = kMissClone->GetMinimum();
				double Max = kMissClone->GetMaximum();
				kMissClone->GetYaxis()->SetRangeUser( 1.1*TMath::Max(-0.5,Min),1.1*TMath::Min(0.5,Max));

				kMissClone->GetYaxis()->SetTickSize(0.02);			

				kMissClone->Draw("e1x0 hist same");					 

				// --------------------------------------------------------------------------------------------------------

				can->SaveAs("myPlots/"+CanvasName+".pdf");			
				delete can;

			} // End of the PMiss slices

		} // End of the loop over the nuclei

	} // End of the loop over the energies

	// -----------------------------------------------------------------------------------------------------------------------------------------

} // End of the program 

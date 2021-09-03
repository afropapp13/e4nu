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
	h->GetYaxis()->SetRangeUser(0.,1.1*h->GetMaximum());		

}

// -------------------------------------------------------------------------------------------------------------------------------------------------

void OverlayLFGProxy(TString Id = "") {

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
	PlotName.push_back("PMiss"); Color.push_back(kBlue+2); PlotLatex.push_back("P_{Miss}");
	PlotName.push_back("kMiss"); Color.push_back(kOrange+7); PlotLatex.push_back("k_{Miss}");
	PlotName.push_back("PnProxy"); Color.push_back(kGreen+2);PlotLatex.push_back("P_{n,proxy}");		

	const int NPlots = PlotName.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Slices

	std::vector<TString> Slice; std::vector<TString> SliceLabel; 
	Slice.push_back("_0"); SliceLabel.push_back("");

	const int NSlices = Slice.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// TFiles to be opened & plots to be used

	TFile* Files[NEnergy][NNuclei];
	TH1D* Plots[NEnergy][NNuclei][NPlots][NSlices];	

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

			for (int WhichSlice = 0; WhichSlice < NSlices; WhichSlice++) {

				// Create the relevant canvas & Legend

				TString CanvasName= Id+"LFG_Pn_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice];
				TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);	
				can->SetBottomMargin(0.42);

				TLegend* leg = new TLegend(0.7,0.7,0.85,0.89);
				leg->SetTextFont(font);
				leg->SetTextSize(size);
				leg->SetBorderSize(0);	

				TLatex *text = new TLatex();
				text->SetTextFont(font);
				text->SetTextSize(size);				

				// -----------------------------------------------------------------------------------------------	

				// Loop over the plots of interest

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice] = (TH1D*)(Files[WhichEnergy][WhichNucleus]->Get(PlotName[WhichPlot]+Slice[WhichSlice]));
					PrettyPlot(font,size,Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],Color[WhichPlot]);
					leg->AddEntry(Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],PlotLatex[WhichPlot],"p");
					Plots[WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("p0 hist same");

				}	

				leg->Draw();
				text->DrawLatexNDC(0.4,0.92,NucleusLatex[WhichNucleus] + ", " + EnergyDoubleString[WhichEnergy] + " GeV");

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

				kMissClone->Draw("p0 hist same");					 
				PnProxyClone->Draw("p0 hist same");


				// --------------------------------------------------------------------------------------------------------

				can->SaveAs("myPlots/"+CanvasName+".pdf");			
				delete can;

			} // End of the loop over the slices 

		} // End of the loop over the nuclei

	} // End of the loop over the energies

	// -----------------------------------------------------------------------------------------------------------------------------------------

} // End of the program 

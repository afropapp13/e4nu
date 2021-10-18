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

void PrettyPlot(int font, double size,TH1D* h, int Color, TString XLabel, double &SF, int slice) {

	//h->GetXaxis()->SetRangeUser(0.01,0.79);

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetTitleFont(font);
	h->GetXaxis()->SetTitleSize(size);	
	h->GetXaxis()->SetLabelFont(font);
	h->GetXaxis()->SetLabelSize(size);
	h->GetXaxis()->SetTitleOffset(1.1);	
	h->GetXaxis()->SetTitle(XLabel+" [GeV/c]");		

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleFont(font);
	h->GetYaxis()->SetTitleSize(size);	
	h->GetYaxis()->SetLabelFont(font);
	h->GetYaxis()->SetLabelSize(size);
	h->GetYaxis()->SetTitleOffset(1.1);		
	h->GetYaxis()->SetTitle("Peak Normalized To 1");
	h->GetYaxis()->SetNdivisions(8);		

	h->SetLineColor(Color);
	h->SetMarkerColor(Color);
	h->SetMarkerSize(1.5);
	h->SetMarkerStyle(20);
	h->Rebin();	

	double Scale = 1./h->GetMaximum();

	// Use the scaling factors from the 0th slice with all the events
	if (slice == 0) { SF = Scale; }

	h->Scale(SF);
	//h->GetYaxis()->SetRangeUser(0.01,1.1);
	h->GetYaxis()->SetRangeUser(0.,1.2*h->GetMaximum());		

}

// -------------------------------------------------------------------------------------------------------------------------------------------------

void Overlay(TString Id = "") {

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

	// Plots & colors to be used for PL plotting purposes	

	std::vector<TString> PLPlotName; std::vector<int> PLColor; std::vector<TString> PLPlotLatex;
	PLPlotName.push_back("PLFromPMiss"); PLColor.push_back(kOrange+7); PLPlotLatex.push_back("P_{L,Miss}");
	PLPlotName.push_back("PL"); PLColor.push_back(kBlue+2); PLPlotLatex.push_back("P_{L}");

	const int NPLPlots = PLPlotName.size();		

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Slices (all events, PT < 0.3 GeV/c, PT > 0.3 GeV/c)

	std::vector<TString> Slice; std::vector<TString> SliceLabel; 
	Slice.push_back("_0"); SliceLabel.push_back("All events");
	Slice.push_back("_GoodBad_0"); SliceLabel.push_back("E_{Cal}^{Reso} < 10 %");	
	Slice.push_back("_GoodBad_1"); SliceLabel.push_back("E_{Cal}^{Reso} > 10 %");

	const int NSlices = Slice.size();	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// TFiles to be opened & plots to be used

	const int NSamples = 2; // 0 = Data, 1 = SuSav2 Rad MC

	TFile* Files[NSamples][NEnergy][NNuclei];

	TH1D* Plots[NSamples][NEnergy][NNuclei][NPlots][NSlices];		
	TH1D* PLPlots[NSamples][NEnergy][NNuclei][NPlots][NSlices];

	double SF[NSamples][NEnergy][NNuclei][NPlots];		
	double PLSF[NSamples][NEnergy][NNuclei][NPlots];

	// -----------------------------------------------------------------------------------------------------------------------------------------	

	// Loop over the energies

	for (int WhichEnergy = 0; WhichEnergy < NEnergy; WhichEnergy++) {

		// Loop over the nuclei

		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus++) {

			if (EnergyDouble[WhichEnergy] == 1.161 && !(Nucleus[WhichNucleus] == "C12") ) { continue; }

			// -----------------------------------------------------------------------------------------------

			// Open the relevant files

			Files[0][WhichEnergy][WhichNucleus] = TFile::Open("mySamples/"+Id+"LFNeutrinos_genie_e2a_ep_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+"_neutrino6_united4_radphot_test_SuSav2_Rad.root","readonly");
			Files[1][WhichEnergy][WhichNucleus] = TFile::Open("mySamples/"+Id+"LFNeutrinos_data_e2a_ep_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+"_neutrino6_united4_radphot_test.root","readonly");

			// -----------------------------------------------------------------------------------------------	

			// All events, EcalReso < 10%, EcalReso > 10% 

			for (int WhichSlice = 0; WhichSlice < NSlices; WhichSlice++) {	

				// PMiss, kMiss, PnProxy

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

					// Create the relevant canvas & Legend

					TString CanvasName= Id+"Overlay_LFG_Pn_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice]+"_"+PlotName[WhichPlot];
					TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);	
					can->SetBottomMargin(0.12);

					TLegend* leg = new TLegend(0.65,0.7,0.85,0.89);
					leg->SetTextFont(font);
					leg->SetTextSize(size);
					leg->SetBorderSize(0);
					leg->SetNColumns(2);															

					// -----------------------------------------------------------------------------------------------

					// Loop over data/MC

					for (int WhichSample = 0; WhichSample < NSamples; WhichSample++) {					

						Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice] = (TH1D*)(Files[WhichSample][WhichEnergy][WhichNucleus]->Get(PlotName[WhichPlot]+Slice[WhichSlice]));
						PrettyPlot(font,size,Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],Color[WhichPlot],PlotLatex[WhichPlot],SF[WhichSample][WhichEnergy][WhichNucleus][WhichPlot],WhichSlice);
						Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->SetTitle(NucleusLatex[WhichNucleus] + ", " + EnergyDoubleString[WhichEnergy] + " GeV, " + SliceLabel[WhichSlice]);

						if (WhichSample == 0) {
							leg->AddEntry(Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],"Data","p");
							Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("e1x0 same");
						}

						if (WhichSample == 1) {
							leg->AddEntry(Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],"MC","l");

							double max = TMath::Max(Plots[0][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->GetMaximum(),Plots[1][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->GetMaximum());

							Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("c hist same");
							Plots[0][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->GetYaxis()->SetRangeUser(0,1.1*max);
							Plots[0][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("e1x0 same");


						}						

					}	

					leg->Draw();					

					can->SaveAs("myPlots/"+CanvasName+".pdf");			
					delete can;							

				}	
			
				// -----------------------------------------------------------------------------------------------

				// PL, PLFromPMiss

				for (int WhichPlot = 0; WhichPlot < NPLPlots; WhichPlot++) {

					// Create the relevant canvas & Legend

					TString CanvasName= Id+"Overlay_PL_"+Nucleus[WhichNucleus]+"_"+EnergyTString[WhichEnergy]+Slice[WhichSlice]+"_"+PLPlotName[WhichPlot];
					TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);	
					can->SetBottomMargin(0.12);

					TLegend* leg = new TLegend(0.25,0.7,0.45,0.89);
					leg->SetTextFont(font);
					leg->SetTextSize(size);
					leg->SetBorderSize(0);	
					leg->SetNColumns(2);														

					// -----------------------------------------------------------------------------------------------

					// Loop over data/MC

					for (int WhichSample = 0; WhichSample < NSamples; WhichSample++) {					

						Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice] = (TH1D*)(Files[WhichSample][WhichEnergy][WhichNucleus]->Get(PLPlotName[WhichPlot]+Slice[WhichSlice]));
						PrettyPlot(font,size,Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],PLColor[WhichPlot],PLPlotLatex[WhichPlot],PLSF[WhichSample][WhichEnergy][WhichNucleus][WhichPlot],WhichSlice);
						Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->SetTitle(NucleusLatex[WhichNucleus] + ", " + EnergyDoubleString[WhichEnergy] + " GeV, " + SliceLabel[WhichSlice]);

						if (WhichSample == 0) {
							leg->AddEntry(Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],"Data","p");
							Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("e1x0 same");
						}

						if (WhichSample == 1) {
							leg->AddEntry(Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice],"MC","l");
							Plots[WhichSample][WhichEnergy][WhichNucleus][WhichPlot][WhichSlice]->Draw("c hist same");
						}						

					}	

					leg->Draw();					

					can->SaveAs("myPlots/"+CanvasName+".pdf");			
					delete can;							

				}	
			
			} // End of the loop over the slices

			// -----------------------------------------------------------------------------------------------			

		} // End of the loop over the nuclei

	} // End of the loop over the energies

	// -----------------------------------------------------------------------------------------------------------------------------------------

} // End of the program 

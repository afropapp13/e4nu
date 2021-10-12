#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>

using namespace std;

void OverlayECal_FigExtData11() {

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	TString TxtName = "../AccCorrXSec/myXSec/UndetPionSub_XSecs.txt";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);
	myTxtFile << std::fixed << std::setprecision(3);
	myTxtFile << "Undetected Pion Subtraction Ext Data Fig.7;" << endl << endl;
	myTxtFile << "Bin #;Low bin edge; High bin edge;Weighted Events/GeV;Uncertainty" << endl << endl;

	// ---------------------------------------------------------------------------------------------------------------------------------------------
 
	TString target = "56Fe";

	double xmin = 1.8;
	double xmax = 5;
	int N_E_bins;

	TString version = "v3_0_6/";
	double SizeMarker = 2.;

	int FontStyle = 132;
	double TextSize = 0.07;
	int Ndivisions = 5;

	TLatex *lat1 = new TLatex();
	TLatex *lat2 = new TLatex();

	lat1->SetNDC();
	lat2->SetNDC();

	TString PathToFile = "/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/4_461/Pinned_Data_Final/NoxBCut/";
	TFile* file_in = new TFile(PathToFile+target+"_4_461_Pinned_Data_Final_Plots_FSI_em.root");

	gStyle->SetOptStat(0);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetNdivisions(Ndivisions,"Y");
	gStyle->SetNdivisions(Ndivisions,"X");
	gStyle->SetHistLineWidth(5);
	TGaxis::SetMaxDigits(3);

	//------------- e- only-------------------------------------------------------------

	TH1D *h1_Erec_0pi_e = (TH1D*)file_in->Get("h1_E_rec_0pi");
	TH1D *h1_E_rec = (TH1D*)file_in->Get("h1_E_rec");
	TH1D *h1_Erec_0pi_sub_e = (TH1D*)file_in->Get("h_Erec_subtruct_piplpimi_noprot_4pi");
	TH1D *h1_Erec_1pi_e_weight = (TH1D*)file_in->Get("h1_E_rec_1pi_weight");
	TH1D *h1_Erec_2pi_e_weight = (TH1D*)file_in->Get("h1_E_rec_2pi_weight");
	TH1D *h1_Erec_3pi_e_weight = (TH1D*)file_in->Get("h1_E_rec_3pi_weight");
	TH1D *h1_Erec_4pi_e_weight = (TH1D*)file_in->Get("h1_E_rec_4pi_weight");
	TH1D *h1_Erec_1pi_e = (TH1D*)file_in->Get("h1_E_rec_1pi");

	h1_Erec_0pi_e->SetLineColor(38);
	h1_Erec_0pi_sub_e->SetLineColor(38);
	h1_Erec_1pi_e_weight->SetLineColor(38);
	h1_Erec_2pi_e_weight->SetLineColor(38);
	h1_Erec_3pi_e_weight->SetLineColor(38);
	h1_Erec_4pi_e_weight->SetLineColor(38);
	h1_Erec_1pi_e->SetLineColor(38);
	h1_E_rec->SetLineColor(38); 

	N_E_bins=h1_Erec_0pi_e->GetNbinsX();
//	double Scale = 1E-6;

	// ---------------------------------------------------------------------------------------------------- 

	// Account for the different bin widths

	for(int i=1;i<=N_E_bins;i++) {

		h1_Erec_0pi_e->SetBinContent(i,h1_Erec_0pi_e->GetBinContent(i)/h1_Erec_0pi_e->GetBinWidth(i));
		h1_Erec_0pi_sub_e->SetBinContent(i,h1_Erec_0pi_sub_e->GetBinContent(i)/h1_Erec_0pi_sub_e->GetBinWidth(i));
		h1_Erec_1pi_e_weight->SetBinContent(i,h1_Erec_1pi_e_weight->GetBinContent(i)/h1_Erec_1pi_e_weight->GetBinWidth(i));
		h1_Erec_2pi_e_weight->SetBinContent(i,h1_Erec_2pi_e_weight->GetBinContent(i)/h1_Erec_2pi_e_weight->GetBinWidth(i));
		h1_Erec_3pi_e_weight->SetBinContent(i,h1_Erec_3pi_e_weight->GetBinContent(i)/h1_Erec_3pi_e_weight->GetBinWidth(i));
		h1_Erec_4pi_e_weight->SetBinContent(i,h1_Erec_4pi_e_weight->GetBinContent(i)/h1_Erec_4pi_e_weight->GetBinWidth(i));
		h1_Erec_1pi_e->SetBinContent(i,h1_Erec_1pi_e->GetBinContent(i)/h1_Erec_1pi_e->GetBinWidth(i));
		h1_E_rec->SetBinContent(i,h1_E_rec->GetBinContent(i)/h1_E_rec->GetBinWidth(i));


		h1_Erec_0pi_e->SetBinError(i,h1_Erec_0pi_e->GetBinError(i)/h1_Erec_0pi_e->GetBinWidth(i));
		h1_Erec_0pi_sub_e->SetBinError(i,h1_Erec_0pi_sub_e->GetBinError(i)/h1_Erec_0pi_sub_e->GetBinWidth(i));
		h1_Erec_1pi_e_weight->SetBinError(i,h1_Erec_1pi_e_weight->GetBinError(i)/h1_Erec_1pi_e_weight->GetBinWidth(i));
		h1_Erec_2pi_e_weight->SetBinError(i,h1_Erec_2pi_e_weight->GetBinError(i)/h1_Erec_2pi_e_weight->GetBinWidth(i));
		h1_Erec_3pi_e_weight->SetBinError(i,h1_Erec_3pi_e_weight->GetBinError(i)/h1_Erec_3pi_e_weight->GetBinWidth(i));
		h1_Erec_4pi_e_weight->SetBinError(i,h1_Erec_4pi_e_weight->GetBinError(i)/h1_Erec_4pi_e_weight->GetBinWidth(i));
		h1_Erec_1pi_e->SetBinError(i,h1_Erec_1pi_e->GetBinError(i)/h1_Erec_1pi_e->GetBinWidth(i));
		h1_E_rec->SetBinError(i,h1_E_rec->GetBinError(i)/h1_E_rec->GetBinWidth(i));

	}

	// ---------------------------------------------------------------------------------------------------- 

	TCanvas *c1=new TCanvas("c1","",1024,768);
	c1->SetTopMargin(0.05);
	c1->SetRightMargin(0.01);
	c1->SetLeftMargin(0.15);
	c1->SetBottomMargin(0.16);
	c1->cd();

	// Scale to get rid of the exponent

	double Scale = 1E6;

	for(int i=1;i<=N_E_bins;i++) {
	
		h1_E_rec->SetBinContent(i,h1_E_rec->GetBinContent(i)/Scale);
		h1_Erec_0pi_e->SetBinContent(i,h1_Erec_0pi_e->GetBinContent(i)/Scale);
		h1_Erec_0pi_sub_e->SetBinContent(i,h1_Erec_0pi_sub_e->GetBinContent(i)/Scale);

		h1_Erec_1pi_e_weight->SetBinContent(i,h1_Erec_1pi_e_weight->GetBinContent(i)/ Scale);
		h1_Erec_2pi_e_weight->SetBinContent(i,h1_Erec_2pi_e_weight->GetBinContent(i)/ Scale);
		h1_Erec_3pi_e_weight->SetBinContent(i,h1_Erec_3pi_e_weight->GetBinContent(i)/ Scale);
		h1_Erec_4pi_e_weight->SetBinContent(i,h1_Erec_4pi_e_weight->GetBinContent(i)/ Scale);
		h1_Erec_1pi_e->SetBinContent(i,h1_Erec_1pi_e->GetBinContent(i)/ Scale);

		h1_E_rec->SetBinError(i,h1_E_rec->GetBinError(i)/ Scale);
		h1_Erec_0pi_e->SetBinError(i,h1_Erec_0pi_e->GetBinError(i)/ Scale);
		h1_Erec_0pi_sub_e->SetBinError(i,h1_Erec_0pi_sub_e->GetBinError(i)/ Scale);
		
		h1_Erec_1pi_e_weight->SetBinError(i,h1_Erec_1pi_e_weight->GetBinError(i)/ Scale);
		h1_Erec_2pi_e_weight->SetBinError(i,h1_Erec_2pi_e_weight->GetBinError(i)/ Scale);
		h1_Erec_3pi_e_weight->SetBinError(i,h1_Erec_3pi_e_weight->GetBinError(i)/ Scale);
		h1_Erec_4pi_e_weight->SetBinError(i,h1_Erec_4pi_e_weight->GetBinError(i)/ Scale);
		h1_Erec_1pi_e->SetBinError(i,h1_Erec_1pi_e->GetBinError(i)/ Scale);

	}

	// ---------------------------------------------------------------------------------------------------- 

	gStyle->SetErrorX(kFALSE);

	h1_Erec_1pi_e_weight->SetAxisRange(xmin, xmax, "X"); 
	h1_Erec_1pi_e_weight->GetXaxis()->SetTitle("E_{QE}[GeV]");
	h1_Erec_1pi_e_weight->GetXaxis()->SetTitleOffset(1.);

	h1_Erec_1pi_e_weight->GetXaxis()->CenterTitle();
	h1_Erec_1pi_e_weight->GetXaxis()->SetTitleFont(FontStyle);
	h1_Erec_1pi_e_weight->GetXaxis()->SetLabelFont(FontStyle);
	h1_Erec_1pi_e_weight->GetXaxis()->SetTitleSize(TextSize);
	h1_Erec_1pi_e_weight->GetXaxis()->SetLabelSize(TextSize);
	h1_Erec_1pi_e_weight->GetXaxis()->SetNdivisions(Ndivisions);

	h1_Erec_1pi_e_weight->GetYaxis()->CenterTitle();
	h1_Erec_1pi_e_weight->GetYaxis()->SetTitleFont(FontStyle);
	h1_Erec_1pi_e_weight->GetYaxis()->SetLabelFont(FontStyle);
	h1_Erec_1pi_e_weight->GetYaxis()->SetTitleSize(TextSize);
	h1_Erec_1pi_e_weight->GetYaxis()->SetLabelSize(TextSize);
	h1_Erec_1pi_e_weight->GetYaxis()->SetNdivisions(Ndivisions);

	h1_Erec_1pi_e_weight->SetLineColor(46);

	h1_Erec_1pi_e_weight->GetYaxis()->SetTitle("Weighted Events / GeV");
	h1_Erec_1pi_e_weight->SetMarkerStyle(20);

	h1_Erec_1pi_e->SetMarkerStyle(20);
	h1_Erec_2pi_e_weight->SetMarkerStyle(20);
	h1_Erec_3pi_e_weight->SetMarkerStyle(20);
	h1_Erec_1pi_e_weight->SetMarkerSize(SizeMarker);
	h1_Erec_1pi_e->SetMarkerSize(SizeMarker);
	h1_Erec_2pi_e_weight->SetMarkerSize(SizeMarker);
	h1_Erec_3pi_e_weight->SetMarkerSize(SizeMarker);
	h1_Erec_1pi_e_weight->SetMarkerColor(46);
	h1_Erec_1pi_e->SetMarkerColor(38);
	h1_Erec_2pi_e_weight->SetMarkerColor(42);
	h1_Erec_3pi_e_weight->SetMarkerColor(28);

	h1_Erec_1pi_e_weight->GetYaxis()->SetRangeUser(0.,1.1*h1_Erec_1pi_e->GetMaximum());
	h1_Erec_1pi_e_weight->Draw("e");

	// -------------------------------------------------------------------------------------------

	// Nature data release

	int nbins = h1_Erec_1pi_e_weight->GetXaxis()->GetNbins();

	double Min = 1.8, Max = 5.1;
	myTxtFile << "Pad (a)" << endl;

	myTxtFile << "Undetected 1pi/gamma" << endl;

	int counter = 0;

	for (int i = 0; i < nbins; i++) {

		double BinContent = h1_Erec_1pi_e_weight->GetBinContent(i+1);
		double BinError = h1_Erec_1pi_e_weight->GetBinError(i+1);
		double BinLowEdge = h1_Erec_1pi_e_weight->GetBinLowEdge(i+1);			
		double BinWidth = h1_Erec_1pi_e_weight->GetBinWidth(i+1);
		double PreviousBinWidth = h1_Erec_1pi_e_weight->GetBinWidth(i);
		double BinHighEdge = BinLowEdge + BinWidth;

		if (BinLowEdge > Min-PreviousBinWidth && BinLowEdge < Max) {
			myTxtFile << counter+1 << ";" << BinLowEdge << ";" << BinHighEdge << ";" << BinContent << ";" << BinError << endl;
			counter++;
		}

	}

	myTxtFile << endl;

	// -------------------------------------------------------------------------------------------

	h1_Erec_1pi_e->Draw("e same");

	myTxtFile << "Detected 1pi/gamma" << endl;

	counter = 0;

	for (int i = 0; i < nbins; i++) {

		double BinContent = h1_Erec_1pi_e->GetBinContent(i+1);
		double BinError = h1_Erec_1pi_e->GetBinError(i+1);
		double BinLowEdge = h1_Erec_1pi_e->GetBinLowEdge(i+1);			
		double BinWidth = h1_Erec_1pi_e->GetBinWidth(i+1);
		double PreviousBinWidth = h1_Erec_1pi_e->GetBinWidth(i);
		double BinHighEdge = BinLowEdge + BinWidth;

		if (BinLowEdge > Min-PreviousBinWidth && BinLowEdge < Max) {
			myTxtFile << counter+1 << ";" << BinLowEdge << ";" << BinHighEdge << ";" << BinContent << ";" << BinError << endl;
			counter++;
		}

	}

	myTxtFile << endl;

	// -------------------------------------------------------------------------------------------

	h1_Erec_2pi_e_weight->SetLineColor(42);
	h1_Erec_2pi_e_weight->Draw("P e Same");

	myTxtFile << "Undetected 2pi/gamma" << endl;

	counter = 0;

	for (int i = 0; i < nbins; i++) {

		double BinContent = h1_Erec_2pi_e_weight->GetBinContent(i+1);
		double BinError = h1_Erec_2pi_e_weight->GetBinError(i+1);
		double BinLowEdge = h1_Erec_2pi_e_weight->GetBinLowEdge(i+1);			
		double BinWidth = h1_Erec_2pi_e_weight->GetBinWidth(i+1);
		double PreviousBinWidth = h1_Erec_2pi_e_weight->GetBinWidth(i);
		double BinHighEdge = BinLowEdge + BinWidth;

		if (BinLowEdge > Min-PreviousBinWidth && BinLowEdge < Max) {
			myTxtFile << counter+1 << ";" << BinLowEdge << ";" << BinHighEdge << ";" << BinContent << ";" << BinError << endl;
			counter++;
		}

	}

	myTxtFile << endl;

	// -------------------------------------------------------------------------------------------

	h1_Erec_3pi_e_weight->SetLineColor(28);
	//h1_Erec_3pi_e_weight->Draw("P e Same");

	//h1_Erec_4pi_e_weight->SetLineColor(6);

	lat1->SetTextSize(TextSize);
	lat1->SetTextFont(FontStyle);
	
	lat1->SetTextColor(38);
	lat1->DrawLatex(0.56, 0.85, "Detected 1 #pi^{#pm}/#gamma");
	lat1->SetTextSize(TextSize);
	lat1->SetTextColor(46);
	lat1->DrawLatex(0.56, 0.78, "Undetected 1#pi^{#pm}/#gamma (-)");
	lat1->SetTextColor(42);
	lat1->DrawLatex(0.56, 0.71, "Undetected 2#pi^{#pm}/#gamma (+)");
	lat1->SetTextColor(26);
	//lat1->DrawLatex(0.58, 0.64, "Undetected 3#pi^{#pm}/#gamma (+)");
	lat1->SetTextColor(6);

	lat1->SetTextColor(1);
	lat1->SetTextSize(TextSize);
	lat1->SetTextFont(FontStyle);
	lat1->SetTextSize(TextSize);
	lat1->DrawLatex(0.17, 0.85, "^{56}Fe");

	TLatex *EQEOne = new TLatex(); 
	EQEOne->SetTextFont(FontStyle); 
	EQEOne->SetTextColor(kBlack); 
	EQEOne->SetTextSize(TextSize);
	EQEOne->DrawLatexNDC(0.91,0.89,"(a)");

	c1->SaveAs("../../../myPlots/pdf/NoxBCut/"+version+target+"/FigExtData11_EQE_Inclusive_Subtractions_"+target+"_SuSav2.pdf");
	c1->SaveAs("../../../myPlots/pdf/NoxBCut/"+version+target+"/FigExtData11_EQE_Inclusive_Subtractions_"+target+"_SuSav2.eps");

	// -------------------------------------------------------------------------------------------

	TCanvas *c2 = new TCanvas("c2","",1024,768);
	c2->cd();

	h1_E_rec->SetAxisRange(xmin, xmax, "X"); 
	h1_E_rec->GetXaxis()->SetTitle("E_{QE} [GeV]");
	h1_E_rec->GetYaxis()->SetTitle("Weighted Events / GeV");
	h1_E_rec->GetXaxis()->SetTitleOffset(1.);
	h1_E_rec->SetMinimum(0);
	h1_E_rec->UseCurrentStyle();
	h1_E_rec->SetMarkerStyle(20);
	h1_E_rec->SetMarkerSize(SizeMarker);
	h1_E_rec->SetMarkerColor(38);
	//h1_E_rec->SetMinimum(3*h1_Erec_3pi_e_weight->GetMinimum());

	h1_E_rec->GetXaxis()->CenterTitle();
	h1_E_rec->GetXaxis()->SetTitleFont(FontStyle);
	h1_E_rec->GetXaxis()->SetLabelFont(FontStyle);
	h1_E_rec->GetXaxis()->SetTitleSize(TextSize);
	h1_E_rec->GetXaxis()->SetLabelSize(TextSize);

	h1_E_rec->GetYaxis()->CenterTitle();
	h1_E_rec->GetYaxis()->SetTitleFont(FontStyle);
	h1_E_rec->GetYaxis()->SetLabelFont(FontStyle);
	h1_E_rec->GetYaxis()->SetTitleSize(TextSize);
	h1_E_rec->GetYaxis()->SetLabelSize(TextSize);
	h1_E_rec->GetYaxis()->SetTitleOffset(1.);

//	h1_E_rec->SetAxisRange(3*h1_Erec_3pi_e_weight->GetMinimum(),1.05*h1_E_rec->GetMaximum(),"Y");
	h1_E_rec->SetAxisRange(0,1.17*h1_E_rec->GetMaximum(),"Y");
	h1_E_rec->Draw("e");

	// ---------------------------------------------------------------------

	myTxtFile << "Pad (b)" << endl;
	myTxtFile << "No cuts" << endl;

	counter = 0;

	for (int i = 0; i < nbins; i++) {

		double BinContent = h1_E_rec->GetBinContent(i+1);
		double BinError = h1_E_rec->GetBinError(i+1);
		double BinLowEdge = h1_E_rec->GetBinLowEdge(i+1);			
		double BinWidth = h1_E_rec->GetBinWidth(i+1);
		double PreviousBinWidth = h1_E_rec->GetBinWidth(i);
		double BinHighEdge = BinLowEdge + BinWidth;

		if (BinLowEdge > Min-PreviousBinWidth && BinLowEdge < Max) {
			myTxtFile << counter+1 << ";" << BinLowEdge << ";" << BinHighEdge << ";" << BinContent << ";" << BinError << endl;
			counter++;
		}

	}

	myTxtFile << endl;

	// ---------------------------------------------------------------------
	
	h1_Erec_0pi_e->SetAxisRange(xmin, xmax, "X"); 
	h1_Erec_0pi_e->SetLineColor(46);
	h1_Erec_0pi_e->GetXaxis()->SetTitleOffset(1.2);
	h1_Erec_0pi_e->SetMarkerStyle(20);
	h1_Erec_0pi_e->SetMarkerSize(SizeMarker);
	h1_Erec_0pi_e->SetMarkerColor(46);
	h1_Erec_0pi_e->Draw("e same");

	myTxtFile << "No detected pi/gamma" << endl;

	counter = 0;

	for (int i = 0; i < nbins; i++) {

		double BinContent = h1_Erec_0pi_e->GetBinContent(i+1);
		double BinError = h1_Erec_0pi_e->GetBinError(i+1);
		double BinLowEdge = h1_Erec_0pi_e->GetBinLowEdge(i+1);			
		double BinWidth = h1_Erec_0pi_e->GetBinWidth(i+1);
		double PreviousBinWidth = h1_Erec_0pi_e->GetBinWidth(i);
		double BinHighEdge = BinLowEdge + BinWidth;

		if (BinLowEdge > Min-PreviousBinWidth && BinLowEdge < Max) {
			myTxtFile << counter+1 << ";" << BinLowEdge << ";" << BinHighEdge << ";" << BinContent << ";" << BinError << endl;
			counter++;
		}

	}

	myTxtFile << endl;

	// ---------------------------------------------------------------------

	h1_Erec_0pi_sub_e->SetLineColor(42);
	h1_Erec_0pi_sub_e->SetMarkerStyle(20);
	h1_Erec_0pi_sub_e->SetMarkerSize(SizeMarker);
	h1_Erec_0pi_sub_e->SetMarkerColor(42);

	h1_Erec_0pi_sub_e->Draw("e Same");

	myTxtFile << "Subtract undetected pi/gamma" << endl;

	counter = 0;

	for (int i = 0; i < nbins; i++) {

		double BinContent = h1_Erec_0pi_sub_e->GetBinContent(i+1);
		double BinError = h1_Erec_0pi_sub_e->GetBinError(i+1);
		double BinLowEdge = h1_Erec_0pi_sub_e->GetBinLowEdge(i+1);			
		double BinWidth = h1_Erec_0pi_sub_e->GetBinWidth(i+1);
		double PreviousBinWidth = h1_Erec_0pi_sub_e->GetBinWidth(i);
		double BinHighEdge = BinLowEdge + BinWidth;

		if (BinLowEdge > Min-PreviousBinWidth && BinLowEdge < Max) {
			myTxtFile << counter+1 << ";" << BinLowEdge << ";" << BinHighEdge << ";" << BinContent << ";" << BinError << endl;
			counter++;
		}

	}

	myTxtFile << endl;

	// ---------------------------------------------------------------------

	lat1->SetTextColor(38);
	lat1->DrawLatex(0.51, 0.89, "No cuts");

	lat1->SetTextColor(46);
	lat1->DrawLatex(0.51, 0.82, "No detected #pi^{#pm}/#gamma");
	 
	lat1->SetTextColor(42);
	lat1->DrawLatex(0.51, 0.75, "Subtract undetected #pi^{#pm}/#gamma");
		
	lat1->SetTextColor(1);
	lat1->SetTextSize(0.09);
	lat1->SetTextSize(TextSize);
	lat1->DrawLatex(0.17, 0.85, "^{56}Fe");

	TLatex *EQETwo = new TLatex(); 
	EQETwo->SetTextFont(FontStyle); 
	EQETwo->SetTextColor(kBlack); 
	EQETwo->SetTextSize(TextSize);
	EQETwo->DrawLatexNDC(0.91,0.89,"(b)");

	c2->SaveAs("../../../myPlots/pdf/NoxBCut/"+version+target+"/FigExtData11_EQE_Exclusive_Subtractions_"+target+"_SuSav2.pdf");
	c2->SaveAs("../../../myPlots/pdf/NoxBCut/"+version+target+"/FigExtData11_EQE_Exclusive_Subtractions_"+target+"_SuSav2.eps");

	// ---------------------------------------------------------------------

	TCanvas *ctot = new TCanvas("ctot","",2048,768);
	ctot->cd();

	TPad* left = new TPad("left","left",0.,0.,0.5,1.,21);
	left->SetFillColor(kWhite); 
	left->Draw();
	left->cd();

	c1->DrawClonePad();

	ctot->cd();
	TPad* right = new TPad("right","right",0.5,0.,1.,1.,21);
	right->SetFillColor(kWhite); 
	right->Draw();
	right->cd();

	c2->DrawClonePad();

	ctot->SaveAs("../../../myPlots/pdf/NoxBCut/"+version+target+"/Merged_FigExtData7_EQE_Exclusive_Subtractions_"+target+"_SuSav2.pdf");
	ctot->SaveAs("../../../myPlots/pdf/NoxBCut/"+version+target+"/Merged_FigExtData7_EQE_Exclusive_Subtractions_"+target+"_SuSav2.eps");

	// --------------------------------------------------------------------- 

}

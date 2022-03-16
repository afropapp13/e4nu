#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <TF1.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

	//----------------------------------------//

	// All units are in GeV

	const double m_pimi = 0.139570, m_pipl = 0.139570, m_pion = 0.139570;
	const double m_prot = 0.9382720813, m_neut = 0.939565;
	const double H3_bind_en = 0.008481, He4_bind_en = 0.0283, C12_bind_en = 0.09215, B_bind_en = 0.0762;
	const double He3_bind_en = 0.0077, D2_bind_en = 0.00222, Fe_bind_en = 0.49226, Mn_bind_en = 0.4820764;
	const double e_mass = 0.000510998;

	const double fine_struc_const = 0.007297;
	const double ns_to_s = 1.0E-9;
	const Double_t c = 2.99792E+10;

	static const Float_t par_EcUVW[6][3] = {{60, 360, 400}, {55, 360, 400}, {50, 363, 400}, {52, 365, 396}, {60, 360, 398}, {50, 362, 398}};

	const double PhiOpeningAngle = 6;

	const double MinThetaProton = 12.;
	const double MinThetaPiPlus = 12.;
	const double MinThetaPiMinus = 0.;
	const double MinThetaGamma = 8.;

	const double CenterFirstSector = 30;
	const double CenterSecondSector = 90;
	const double CenterThirdSector = 150;
	const double CenterFourthSector = 210;
	const double CenterFifthSector = 270;
	const double CenterSixthSector = 330;

	const int RotCounterLimit = 100;
	
	//----------------------------------------//
	
	// Interaction labels
	
	const std::vector<int> InteBreakColors{kBlack,kBlue-5,kYellow+1,kOrange+7,kRed+1,kBlue};	
	const std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};
	const int NInte = InteractionLabels.size();	
	
	//----------------------------------------//
	
	const int NInt = 6; // All Interactions = 0, QE = 1, MEC = 2, RES = 3, DIS = 4, COH = 5
	
	const std::vector<double> TwoDArrayNBinsDeltaPT{0.0,0.2,0.4,10.0};	
	const std::vector<double> TwoDArrayNBinsDeltaAlphaT{0.0,45.0,90.0,135.0,180.0};
	const std::vector<double> TwoDArrayNBinsDeltaPtx{-5.,-0.15,0.15,5.};
	const std::vector<double> TwoDArrayNBinsDeltaPty{-5.,-0.15,0.15,5.};
	
	//----------------------------------------//
	
	const std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices{ 
												// DeltaPT < 0.2 GeV/c
												{0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.},
												// 0.2 < DeltaPT < 0.4 GeV/c
												{0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.},
												// DeltaPT > 0.4 GeV/c
												{0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.}		
											};

	//----------------------------------------//
	
	const std::vector< std::vector<double> > TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices{ 
												// DeltaAlphaT < 45 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.8,1.},
												// 45 < DeltaAlphaT < 90 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.8,1.},
												// 90 < DeltaAlphaT < 135 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.43,0.52,0.6,0.65,0.7,0.75,0.8,0.84,0.88,0.92,0.97,1.12,1.2},
												// DeltaAlphaT > 135 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.63,0.71,0.79,0.87,0.95,1.1}		
											};	
											
	//----------------------------------------//
	
	const std::vector< std::vector<double> > TwoDArrayNBinsDeltaPtyInDeltaPtxSlices{ 
											// DeltaPtx < -0.15 GeV/c
											{-0.8,-0.72,-0.64,-0.58,-0.5,-0.42,-0.36,-0.28,-0.2,-0.15,-0.1,-0.05,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5},
											// -0.15 < DeltaPtx < 0.15 GeV/c
											{-0.8,-0.74,-0.67,-0.6,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5},
											// DeltaPtx > 0.15 GeV/c
											{-0.8,-0.71,-0.62,-0.55,-0.48,-0.41,-0.34,-0.27,-0.2,-0.15,-0.1,-0.05,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5}
										};

	//----------------------------------------//
	
	const std::vector< std::vector<double> > TwoDArrayNBinsDeltaPtxInDeltaPtySlices{ 
											// DeltaPty < -0.15 GeV/c
											{-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75},
											// -0.15 < DeltaPty < 0.15 GeV/c
											{-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75},
											// DeltaPty > 0.15 GeV/c
											{-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75}
										};																																		
	
	//----------------------------------------//		

#endif

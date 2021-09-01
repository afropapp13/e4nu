
//Class created by Afroditi Papadopoulou (apapadop@mit.edu)

// _________________________________________________________________________________________________________________________________________________

#ifndef STV_TOOLS_CXX
#define STV_TOOLS_CXX

#include <TMath.h>

#include "STV_Tools.h"

// __________________________________________________________________________________________________________________________________________________

STV_Tools::STV_Tools(TVector3 MuonVector,TVector3 ProtonVector, double MuonEnergy, double ProtonEnergy) {

	double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938272, NeutronMass_GeV = 0.939565; // GeV
	double DeltaM2 = TMath::Power(NeutronMass_GeV,2.) - TMath::Power(ProtonMass_GeV,2.);	
	double BE = 0.04; // GeV	
			
	// STV Calculation		
			
	TVector3 MuonVectorTrans;
	MuonVectorTrans.SetXYZ(MuonVector.X(),MuonVector.Y(),0.);
	double MuonVectorTransMag = MuonVectorTrans.Mag();

	TVector3 MuonVectorLong;
	MuonVectorLong.SetXYZ(0.,0.,MuonVector.Z());
	double MuonVectorLongMag = MuonVectorLong.Mag();
	
	TLorentzVector MuonLorentzVector(MuonVector,MuonEnergy);
//	double MuonKE = MuonEnergy - MuonMass_GeV;	
			
	TVector3 ProtonVectorTrans;
	ProtonVectorTrans.SetXYZ(ProtonVector.X(),ProtonVector.Y(),0.);
	double ProtonVectorTransMag = ProtonVectorTrans.Mag();

	TVector3 ProtonVectorLong;
	ProtonVectorLong.SetXYZ(0.,0.,ProtonVector.Z());
	double ProtonVectorLongMag = ProtonVectorLong.Mag();
	
	TLorentzVector ProtonLorentzVector(ProtonVector,ProtonEnergy);
	double ProtonKE = ProtonEnergy - ProtonMass_GeV;		

	TVector3 PtVector = MuonVectorTrans + ProtonVectorTrans;

	fPt = PtVector.Mag();

	fDeltaAlphaT = TMath::ACos( (- MuonVectorTrans * PtVector) / ( MuonVectorTransMag * fPt ) ) * 180./TMath::Pi();
	if (fDeltaAlphaT > 180.) { fDeltaAlphaT -= 180.; }
	if (fDeltaAlphaT < 0.) { fDeltaAlphaT += 180.; }

	fDeltaPhiT = TMath::ACos( (- MuonVectorTrans * ProtonVectorTrans) / ( MuonVectorTransMag * ProtonVectorTransMag ) ) * 180./TMath::Pi();
	if (fDeltaPhiT > 180.) { fDeltaPhiT -= 180.; }
	if (fDeltaPhiT < 0.) { fDeltaPhiT += 180.; }

	// -------------------------------------------------------------------------------------------------------------------------

	// Calorimetric Energy Reconstruction

	fECal = MuonEnergy + ProtonKE + BE; // GeV

	// QE Energy Reconstruction

	double EQENum = 2 * (NeutronMass_GeV - BE) * MuonEnergy - (BE*BE - 2 * NeutronMass_GeV *BE + MuonMass_GeV * MuonMass_GeV + DeltaM2);
	double EQEDen = 2 * ( NeutronMass_GeV - BE - MuonEnergy + MuonVector.Mag() * MuonVector.CosTheta() );
	fEQE = EQENum / EQEDen;

	// Reconstructed Q2

	TLorentzVector nuLorentzVector(0.,0.,fECal,fECal);
	TLorentzVector qLorentzVector = nuLorentzVector - MuonLorentzVector;
	fQ2 = - qLorentzVector.Mag2();
	
	// -------------------------------------------------------------------------------------------------------------------------
	
	// Light Cone Variables
	
	TLorentzVector MissLorentzVector = MuonLorentzVector + ProtonLorentzVector - nuLorentzVector;
	
	fEMiss = TMath::Abs(MissLorentzVector.E());
	fPMiss = (MissLorentzVector.Vect()).Mag();
	
//	fPMissMinus = fEMiss - MissLorentzVector.Z();

	// Suggestion from Jackson to avoid Ecal assumption
	fPMissMinus = (MuonEnergy - MuonVector.Z()) + (ProtonEnergy - ProtonVector.Z());
	
	double fkMissNum = ( TMath::Power(fPt,2.) + TMath::Power(ProtonMass_GeV,2.) );
	double fkMissDen = ( fPMissMinus * (2*ProtonMass_GeV - fPMissMinus) );
	
	double fkMiss2 = TMath::Power(ProtonMass_GeV,2.) * fkMissNum / fkMissDen - TMath::Power(ProtonMass_GeV,2.); // Jackson's GlueX note				

	fkMiss = sqrt(fkMiss2);

	fA = fPMissMinus / ProtonMass_GeV;

	// -------------------------------------------------------------------------------------------------------------------------

	// Minerva longitudinal & total variables

	// For the calcukation of the masses
	//https://journals.aps.org/prc/pdf/10.1103/PhysRevC.95.065501

	double MA = 22 * NeutronMass_GeV + 18 * ProtonMass_GeV - 0.34381; // GeV
	double MAPrime = MA - NeutronMass_GeV + 0.00556; // GeV, constant obtained from table I 

	// For the calculation of p_n, back to the Minerva PRL
	// https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.022504

	double R = MA + MuonVectorLongMag + ProtonVectorLongMag - MuonEnergy - ProtonEnergy; // Equation 8

	// Equation 7

	fPL = 0.5 * R - (MAPrime * MAPrime + fPt * fPt) / (2 * R);

	fPn = TMath::Sqrt( fPt * fPt + fPL * fPL );

	// https://journals.aps.org/prd/pdf/10.1103/PhysRevD.101.092001

//	fPtx = fPt * TMath::Sin(fDeltaAlphaT * TMath::Pi() / 180.);
//	fPty = fPt * TMath::Cos(fDeltaAlphaT * TMath::Pi() / 180.);

	TVector3 UnitZ(0,0,1);
	fPtx = ( UnitZ.Cross(MuonVectorTrans) ).Dot(PtVector) / MuonVectorTransMag;
	fPty = - (MuonVectorTrans).Dot(PtVector) / MuonVectorTransMag;

	// -------------------------------------------------------------------------------------------------------------------------

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnkMiss() {

	return fkMiss;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnEMiss() {

	return fEMiss;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPMissMinus() {

	return fPMissMinus;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPMiss() {

	return fPMiss;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPt() {

	return fPt;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPtx() {

	return fPtx;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPty() {

	return fPty;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPL() {

	return fPL;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnPn() {

	return fPn;
}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnDeltaAlphaT() {

	return fDeltaAlphaT;

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnDeltaPhiT() {

	return fDeltaPhiT;

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnECal() {

	return fECal;

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnEQE() {

	return fEQE;

}

// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnQ2() {

	return fQ2;

}
// __________________________________________________________________________________________________________________________________________________

double STV_Tools::ReturnA() {

	return fA;

}

// __________________________________________________________________________________________________________________________________________________

#endif

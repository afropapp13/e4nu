#!/bin/bash

xB=NoxBCut
#xB=xBCut

LocalPath=/home/afroditi/Dropbox/PhD/myCode/43rd_LightFrontNeutrinos/myFiles

GENIEAccount=apapadop@geniegpvm02.fnal.gov
GENIEPath=/genie/app/users/apapadop/e4v_LightFrontNeutrinos

## SuSAv2, w/ rotations, w/ radiation, w/ fiducials, w/ smearing
#GENIEOnlineExtension=SuSav2_Rad
#GENIELocalExtension=SuSav2_RadCorr_LFGM_UpdatedSchwinger

## SuSAv2, w/ rotations, w/o radiation, w/ fiducials, w/ smearing
#GENIEOnlineExtension=SuSav2
#GENIELocalExtension=SuSav2_NoRadCorr_LFGM

## G2018, w/ rotations, w/o radiation, w/ fiducials, w/ smearing
#GENIEOnlineExtension=G18_10a_02_11a
#GENIELocalExtension=hA2018_Final_NoRadCorr_LFGM

## SuSAv2, 1p0pi, w/ radiation, w/ fiducials, w/ smearing
#GENIEOnlineExtension=SuSav2_Rad
#GENIELocalExtension=SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger

## SuSAv2, 1p0pi, w/o radiation, w/ fiducials, w/ smearing
#GENIEOnlineExtension=SuSav2
#GENIELocalExtension=SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc

## G2018, 1p0pi, w/o radiation, w/ fiducials, w/ smearing
#GENIEOnlineExtension=G18_10a_02_11a
#GENIELocalExtension=hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc

## SuSAv2, 1p0pi, w/o radiation, w/o fiducials, w/o smearing
#GENIEOnlineExtension=SuSav2
#GENIELocalExtension=SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc

## G2018, 1p0pi, w/o radiation, w/o fiducials, w/o smearing
#GENIEOnlineExtension=G18_10a_02_11a
#GENIELocalExtension=hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc

## SuSAv2, 1p0pi, w/o radiation, w/o fiducials, w/ smearing
GENIEOnlineExtension=SuSav2
GENIELocalExtension=SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing

## G2018, 1p0pi, w/o radiation, w/o fiducials, w/ smearing
#GENIEOnlineExtension=G18_10a_02_11a
#GENIELocalExtension=hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing

#-----------------------------------------------

##### Genie Samples

echo "$LocalPath/1_161/${GENIELocalExtension}/$xB/12C_1_161_${GENIELocalExtension}_Plots_FSI_em.root"

scp $GENIEAccount:$GENIEPath/LFNeutrinos_genie_e2a_ep_C12_1161_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/1_161/${GENIELocalExtension}/$xB/12C_1_161_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/LFNeutrinos_genie_e2a_ep_4He_2261_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/2_261/${GENIELocalExtension}/$xB/4He_2_261_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/LFNeutrinos_genie_e2a_ep_C12_2261_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/2_261/${GENIELocalExtension}/$xB/12C_2_261_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/LFNeutrinos_genie_e2a_ep_56Fe_2261_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/2_261/${GENIELocalExtension}/$xB/56Fe_2_261_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/LFNeutrinos_genie_e2a_ep_4He_4461_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/4_461/${GENIELocalExtension}/$xB/4He_4_461_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/LFNeutrinos_genie_e2a_ep_C12_4461_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/4_461/${GENIELocalExtension}/$xB/12C_4_461_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/LFNeutrinos_genie_e2a_ep_56Fe_4461_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/4_461/${GENIELocalExtension}/$xB/56Fe_4_461_${GENIELocalExtension}_Plots_FSI_em.root

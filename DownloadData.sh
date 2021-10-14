################################################

xB=NoxBCut
#xB=xBCut

JLabAccount=apapadop@ftp.jlab.org
JLabPath=/u/home/apapadop/e4v_Inclusive
DataExtension=Pinned_Data_Final
LocalPath=/home/afroditi/Dropbox/PhD/myCode/44th_e4v_Inclusive/myFiles

##### Data Samples

scp $JLabAccount:$JLabPath/LFNeutrinos_data_e2a_ep_C12_1161_neutrino6_united4_radphot_test.root $LocalPath/1_161/$DataExtension/$xB/12C_1_161_${DataExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/LFNeutrinos_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/4He_2_261_${DataExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/LFNeutrinos_data_e2a_ep_C12_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/12C_2_261_${DataExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/LFNeutrinos_data_e2a_ep_56Fe_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/56Fe_2_261_${DataExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/LFNeutrinos_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/4He_4_461_${DataExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/LFNeutrinos_data_e2a_ep_C12_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/12C_4_461_${DataExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/LFNeutrinos_data_e2a_ep_56Fe_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/56Fe_4_461_${DataExtension}_Plots_FSI_em.root

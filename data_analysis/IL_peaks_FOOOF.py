#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 10:31:41 2022

@author: sakkol
"""
## Multi FOOOFing
import numpy as np
from scipy.io import loadmat # savemat
from fooof import FOOOF
#from fooof import FOOOFGroup, FOOOF
import mat73
import os

############################### DEFINE NAMES ##################################
sbj_ID = 'NS164'
blockname = 'IsochronousListening'
############################### DEFINE NAMES ##################################

## Arrange names
proj_folder = '/media/sakkol/HDD1/HBML/PROJECTS_DATA/IsochronousListening/'
ecog_avg_name = '%s%s/iEEG_data/%s/%s_ecog_avg.mat' %(proj_folder,sbj_ID,blockname,blockname)
info_name = '%s%s/iEEG_data/%s/%s_info.mat' %(proj_folder,sbj_ID,blockname,blockname)
trial_nos_name = '%s%s/iEEG_data/%s/%s_trial_nos.mat' %(proj_folder,sbj_ID,blockname,blockname)
elec_TF_folder = '%s%s/iEEG_data/%s/elec_TF/' %(proj_folder,sbj_ID,blockname)
picfolder = '%s%s/results/FOOOF/' %(proj_folder,sbj_ID)
blockfolder = '%s%s/results/FOOOF/' %(proj_folder,sbj_ID)


##Unpack data from dictionary, and squeeze into numpy arrays
data_dict = mat73.loadmat(ecog_avg_name)
info_dict = mat73.loadmat(info_name)
trialnos_dict = loadmat(trial_nos_name)
Label=data_dict['ecog_avg']['ftrip']['label']

loop1={'iso','a'}
loop2={4,3,5}
loop3={'sentence','scrambled'}

for el in range(0, len(Label)-1):
    power_spectra_name = '%s%s/iEEG_data/%s/elec_TF/%s_%s_mtmconvol_fourier.mat' %(proj_folder,sbj_ID,blockname,Label[el][0],blockname)
    data_dict = mat73.loadmat(power_spectra_name)
    
    freqs = np.squeeze(data_dict['epoched_wlt']['freq']).astype('float')
    
    for l1 in loop1:
        
        for l2 in loop2:
            
            for l3 in loop3:
                
                print([l1+'_'+str(l2)+'_'+l3])
                curr_tr=trialnos_dict['trial_nos'][l1+'_'+str(l2)+'_'+l3][0][0][0] - 1
                psds = np.squeeze(np.nanmean(np.absolute(data_dict['epoched_wlt']['fourierspctrm'][curr_tr,:,:]),axis=(0,2))).astype('float')
    
                # Initialize FOOOF object
                fm = FOOOF()
                
                # Fit the FOOOF model on all PSDs, and report
                fm.report(freqs, psds, [0, 5])
                
                # Save out a specific FOOOF measure of interest - for example, slopes
                # These are not very useful (SA)
                #FOOOF_results = fm.get_results()
                #savemat(fmsave+'.mat', {'FOOOF_results' : FOOOF_results})
                
                # Save out fooof results to json file
                #  There is a utility file to load this json file directly into Matlab
                savefolder='%s%s/results/FOOOF/%s' %(proj_folder,sbj_ID,l1+'_'+str(l2)+'_'+l3)
                if not os.path.exists(savefolder):
                    os.mkdir(savefolder)
                savename='%s/%s_%s_slopes' %(savefolder,Label[el][0],blockname)
                fm.save(savename+'.json', save_results=True)
                fm.save_report(savename+'.png')
    
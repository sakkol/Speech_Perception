#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 15:35:31 2022

@author: sakkol
"""
import gensim
from scipy.io import savemat, loadmat


model = gensim.models.keyedvectors.load_word2vec_format('/media/sakkol/HDD1/HBML/PROJECTS_DATA/word2vec/GoogleNews-vectors-negative300.bin', binary=True) 

tmp = loadmat('/media/sakkol/HDD1/HBML/PROJECTS_DATA/word2vec/WordList.mat')
WordList = tmp['WordList']
WordList = list(WordList[w][0][0] for w in range(len(WordList)))

word2vec_out=list()
for w in WordList:
    word2vec_out.append(model[w])

savemat('/media/sakkol/HDD1/HBML/PROJECTS_DATA/word2vec/word2vec_out.mat',{'word2vec_out':word2vec_out,'WordList':WordList})

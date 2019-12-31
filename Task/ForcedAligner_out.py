#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:55:20 2019

@author: sakkol
"""
import os
from pathlib import Path
from p2fa import align

txtdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-English/txt'
wavdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-English/wav'
TGdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-English/TextGrid'
matdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-English/mat'

txtdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-English/txt'
wavdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-English/wav'
TGdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-English/TextGrid'
matdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-English/mat'

# Define language here: Spanish or English
lang = 'Spanish'
lang = 'English'

# Get items in directory with text files here
items = os.listdir(txtdirname)

# Gather text file names from directory
txtfilelist = []
for names in items:
    if names.endswith(".txt"):
        txtfilelist.append(names)


# Iterate through text files, performing aligner and saving Praat as we go
# for file in txtfilelist:
for ff in range(0,len(txtfilelist)):
    
    curr_text = txtfilelist[ff]
    curr_wav = curr_text.replace('txt', "wav")
    curr_TextGird = curr_text.replace('txt', "TextGrid")

    print([str(ff)+': '+curr_text])
    
    phoneme_alignments, word_alignments = align.align(wavdirname+'/'+curr_wav, txtdirname+'/'+curr_text,TGdirname+'/'+curr_TextGird)
            
    with open(Path(matdirname,curr_text), 'w') as fp:
        fp.writelines('Phonemes\n')
        for pp in range(1,len(phoneme_alignments)):
            fp.write(' '.join('%s' % x for x in phoneme_alignments[pp]))
            fp.writelines('\n')
        fp.writelines('Words\n')
        for pp in range(1,len(word_alignments)):
            fp.write(' '.join('%s' % x for x in word_alignments[pp]))
            fp.writelines('\n')
            


## For Spanish (unfortunately P2FA doesn't have Spanish support)
txtdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-Spanish/txt'
wavdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-Spanish/wav'
TGdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-Spanish/TextGrid'
matdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-Spanish/mat'

txtdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-Spanish/txt'
wavdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-Spanish/wav'
TGdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-Spanish/TextGrid'
matdirname = r'/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-Spanish/mat'

# Define language here: Spanish or English
lang = 'Spanish'
lang = 'English'

# Get items in directory with text files here
items = os.listdir(txtdirname)

# Gather text file names from directory
txtfilelist = []
for names in items:
    if names.endswith(".txt"):
        txtfilelist.append(names)


# Iterate through text files, performing aligner and saving Praat as we go
# for file in txtfilelist:
for ff in range(0,len(txtfilelist)):
    
    curr_text = txtfilelist[ff]
    curr_wav = curr_text.replace('txt', "wav")
    curr_TextGird = curr_text.replace('txt', "TextGrid")

    print([str(ff)+': '+curr_text])
    
    phoneme_alignments, word_alignments = align.align(wavdirname+'/'+curr_wav, txtdirname+'/'+curr_text,TGdirname+'/'+curr_TextGird)
            
    with open(Path(matdirname,curr_text), 'w') as fp:
        fp.writelines('Phonemes\n')
        for pp in range(1,len(phoneme_alignments)):
            fp.write(' '.join('%s' % x for x in phoneme_alignments[pp]))
            fp.writelines('\n')
        fp.writelines('Words\n')
        for pp in range(1,len(word_alignments)):
            fp.write(' '.join('%s' % x for x in word_alignments[pp]))
            fp.writelines('\n')
            
            
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:37:15 2020

@author: sakkol
"""

"""Synthesizes speech from the input string of text or ssml.

Note: ssml must be well-formed according to:
    https://www.w3.org/TR/speech-synthesis/
"""
#from google.cloud import storage
import os
from pathlib import Path
from google.cloud import texttospeech

## Instantiates a client
client = texttospeech.TextToSpeechClient.from_service_account_json(Path(r"/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/My First Project-4597d3a51499.json"))


## DEFINE ALL PARAMETERS HERE
# Either folder name with text files
CORPUS = 'CST-Repeated'
dirname = r'/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/word-texts/English'
dirname = Path(dirname)
outdirname = r'/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/word_sounds/English'
outdirname = Path(outdirname)

# Define output voice gender
# Either 'M' or 'F' ; Male preferred
GENDER = 'F'
GENDER = 'M'

# Define language here: Spanish or English
lang = 'Spanish'
lang = 'English'

# Speed up by %130
# rate = 0.9
rate = 1.6





# Alter parameter dependent variables
if lang == 'Spanish':
    voice_name = 'es-ES-Standard-A'
    lang_code = 'es-ES'
    # Pitch up female by 3 semitones
    pitch = 3.0
elif lang == 'English':
    lang_code = 'en-US'
    if GENDER == 'F':
        voice_name = 'en-US-Wavenet-E' # D for male, E for female (Wavenet voices)
        # Pitch up female by 3 semitones
        pitch = 3.0
    elif GENDER == 'M':
        voice_name = 'en-US-Wavenet-D'
        # Pitch down male by 3 semitones
        pitch = -3.0
else:
    raise ValueError('Invalid gender entered')


# Build the voice request, select the language code ("en-US") and the ssml
# voice gender ("neutral")
voice = texttospeech.types.VoiceSelectionParams(
    language_code=lang_code,
    name=voice_name)

# Select the type of audio file you want returned
audio_config = texttospeech.types.AudioConfig(
    audio_encoding = texttospeech.enums.AudioEncoding.LINEAR16, #LPCM WAV
    speaking_rate = rate,
    pitch = pitch,
    sample_rate_hertz = 44100)

# Get items in directory with text files here
items = os.listdir(dirname)


# Gather text file names from directory
filelist = []
for names in items:
    if names.endswith(".txt"):
        filelist.append(names)

# Iterate through text files, performing TTS and saving audio as we go
for file in filelist:
    # Read in text
    if lang_code == 'es-ES':
        with open(Path(dirname,file), 'r', encoding='UTF-8') as text_file:
            text = text_file.read()
    else:
        with open(Path(dirname,file), 'r') as text_file:
            text = text_file.read()

    print(text)

    # Set the text input to be synthesized
    synthesis_input = texttospeech.types.SynthesisInput(text=text)

    # Perform the text-to-speech request on the text input with the selected
    # voice parameters and audio file type
    response = client.synthesize_speech(synthesis_input, voice, audio_config)

    name = file[0:-4]

    # The response's audio_content is binary.
    with open(Path(outdirname,name+'-'+GENDER+'.wav'), 'wb') as out:
        # Write the response to the output file.
        out.write(response.audio_content)
        print('Audio content written to file "'+name+'-'+GENDER+'.wav"')
        print('~~~')
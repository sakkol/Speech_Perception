"""Synthesizes speech from the input string of text or ssml.

Note: ssml must be well-formed according to:
    https://www.w3.org/TR/speech-synthesis/
"""
from google.cloud import texttospeech
import os
from pathlib import Path
from google.cloud import storage

# Explicitly use service account credentials by specifying the private key
# file.
storage_client = storage.Client.from_service_account_json('My First Project-dc1f01a6b01a.json')

# Make an authenticated API request
buckets = list(storage_client.list_buckets())
print(buckets)

# DEFINE ALL PARAMETERS HERE
# Define output voice gender
# Either 'M' or 'F'
GENDER = 'M'
# Either folder name with text files
CORPUS = 'CST-Repeated'
dirname = r'C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentences\Selected_Texts'
dirname = Path(dirname)
outdirname = r'C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentences\Speech'
outdirname = Path(outdirname)

# Alter parameter dependent variables
if GENDER == 'M':
    voice_name = 'en-US-Wavenet-D'
    # Pitch down male by 3 semitones
    pitch = -3.0
elif GENDER == 'F':
    voice_name = 'en-US-Wavenet-E'
    # Pitch up female by 3 semitones
    pitch = 3.0
else:
    raise ValueError('Invalid gender entered')
# Slow down speech to 90%
rate = 0.9 

## Instantiates a client
client = texttospeech.TextToSpeechClient()

# Build the voice request, select the language code ("en-US") and the ssml
# voice gender ("neutral")
voice = texttospeech.types.VoiceSelectionParams(
    language_code='en-US',
    name=voice_name) # D for male, E for female (Wavenet voices)

# Select the type of audio file you want returned
audio_config = texttospeech.types.AudioConfig(
    audio_encoding = texttospeech.enums.AudioEncoding.LINEAR16, #LPCM WAV
    speaking_rate = rate,
    pitch = pitch)

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
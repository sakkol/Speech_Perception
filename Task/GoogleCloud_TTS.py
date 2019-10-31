"""Synthesizes speech from the input string of text or ssml.

Note: ssml must be well-formed according to:
    https://www.w3.org/TR/speech-synthesis/
"""
#from google.cloud import storage
import os
from pathlib import Path

# set GOOGLE_APPLICATION_CREDENTIALS="C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\My First Project-dc1f01a6b01a.json"
# os.environ["GOOGLE_APPLICATION_CREDENTIALS"]=r"C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\My First Project-dc1f01a6b01a.json"
# Explicitly use service account credentials by specifying the private key
# file.
#storage_client = storage.Client.from_service_account_json('My First Project-dc1f01a6b01a.json')
#buckets = list(storage_client.list_buckets())
#print(buckets)

from google.cloud import texttospeech

## Instantiates a client
client = texttospeech.TextToSpeechClient.from_service_account_json(Path(r"/home/sakkol/Documents/Spanish_Matrix_Sentence/Version_1/My First Project-4597d3a51499.json"))


## DEFINE ALL PARAMETERS HERE
# Either folder name with text files
CORPUS = 'CST-Repeated'
dirname = r'/home/sakkol/Documents/Spanish_Matrix_Sentence/Version_2/Selected_20000'
dirname = Path(dirname)
outdirname = r'/home/sakkol/Documents/Spanish_Matrix_Sentence/Version_2/Google_TTS_Rate_20000_1.3'
outdirname = Path(outdirname)

# Define output voice gender
# Either 'M' or 'F'
GENDER = 'F'
GENDER = 'M'

# Define language here: Spanish or English
lang = 'Spanish'
lang = 'English'

# Slow down speech to 90% or Speed up by %130
rate = 0.9
rate = 1.3





# Alter parameter dependent variables
if lang == 'Spanish':
    voice_name = 'es-ES-Standard-A'
    lang_code = 'es-ES'
    # Pitch up female by 3 semitones
    pitch = 3.0
elif lang == 'English':
    lang_code = 'en-US'
    if GENDER == 'F':
        voice_name = 'en-US-Wavenet-E'
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
    if lang_code == 'es-ES':
        with open(Path(dirname,file), 'r', encoding='latin-1') as text_file:
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
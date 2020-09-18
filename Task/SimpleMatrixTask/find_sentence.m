function filename = find_sentence(sentence,main_stim_loc,speech_rate)
% returns exact location of speech sentence

% prepare sentence
repl_sentence = replace(sentence,' ','_');
repl_sentence = erase(repl_sentence,'.');

if contains(main_stim_loc,'Spanish')
    sentence_file = [repl_sentence '-F.wav'];
else
    sentence_file = [repl_sentence '-M.wav'];
end

filename = fullfile(main_stim_loc,['Sentences_Rate' num2str(speech_rate)],sentence_file);

if ~isfile(filename)
    error_msg = sprintf('The file: "%s" could not be found!',filename);
    errordlg(error_msg, 'File not found!');
end

end


function filename = find_word(word,main_stim_loc,language)
% returns exact location of speech word

% prepare word
repl_word = replace(word,' ','_');
repl_word = erase(repl_word,'.');

if strcmpi(language,'Spanish')
    word_file = [repl_word '-F.wav'];
else
    word_file = [repl_word '-M.wav'];
end

filename = fullfile(main_stim_loc,language,word_file);

if ~isfile(filename)
    error_msg = sprintf('The file: "%s" could not be found!',filename);
    errordlg(error_msg, 'File not found!');
end

end


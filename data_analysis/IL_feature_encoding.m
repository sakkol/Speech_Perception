%% Load the mTRF features

% load lexico-semantic features: word2vec
% load trial info and loop trials
% prepare word onset
% preparesentence onset
% 
Y = tsne(word2vec_out);
gscatter(Y(:,1),Y(:,2),zz)
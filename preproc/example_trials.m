load('TDT_data.mat');

% envelope of speech channel (Wav1)
% 2 examples from block B52
Wav1_beg_idx = 145700;
Wav1_end_idx = 148000;
Wav1_beg_idx = 368400;
Wav1_end_idx = 370500;

% 2 examples from block B65
Wav1_beg_idx = 61750;
Wav1_end_idx = 64500;
Wav1_beg_idx = 439750;
Wav1_end_idx = 442500;

% chop only that part
Wav1_part = double(data.streams.Wav1.data(1,Wav1_beg_idx:Wav1_end_idx));
Wav1_to_plot = resample(Wav1_part,floor(data.streams.eS1r.fs),floor(data.streams.Wav1.fs));

% input channel from presentation computer
Accl_beg_idx = (Wav1_beg_idx/data.streams.Wav1.fs)*data.streams.Accl.fs;
Accl_end_idx = (Wav1_end_idx/data.streams.Wav1.fs)*data.streams.Accl.fs;
Accl_part = double(data.streams.Accl.data(1,Accl_beg_idx:Accl_end_idx));
Accl_to_plot = resample(Accl_part,floor(data.streams.eS1r.fs),floor(data.streams.Accl.fs));

% e-stim channel
eS1r_beg_idx = (Wav1_beg_idx/data.streams.Wav1.fs)*data.streams.eS1r.fs;
eS1r_end_idx = (Wav1_end_idx/data.streams.Wav1.fs)*data.streams.eS1r.fs;

% filtered speech signal?
Filh_beg_idx = (Wav1_beg_idx/data.streams.Wav1.fs)*data.streams.Filh.fs;
Filh_end_idx = (Wav1_end_idx/data.streams.Wav1.fs)*data.streams.Filh.fs;
Filh_part = double(data.streams.Filh.data(1,Filh_beg_idx:Filh_end_idx));
Filh_to_plot = resample(Filh_part,floor(data.streams.eS1r.fs),floor(data.streams.Filh.fs));


figure('Units','normalized','Position', [0 0  1 1]);hold on

plot(220+50*Accl_to_plot)
plot(Filh_to_plot/50+100)
plot(Wav1_to_plot/6)
plot(data.streams.eS1r.data(1,eS1r_beg_idx:eS1r_end_idx)/40-75)

axis tight

legend('What goes into BNC (clean speech)','Wav1 channel','Filh channel','eS1r channel: electrical stimulation output')
title('Stimulation trial - Example2')
ylabel('(Amplitudes were changed for visualization)')
grid on

print(fullfile('/media/sakkol/HDD1/HBML/PROJECTS_DATA/Speech_Perception/NS148/iEEG_data/B65/PICS','example2.jpg'),'-djpeg','-r300')
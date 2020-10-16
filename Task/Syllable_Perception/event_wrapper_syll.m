function [events_table] = event_wrapper_syll(trial_count)

load('syll_table.mat','syll_table')

mixed_table = repmat(syll_table,[trial_count/32,1]);
rng('shuffle')
mixed_table = mixed_table(randperm(size(mixed_table,1)),:);
mixed_table.estim(1:trial_count/2)=1;
mixed_table.estim((trial_count/2)+1:end)=0;
mixed_table = mixed_table(randperm(size(mixed_table,1)),:);

events_table=table;
events_table.trials(1:trial_count)={''};
events_table.cfgs(1:trial_count)={''};
events_table.cond_info(1:trial_count)={''};

for t=1:trial_count
    
    cfg=[];
    cfg.language='English';
    cfg.SNR=20;
    cfg.LvsR='LcleanRnoise';
    cfg.noise='pink';
    cfg.part1.length=0.5;
    
    cfg.part2.chronicity='iso';
    cfg.part2.frequency=2;
    cfg.part2.word1=mixed_table.Stim{t};
    cfg.part2.estim=mixed_table.estim(t);
    cfg.part2.delay=0;
    
    cfg.part3.length=0.8;
    [events_table.trials{t}, events_table.cfgs{t}] = trial_creator(cfg);
    
    events_table.cond_info{t} = mixed_table(t,:);
end
end
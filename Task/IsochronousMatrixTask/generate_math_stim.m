function [math_stimuli,accuracy,operand] = generate_math_stim(trialN)
% Generating the math calculation stimuli for free recall version of
% isochronous matrix task.
%% Set several variables
percent_accr = 0.3;
percent_addition = 0.5;

%% Arrange numbers and addition/subtraction
math_stimuli = cell(trialN,10);

% based on percent accuracy
accuracy = zeros(trialN,2);
accuracy(randperm(trialN,floor(trialN*percent_accr)),1) = 1;
accuracy(randperm(trialN,floor(trialN*percent_accr)),2) = 1;

% which operand to use
operation_rand = zeros(trialN,2);
operation_rand(randperm(trialN,floor(trialN*percent_addition)),1) = 1;
operation_rand(randperm(trialN,floor(trialN*percent_addition)),2) = 1;
operand=cell(trialN,2);
operand(operation_rand(:,1)==0,1) = {'-'};
operand(operation_rand(:,1)==1,1) = {'+'};
operand(operation_rand(:,2)==0,2) = {'-'};
operand(operation_rand(:,2)==1,2) = {'+'};

% put them in output
math_stimuli(:,1) = cellstr(num2str(randi([10,90],trialN,1)));
math_stimuli(:,3) = cellstr(num2str(randi([1,9],trialN,1)));
math_stimuli(:,6) = cellstr(num2str(randi([10,90],trialN,1)));
math_stimuli(:,8) = cellstr(num2str(randi([1,9],trialN,1)));

math_stimuli(:,[4 9]) = {'='};

math_stimuli(:,2) = operand(:,1);
math_stimuli(:,7) = operand(:,2);

for t = 1:trialN
    
    if accuracy(t,1) == 1
        if strcmp(math_stimuli{t,2},'+')
            math_stimuli{t,5} = num2str(str2double(math_stimuli{t,1}) + str2double(math_stimuli{t,3}));
        elseif strcmp(math_stimuli{t,2},'-')
            math_stimuli{t,5} = num2str(str2double(math_stimuli{t,1}) - str2double(math_stimuli{t,3}));
        end
    elseif accuracy(t,1) == 0
        if strcmp(math_stimuli{t,2},'+')
            tmp = str2double(math_stimuli{t,1}) + str2double(math_stimuli{t,3});
        elseif strcmp(math_stimuli{t,2},'-')
            tmp = str2double(math_stimuli{t,1}) - str2double(math_stimuli{t,3});
        end
        
        tmp2 = randi([1 99]);
        while tmp == tmp2
            tmp2 = randi([1 99]);
        end
        
        math_stimuli{t,5} = num2str(tmp2);
    end
    
    if accuracy(t,2) == 1
        if strcmp(math_stimuli{t,7},'+')
            math_stimuli{t,10} = num2str(str2double(math_stimuli{t,6}) + str2double(math_stimuli{t,8}));
        elseif strcmp(math_stimuli{t,7},'-')
            math_stimuli{t,10} = num2str(str2double(math_stimuli{t,6}) - str2double(math_stimuli{t,8}));
        end
    elseif accuracy(t,2) == 0
        if strcmp(math_stimuli{t,7},'+')
            tmp = str2double(math_stimuli{t,6}) + str2double(math_stimuli{t,8});
        elseif strcmp(math_stimuli{t,7},'-')
            tmp = str2double(math_stimuli{t,6}) - str2double(math_stimuli{t,8});
        end
        
        tmp2 = randi([1 99]);
        while tmp == tmp2
            tmp2 = randi([1 99]);
        end
        
        math_stimuli{t,10} = num2str(tmp2);
    end
    
end


end

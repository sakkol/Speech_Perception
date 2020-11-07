function [math_stimuli,accuracy,first_dig] = generate_math_stim2(trialN)

%% Set several variables
percent_accr = 0.3;
percent_addition = 0.5;
math_mat = readtable('math_material.xlsx');
math_mat.Op1 = erase(cellstr(num2str(math_mat.Op1)),' ');
math_mat.Op3 = erase(cellstr(num2str(math_mat.Op3)),' ');
math_mat.Op4 = erase(cellstr(num2str(math_mat.Op4)),' ');

%% Arrange numbers and addition/subtraction
math_stimuli = cell(trialN,10);

% based on percent accuracy
accuracy = zeros(trialN,2);
accuracy(randperm(trialN,floor(trialN*percent_accr)),1) = 1;
accuracy(randperm(trialN,floor(trialN*percent_accr)),2) = 1;

% first single or first double digit
first_dig = ones(trialN,2);
first_dig(randperm(trialN,floor(trialN*percent_addition)),1) = 2;
first_dig(randperm(trialN,floor(trialN*percent_addition)),2) = 2;

% precreate and shuffle
acc0type1 = math_mat(math_mat.Acc == 0 & math_mat.Type == 1,:);
acc1type1 = math_mat(math_mat.Acc == 1 & math_mat.Type == 1,:);
acc0type2 = math_mat(math_mat.Acc == 0 & math_mat.Type == 2,:);
acc1type2 = math_mat(math_mat.Acc == 1 & math_mat.Type == 2,:);

acc0type1 = acc0type1(randperm(size(acc0type1,1)),:);
acc1type1 = acc1type1(randperm(size(acc1type1,1)),:);
acc0type2 = acc0type2(randperm(size(acc0type2,1)),:);
acc1type2 = acc1type2(randperm(size(acc1type2,1)),:);

acc0type1_num = 1;
acc1type1_num = 1;
acc0type2_num = 1;
acc1type2_num = 1;


% loop
for t = 1:trialN
    
    if accuracy(t,1) == 0
        if first_dig(t,1) == 1
            math_stimuli(t,1:5) = table2cell(acc0type1(acc0type1_num,1:5));acc0type1_num=acc0type1_num+1;
        elseif first_dig(t,1) == 2
            math_stimuli(t,1:5) = table2cell(acc0type2(acc0type2_num,1:5));acc0type2_num=acc0type2_num+1;
        end
    elseif accuracy(t,1) == 1
        if first_dig(t,1) == 1
            math_stimuli(t,1:5) = table2cell(acc1type1(acc1type1_num,1:5));acc1type1_num=acc1type1_num+1;
        elseif first_dig(t,1) == 2
            math_stimuli(t,1:5) = table2cell(acc1type2(acc1type2_num,1:5));acc1type2_num=acc1type2_num+1;
        end
    end
    
    % restart if more than present
    if acc0type1_num > size(acc0type1,1),acc0type1_num=1;end
    if acc0type2_num > size(acc0type2,1),acc0type2_num=1;end
    if acc1type1_num > size(acc1type1,1),acc1type1_num=1;end
    if acc1type2_num > size(acc1type2,1),acc1type2_num=1;end
    
    if accuracy(t,1) == 0
        if first_dig(t,1) == 1
            math_stimuli(t,6:10) = table2cell(acc0type1(acc0type1_num,1:5));acc0type1_num=acc0type1_num+1;
        elseif first_dig(t,1) == 2
            math_stimuli(t,6:10) = table2cell(acc0type2(acc0type2_num,1:5));acc0type2_num=acc0type2_num+1;
        end
    elseif accuracy(t,1) == 1
        if first_dig(t,1) == 1
            math_stimuli(t,6:10) = table2cell(acc1type1(acc1type1_num,1:5));acc1type1_num=acc1type1_num+1;
        elseif first_dig(t,1) == 2
            math_stimuli(t,6:10) = table2cell(acc1type2(acc1type2_num,1:5));acc1type2_num=acc1type2_num+1;
        end
    end
    
    % restart if more than present
    if acc0type1_num > size(acc0type1,1),acc0type1_num=1;end
    if acc0type2_num > size(acc0type2,1),acc0type2_num=1;end
    if acc1type1_num > size(acc1type1,1),acc1type1_num=1;end
    if acc1type2_num > size(acc1type2,1),acc1type2_num=1;end
    
end





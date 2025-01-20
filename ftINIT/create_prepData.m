% this script creates the prepData.mat file for the Human-GEM model

% load Human-GEM model
load('Human-GEM.mat');

% create prepData file 
% inputs: Human-GEM model, essential metabolic tasks & model reactions
prepData = prepHumanModelForftINIT(ihuman, false, 'path/to/HumanGEM/data/metabolicTasks/metabolicTasks_Essential.txt', 'path/to/HumanGEM/model/reactions.tsv');

% save prepData file
save('prepData.mat', 'prepData')
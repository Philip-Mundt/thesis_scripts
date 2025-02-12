% ****************************************************************************************************
% Step 3:   Pre-calculating essential reactions for Human1
% ****************************************************************************************************

% this script creates the prepData.mat file for the Human1 model

% inputs:
% path to Human-GEM folder
human-gem_path = 'path/to/Human-GEM/'; % do not forget the last "/"

% load Human-GEM model
load(strcat(human-gem_path, 'model/reactions.tsv'));

% create prepData file 
% inputs: Human-GEM model, essential metabolic tasks & model reactions
prepData = prepHumanModelForftINIT(ihuman, false, strcat(human-gem_path, 'data/metabolicTasks/metabolicTasks_Essential.txt'), strcat(human-gem_path, 'path/to/HumanGEM/model/reactions.tsv'));

% save prepData file
save('prepData.mat', 'prepData')
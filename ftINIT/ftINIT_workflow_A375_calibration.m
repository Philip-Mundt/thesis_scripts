% loading the original Human-GEM as well as a prepData file
load("Human-GEM.mat");
load("prepData.mat");


% Import the RNAseq data
expr_data = readtable("/vol/cluster-data/pmundt/A375_expression_data.txt");
[~, n] = size(expr_data);
numSamp = n-2; %the first two columns are the genes in ENSEMBL and gene symbols format

% Extract information from the table into a structure called data_struct
% gene names
data_struct.genes = expr_data{:, 1};
% sample (tissue) names
data_struct.tissues = expr_data.Properties.VariableNames(3:n);
% gene TPM values
data_struct.levels = expr_data{:, 3};

% Although the original tINIT implementation always compared gene 
% expression in the tissue of interest to the average of all other 
% provided tissues, ftINIT enables the option to instead compare the 
% expression to a threshold value. Here, we will use this alternative 
% approach with a threshold value of 1 TPM.
thresholds = {0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5};
thresholds_cptn = {'010', '025', '050', '075', '100', '125', '150', '175', '200', '250', '300', '400', '500'};

for i = 1:length(thresholds)
    data_struct.threshold = thresholds{i};


% run ftINIT
% the ftINIT optimization is divided into two steps to reduce computation time.
% ftINIT supports two modes: Mode “1+0” only runs the first step, 
% with the result that most reactions without Gene-Protein-Reaction (GPR) 
% associations are included in the final model, while the mode “1+1” 
% runs both steps. The “1+0” mode is suitable for structural comparisons, 
% while the “1+1” mode can be useful to generate a smaller model
    modes = {'1+0', '1+1'};

    for i_mode = 1:length(modes)

        model = ftINIT(prepData, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps(modes(i_mode)), false, true);
        save(strcat("/homes/pmundt/BA/models/A375_ftINIT_", modes(i_mode), "_thr", thresholds_cptn(i), "_model.mat"), "model")
    
    end

end
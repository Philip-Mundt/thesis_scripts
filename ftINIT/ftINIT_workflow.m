% this script creates an ftINIT from Human-GEM using gene expression data

% ****************************************************************************************************
% Step 4:   Loading reformatted gene expression dataset(s), Human1 model & prepData file 
%           and creating ftINIT models
% ****************************************************************************************************

% inputs:

% list of name(s) of your cell/tissue type(s)
% (expression data table(s) should be named "<cell/tissue type name>_expression_data.txt")
type_list = {'your_celltype', 'your_other_celltype', ...};

% path to folder containing reformatted expression dataset(s)
expr_data_path = 'path/to/expression_data_folder/';

% path to your output folder and cell/tissue type name
output_path = 'path/to/output_folder/';

% threshold (as a list) or multiple thresholds for TPM values 
% (genes with TPM values below this threshold are likely excluded from the model),
% default value: 1
thresholds = {0.5, 0.75, 1, 1.25, 1.5};
% captions for the thresholds as part of the name of the created model (floats will create an error)
thresholds_cptn = {'thr0050', 'thr0075', 'thr0100', 'thr0125', 'thr0150'};

% ftINIT modes (as a list)
% (Mode "1+0” only runs the first step, with the result that most reactions 
% without Gene-Protein-Reaction (GPR) associations are included in the 
% final model, while the mode "1+1” runs both steps. The “1+0” mode is 
% suitable for structural comparisons, while the “1+1” mode can be useful 
% to generate a smaller model)
modes = {'1+0', '1+1'};

# path to Human-GEM model (or other input model)
human-gem_path = 'path/to/Human-GEM.mat'

# path to prepData file created in step 3
prepData_path = 'path/to/prepData.m'


% ****************************************************************************************************
% start of the main part of the script

% loading the original Human-GEM as well as a prepData file
load(human-gem_path);
load(prepData_path);

for i_type = 1:length(type_list)

% Import the RNAseq data
    expr_data = readtable(strcat(expr_data_path, type_list(i_type), "_expression_data.txt"));
    [~, n] = size(expr_data);
    numSamp = n-2; %the first two columns are the genes in ENSEMBL and gene symbols format

% Extract information from the table into a structure called data_struct
% gene names
    data_struct.genes = expr_data{:, 1};
% sample (tissue) names
    data_struct.tissues = expr_data.Properties.VariableNames(3:n);
% gene TPM values
    data_struct.levels = expr_data{:, 3};

    for i_thr = 1:length(thresholds)

% load thresholds
        data_struct.threshold = thresholds{i_thr};

        for i_mode = 1:length(modes)
% run ftINIT & save model
            model = ftINIT(prepData, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps(modes(i_mode)), false, true);
            save(strcat(output_path, type_list(i_type), "_", modes(i_mode), "_", thresholds_cptn(i_thr), "_model.mat"), "model")
        
        end
        
    end

end
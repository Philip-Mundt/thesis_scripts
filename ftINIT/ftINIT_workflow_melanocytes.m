% this script creates an ftINIT from Human-GEM using gene expression data

% inputs:

% path to reformatted expression dataset
expr_data_path = "melanocyte_expression_data.txt";

% path to your output folder and cell/tissue type name
output_path = "/homes/pmundt/BA/models/";
type_name = "melanocyte";

% threshold (as a list) or multiple thresholds for TPM values 
% (genes with TPM values below this threshold are excluded from the model),
% default value: 1
thresholds = {1};
% captions for the name of the created model (floats will create an error)
thresholds_cptn = {0100};

% ftINIT mode (as a list) or both modes
% Mode "1+0” only runs the first step, with the result that most reactions 
% without Gene-Protein-Reaction (GPR) associations are included in the 
% final model, while the mode "1+1” runs both steps. The “1+0” mode is 
% suitable for structural comparisons, while the “1+1” mode can be useful 
% to generate a smaller model
modes = {'1+0', '1+1'};


% *************************************************************************
% start of the script

% loading the original Human-GEM as well as a prepData file
load("Human-GEM.mat");
load("prepData.mat");


% Import the RNAseq data
expr_data = readtable(expr_data_path);
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
for i = 1:length(thresholds)

    data_struct.threshold = thresholds{i};

% run ftINIT
% the ftINIT optimization is divided into two steps to reduce computation time.
    
    for i_mode = 1:length(modes)
        model = ftINIT(prepData, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps(modes(i_mode)), false, true);
        save(strcat(output_path, type_name, "_", modes(i_mode), "_thr", thresholds_cptn(i), "_model.mat"), "model")
    end

end
% loading the original Human-GEM as well as a prepData file
load("Human-GEM.mat");
load("prepData_gurobi.mat");

cancer_types = {
    'BRCA', 'BLCA', 'CESC', 'COAD', 'READ', 'GBM', 'HNSC', 'KIRC', ...
    'KIRP', 'KICH', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PRAD', ...
    'SKCM', 'STAD', 'TGCT', 'THCA', 'UCEC'
    };


for i = 1:length(cancer_types)
        
    % Import the RNAseq data
    expr_data = readtable(strcat("C:\Users\phili\OneDrive - uni-bielefeld.de\Uni\SoSe24\BA\test\cancer_expression_datasets\", cancer_types{i},"_expression_data.txt"));
    [~, n] = size(expr_data);
    numSamp = n-2; %the first two columns are the genes in ENSEMBL and gene symbols format
    
    % Extract information from the table into a structure called data_struct
    % gene names
    data_struct.genes = expr_data{:, 1};
    % sample (tissue) names
    data_struct.tissues = expr_data.Properties.VariableNames(3:n);
    % gene TPM values
    data_struct.levels = expr_data{:, 3:n};
    
    % Although the original tINIT implementation always compared gene 
    % expression in the tissue of interest to the average of all other 
    % provided tissues, ftINIT enables the option to instead compare the 
    % expression to a threshold value. Here, we will use this alternative 
    % approach with a threshold value of 1 TPM.
    data_struct.threshold = 1;
    
    
    % run ftINIT
    % the ftINIT optimization is divided into two steps to reduce computation time.
    % ftINIT supports two modes: Mode “1+0” only runs the first step, 
    % with the result that most reactions without Gene-Protein-Reaction (GPR) 
    % associations are included in the final model, while the mode “1+1” 
    % runs both steps. The “1+0” mode is suitable for structural comparisons, 
    % while the “1+1” mode can be useful to generate a smaller model
    model = ftINIT(prepdata_gurobi, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);
    save(strcat("C:\Users\phili\OneDrive - uni-bielefeld.de\Uni\SoSe24\BA\models\", cancer_types{i}, "_model.mat"), "model")

end
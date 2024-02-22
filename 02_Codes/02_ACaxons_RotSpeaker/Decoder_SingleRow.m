global PlotAnimalData
MainFolder = 'D:\AVspace_final';

timeVect = 0:1/6.0962:7;
base_window = [-1 0];
resp_window = [0.2  1.2];

metric = 'dFF';
selection_method = 'wilcoxon'; % 'ttest'; 'wilcoxon'; 'bootstrap'
alpha_threshold = 0.01;
dFF_threshold = 0.15;

%% load the data
animalID ={'CMad97','CMad98'};
n_animals = length(animalID);

    fprintf(1,'\nloading the data...')
    loadName = [MainFolder filesep '01_Data\03_ResponsiveData\ACaxons_SPKwn_2-80_bl6_dFF_wilcoxon_dot01_0_dot15.mat'];
    load(loadName)
    temp.data{1,1} = RespROIs.data{1,3};
    temp.data{1,2} = RespROIs.data{1,4};
    temp.data{1,3} = RespROIs.data{1,6};
    temp.data{1,4} = RespROIs.data{1,7};
    temp.data{2,1} = RespROIs.data{2,1};
    temp.data{2,2} = RespROIs.data{2,2};
    temp.data{2,3} = RespROIs.data{2,7};
    
    temp.pearsonR{1,1} = RespROIs.pearsonR{1,3};
    temp.pearsonR{1,2} = RespROIs.pearsonR{1,4};
    temp.pearsonR{1,3} = RespROIs.pearsonR{1,6};
    temp.pearsonR{1,4} = RespROIs.pearsonR{1,7};
    temp.pearsonR{2,1} = RespROIs.pearsonR{2,1};
    temp.pearsonR{2,2} = RespROIs.pearsonR{2,2};
    temp.pearsonR{2,3} = RespROIs.pearsonR{2,7};
        
    temp.info.data_details{1,1} = RespROIs.info.data_details{1,3};
    temp.info.data_details{1,2} = RespROIs.info.data_details{1,4};
    temp.info.data_details{1,3} = RespROIs.info.data_details{1,6};
    temp.info.data_details{1,4} = RespROIs.info.data_details{1,7};
    temp.info.data_details{2,1} = RespROIs.info.data_details{2,1};
    temp.info.data_details{2,2} = RespROIs.info.data_details{2,2};
    temp.info.data_details{2,3} = RespROIs.info.data_details{2,7};
    
%     loadName = [MainFolder filesep 'Data\ACaxons\RespData\ACaxons_SPKwn_2-80_cba_dFF_wilcoxon_dot01_0_dot15.mat'];
%     load(loadName)
%     
%     temp.data{3,1} = RespROIs.data{1,7};
%     temp.data{3,2} = RespROIs.data{1,8};
%     
%     temp.pearsonR{3,1} = RespROIs.pearsonR{1,7};
%     temp.pearsonR{3,2} = RespROIs.pearsonR{1,8};
%     
%     temp.info.data_details{3,1} = RespROIs.info.data_details{1,7};
%     temp.info.data_details{3,2} = RespROIs.info.data_details{1,8};
%     
%     temp.info.base_window = RespROIs.info.base_window;
%     temp.info.resp_window = RespROIs.info.resp_window;
%     
%     clear RespROIs;
%     RespROIs = temp;
%     
    RespROIs.info.animalID = animalID;
    
    
    fprintf(1,'done \n')

SaveFolder = [MainFolder filesep 'Plots_RotSPK']; if ~exist(SaveFolder,'dir');mkdir(SaveFolder);end

%%
Decoder_Output = BayesianDecoder_batch_SingleRow(RespROIs,SaveFolder,[-10:10:100]);

% function Batch_dataResponsive_v2(dataset,stim_toAnalyze)
%% Preambule
clearvars -except SMI_LMvis SMI_LMaud SMI_ACvis SMI_ACaud smi_LMvis smi_LMaud smi_ACaud Decoder_Output_AC Decoder_OutpuLMspkt Decoder_OutputLMled Decoder_Output_Bl6_80kHz Decoder_Output_CBA_80kHz Decoder_Output_RotSPK%RF_bl6_20kHz RF_bl6_80kHz RF_CBA RF_LMad RF_LMvis RF_bl6_vis
dataset = 'ACaxons';           % 'ACaxons', 'LMaxons','TEAaxons','V2Laxons'
stim_toAnalyze = 'SPKwn';       % 'SPKwn', 'LED'
soundStim = '2-80';            %'2-20', '2-80' max fqcy (in kHz). Only if ACaxons
bckgd = 'cba';                 % 'bl6', 'cba'. Only if 2-80 kHz sound


% ------------------------------------
% --------------- /!\ ----------------
removeFaceSession = true; % true means analysis on sessions w/o detectable face movement (see motion energy analysis)
% --------------- /!\ ----------------
% ------------------------------------


metric = 'dFF';               % 'dFF'; 'spike'
selection_method = 'wilcoxon';  % 'ttest'; 'wilcoxon'; 'bootstrap'
alpha_threshold = 0.01;
amp_threshold = 0;
dFF_threshold = 0.15;

% - nROIs responsives and nBoutons
basicQuantif = false;

% - are the notches in the response expected?
do_notches = false;

% - bayesian decoder
do_decoder = true;
doSingleSession = true; % some quirks with plotting
plotSingleSession = false;
doMouseagregate = false; % also does Mickey mouse
doAccVsV1 = false;

% - SMI
do_SMI = false;

% - ROI RF
do_RF = false;
Boutons = 1;          % Boutons = 1, boutons; 2, axons; else ROIs
rsquare_th = 0.3;     % what is a good fit?

% - Peak response to determine tuning
do_Peak = false;

% - ROI sensitive to az/ele?
doStimPosSensitive = false;

% - PCA
do_PCA = false;

% - noise correlations
do_noiseCorrelations = false;

% - Elevation
do_elevation = false;

% - Cross validate Boutons
do_XvalBoutons = false; % - will also do peakAz = f(V1RF)
n_resample = 100;
pearsonR_th = 0.3;

% - Grand Average
do_GrandAverages = false;
plotGrandAverage = false; % also for x-val boutons...
normalize = false;        % if false, then don't. if true, to edit coz it will look for x-val boutons in the sesction w/o x-val...

% - tuning curves. uses 'pearsonR_th' and 'do_XvalBoutons'
do_TuningCurves = false;
nShuffles = 100;       % H0: there is no tuning
magn_th   = 0;         % to eliminate ROIs with too small resp - maybe don't do it if X-Val
smi_range = [0 1];     % only plot tuning curve of ROIs within that range
XVal = 1;              % whether to cross validate to determine best position or not


%% edit what to do based on dependencies
if do_TuningCurves
    %     do_SMI = true;
    %     do_RF = true;
end

%% get what machine matlab is runing on
if isfolder('E:\AudSpace')
    MainFolder = 'E:\AudSpace';
    PCflag = 'openLab';
elseif isfolder('C:\Users\camil')
    MainFolder = 'C:\Users\camil\AudSpace';
end

%% load the data
switch dataset
    case 'ACaxons'
        switch soundStim
            case '2-20'
                animalID ={'CMad50','CMad54','CMad56',... % 2019
                    'CMad58','CMad62',...                 % Dec 2020
                    'CMad65','CMad67','CMad68'};          % Jan 2021
                TwentykHz = true;
            case '2-80'
                switch bckgd
                    case 'bl6'
                        animalID ={'CMad97','CMad98'};%,'CMad99'};
                    case 'cba'
                        animalID ={'CMad103','CMad104','CMad106',...
                            'CMad109','CMad110','CMad111'}; %CMad107
                end
                TwentykHz = false;
        end
        
    case 'LMaxons'
        animalID ={'CMad85','CMad86',...  % Jul 2021
            'RD10278',...                 % Jun 2021
            };
        TwentykHz = false;
        
    case 'TEAaxons'
        animalID ={'CMad119', 'CMad120'};
        TwentykHz = false;
        
    case 'V2Laxons'
        animalID ={'CMad85','CMad86',...  % Jul 2021
            'RD10278',...
            'CMad119', 'CMad120'};
        TwentykHz = false;
end

DataFolder = [MainFolder filesep 'Data'  filesep dataset];
if strcmp(dataset,'V2Laxons')
    DataFolder = [MainFolder filesep 'Data'  filesep 'LMaxons'];
end

dt = datestr(now,'yyyymmdd_HHMM');
if strcmp(stim_toAnalyze,'SPKwn')
    stimType = [stim_toAnalyze soundStim 'kHz' '_' bckgd];
else
    stimType = stim_toAnalyze;
end
temp = num2str(alpha_threshold); temp2 = strfind(temp,'.'); alpha = temp(temp2+1:end);
temp = num2str(dFF_threshold); temp2 = strfind(temp,'.'); dFF = temp(temp2+1:end);
SaveFolder = [MainFolder filesep 'Plots' filesep dataset filesep stimType filesep ...
    dt '_' selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF filesep];

if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end

if strcmp(soundStim,'2-80')
    dataname = [dataset '_' stim_toAnalyze '_' soundStim '_' bckgd '_' ...
        metric '_' selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF];
else
    try
        dataname = [dataset '_' stim_toAnalyze '_' ...
            metric '_' selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF];
    catch
        
        dataname = [dataset '_' stim_toAnalyze '_' ...
            selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF];
    end
end
fprintf(1,['\n loading ' dataname '...'])
load([DataFolder filesep 'RespData' filesep dataname '.mat'])
diaryName = 'AnalysisParameters';
diary([SaveFolder diaryName])
fprintf(1,'done \n')
disp(RespROIs.info);

RespROIs.info.stim_toAnalyze = stim_toAnalyze;

if removePupilSession
    doStimPosSensitive = false;
    animalID = {'CMad58','CMad62','CMad65','CMad67','CMad68'};
    for i = 1:5
        fprintf(['mouse' num2str(i)])
        counter = 0;
        for ii = 1:7
            if ~isempty(RespROIs.info.data_details{i+3,ii})
                if i == 2 && (ii == 1 || ii ==2 || ii==6)
                elseif i == 3 && ii == 3
                elseif i == 4 && (ii == 4 || ii == 5)
                else
                    fprintf([' pos' num2str(ii)])
                    counter = counter+1;
                    resp_window{i,counter} = RespROIs.info.resp_window{i+3,ii};
                    data_details{i,counter} = RespROIs.info.data_details{i+3,ii};
                    data{i,counter} = RespROIs.data{i+3,ii};
                    nBoutonsPerROI{i,counter} = RespROIs.nBoutonsPerROI{i+3,ii};
                    pearsonR{i,counter} = RespROIs.pearsonR{i+3,ii};
                    V1az(i,counter) = RespROIs.V1az(i+3,ii);
                    V1ele(i,counter) = RespROIs.V1ele(i+3,ii);
                    idxmax_MedianResp{i,counter} = RespROIs.idxmax_MedianResp{i+3,ii};
                end
            end
        end
        fprintf('\n')
    end
    RespROIs.info.resp_window = resp_window;
    RespROIs.info.animalID = animalID;
    RespROIs.info.data_details = data_details;
    RespROIs.data = data;
    RespROIs.nBoutonsPerROI = nBoutonsPerROI;
    RespROIs.pearsonR = pearsonR;
    RespROIs.V1az = V1az;
    RespROIs.V1ele = V1ele;
    RespROIs.idxmax_MedianResp = idxmax_MedianResp;
    
    SaveFolder = [MainFolder filesep 'Plots_ExludePupil' filesep dataset filesep stimType filesep ...
        dt '_' selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF filesep];
    
    if ~exist(SaveFolder,'dir')
        mkdir(SaveFolder)
    end
    
elseif removeFaceSession
    animalID = {'CMad103','CMad104','CMad97'};
%     animalID = {'CMad103','CMad104','CMad106','CMad109','CMad110','CMad111','CMad98'};
    n_mice = length(animalID);
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    
    for mouse = 1:n_mice-1   
        disp(animalID{mouse})
        counter = 0;
        for ii = 1:n_pos(mouse)
            if ~isempty(RespROIs.info.data_details{mouse,ii})
                if mouse == 1 && (ii == 2 || ii == 3 || ii == 4 || ii == 5 || ii == 6 || ii > 7)
                elseif mouse == 2 && (ii == 3 || ii > 4) % 
                else
                    fprintf([char(RespROIs.info.data_details{mouse,ii}) ' ,'])
                    counter = counter+1;
                    resp_window{mouse,counter} = RespROIs.info.resp_window{mouse,ii};
                    data_details{mouse,counter} = RespROIs.info.data_details{mouse,ii};
                    data{mouse,counter} = RespROIs.data{mouse,ii};
                    nBoutonsPerROI{mouse,counter} = RespROIs.nBoutonsPerROI{mouse,ii};
                    pearsonR{mouse,counter} = RespROIs.pearsonR{mouse,ii};
                    V1az(mouse,counter) = RespROIs.V1az(mouse,ii);
                    V1ele(mouse,counter) = RespROIs.V1ele(mouse,ii);
                    idxmax_MedianResp{mouse,counter} = RespROIs.idxmax_MedianResp{mouse,ii};
                end
% if mouse == 1
%     if (ii == 2 || ii == 3 || ii == 4 || ii == 5 || ii == 6 || ii > 7)
%     fprintf([char(RespROIs.info.data_details{mouse,ii}) ' ,'])
%     counter = counter+1;
%     resp_window{mouse,counter} = RespROIs.info.resp_window{mouse,ii};
%     data_details{mouse,counter} = RespROIs.info.data_details{mouse,ii};
%     data{mouse,counter} = RespROIs.data{mouse,ii};
%     nBoutonsPerROI{mouse,counter} = RespROIs.nBoutonsPerROI{mouse,ii};
%     pearsonR{mouse,counter} = RespROIs.pearsonR{mouse,ii};
%     V1az(mouse,counter) = RespROIs.V1az(mouse,ii);
%     V1ele(mouse,counter) = RespROIs.V1ele(mouse,ii);
%     idxmax_MedianResp{mouse,counter} = RespROIs.idxmax_MedianResp{mouse,ii};
%     end
% elseif mouse == 2 
%     if ii > 4
%     fprintf([char(RespROIs.info.data_details{mouse,ii}) ' ,'])
%     counter = counter+1;
%     resp_window{mouse,counter} = RespROIs.info.resp_window{mouse,ii};
%     data_details{mouse,counter} = RespROIs.info.data_details{mouse,ii};
%     data{mouse,counter} = RespROIs.data{mouse,ii};
%     nBoutonsPerROI{mouse,counter} = RespROIs.nBoutonsPerROI{mouse,ii};
%     pearsonR{mouse,counter} = RespROIs.pearsonR{mouse,ii};
%     V1az(mouse,counter) = RespROIs.V1az(mouse,ii);
%     V1ele(mouse,counter) = RespROIs.V1ele(mouse,ii);
%     idxmax_MedianResp{mouse,counter} = RespROIs.idxmax_MedianResp{mouse,ii};
%     end
% else
%     fprintf([char(RespROIs.info.data_details{mouse,ii}) ' ,'])
%     counter = counter+1;
%     resp_window{mouse,counter} = RespROIs.info.resp_window{mouse,ii};
%     data_details{mouse,counter} = RespROIs.info.data_details{mouse,ii};
%     data{mouse,counter} = RespROIs.data{mouse,ii};
%     nBoutonsPerROI{mouse,counter} = RespROIs.nBoutonsPerROI{mouse,ii};
%     pearsonR{mouse,counter} = RespROIs.pearsonR{mouse,ii};
%     V1az(mouse,counter) = RespROIs.V1az(mouse,ii);
%     V1ele(mouse,counter) = RespROIs.V1ele(mouse,ii);
%     idxmax_MedianResp{mouse,counter} = RespROIs.idxmax_MedianResp{mouse,ii};
% end
            end
        end
        fprintf('\n')
    end
    
    bckgd = 'bl6';
    dataname = [dataset '_' stim_toAnalyze '_' soundStim '_' bckgd '_' ...
        metric '_' selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF];
    load([DataFolder filesep 'RespData' filesep dataname '.mat'])
    
    disp(animalID{n_mice})
    counter = 1;
    for ii = 1:8
        fprintf([char(RespROIs.info.data_details{1,ii})  ' ,'])
        resp_window{n_mice,counter} = RespROIs.info.resp_window{1,ii};
        data_details{n_mice,counter} = RespROIs.info.data_details{1,ii};
        data{n_mice,counter} = RespROIs.data{1,ii};
        nBoutonsPerROI{n_mice,counter} = RespROIs.nBoutonsPerROI{1,ii};
        pearsonR{n_mice,counter} = RespROIs.pearsonR{1,ii};
        V1az(n_mice,counter) = RespROIs.V1az(1,ii);
        V1ele(n_mice,counter) = RespROIs.V1ele(1,ii);
        idxmax_MedianResp{n_mice,counter} = RespROIs.idxmax_MedianResp{1,ii};
        counter = counter+1;
    end
    
%   for ii = 1:8
%         fprintf([char(RespROIs.info.data_details{2,ii})  ' ,'])
%         resp_window{n_mice,counter} = RespROIs.info.resp_window{2,ii};
%         data_details{n_mice,counter} = RespROIs.info.data_details{2,ii};
%         data{n_mice,counter} = RespROIs.data{2,ii};
%         nBoutonsPerROI{n_mice,counter} = RespROIs.nBoutonsPerROI{2,ii};
%         pearsonR{n_mice,counter} = RespROIs.pearsonR{2,ii};
%         V1az(n_mice,counter) = RespROIs.V1az(2,ii);
%         V1ele(n_mice,counter) = RespROIs.V1ele(2,ii);
%         idxmax_MedianResp{n_mice,counter} = RespROIs.idxmax_MedianResp{2,ii};
%         counter = counter+1;
%   end
    fprintf('\n')
    
    RespROIs.info.resp_window = resp_window;
    RespROIs.info.animalID = animalID;
    RespROIs.info.data_details = data_details;
    RespROIs.data = data;
    RespROIs.nBoutonsPerROI = nBoutonsPerROI;
    RespROIs.pearsonR = pearsonR;
    RespROIs.V1az = V1az;
    RespROIs.V1ele = V1ele;
    RespROIs.idxmax_MedianResp = idxmax_MedianResp;
    
    SaveFolder = [MainFolder filesep 'Plots_ExludePupil' filesep dataset filesep stimType filesep ...
        dt '_' selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF filesep];
    
    if ~exist(SaveFolder,'dir')
        mkdir(SaveFolder)
    end
    
    
end

%% some common variable definition
n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
timeVect = 0:1/6.0962:7;

%% quantif nBoutons and nROIs
if basicQuantif
    n_animals = size(RespROIs.info.animalID ,2);
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    colors = distinguishable_colors(n_animals);
    
    % - fraction responsive
    figure; subplot(1,3,[1 2]);hold on
    for i = 1:n_animals
        frResp = RespROIs.fractionResp.nResp(i,~isnan(RespROIs.fractionResp.nResp(i,:)))./RespROIs.fractionResp.nROIs(i,~isnan(RespROIs.fractionResp.nROIs(i,:)));
        %         frResp = RespROIs.fractionResp.nResp(i,RespROIs.fractionResp.nROIs(i,:)~=0)./RespROIs.fractionResp.nROIs(i,RespROIs.fractionResp.nROIs(i,:)~=0);
        frResp(isnan(frResp))=[];
        scatter(i*ones(n_pos(i),1),frResp,...
            'filled','markerfacecolor',colors(i,:),'markerfacealpha',0.5)
        sh_mean = sum(RespROIs.fractionResp.nResp_sh(i,~isnan(RespROIs.fractionResp.nResp_sh(i,:))))./sum(RespROIs.fractionResp.nROIs(i,~isnan(RespROIs.fractionResp.nROIs(i,:))));
        %         sh_mean = sum(RespROIs.fractionResp.nResp_sh(i,RespROIs.fractionResp.nROIs(i,:)~=0))./sum(RespROIs.fractionResp.nROIs(i,RespROIs.fractionResp.nROIs(i,:)~=0));
        plot([i-.2 i+.2],[sh_mean sh_mean],'k:')
    end
    temp  = RespROIs.fractionResp.nResp./RespROIs.fractionResp.nROIs;
    allSessions = reshape(temp,size(temp,1)*size(temp,2),1);
    sem = std(allSessions,'omitnan')./sqrt(length(allSessions)-sum(isnan(allSessions)));
    plot([n_animals+1-0.2  n_animals+1.2],[nanmean(allSessions) nanmean(allSessions)],...
        'k','linewidth',3)
    plot([n_animals+1 n_animals+1],[nanmean(allSessions)-sem nanmean(allSessions)+sem],'k')
    ylabel('fraction responsive');
    xlabel('animal ID'); xlim([0 n_animals+2]); xticks([1:n_animals+1]); xticklabels(RespROIs.info.animalID);xtickangle(45)
    
    subplot(1,3,3);hold on
    allSessions(isnan(allSessions)) = [];
    scatter(0.8*ones(size(allSessions)),allSessions,'filled','markerfacecolor',[0 0 0],'markerfacealpha',0.2)
    plot([1 1.2],[mean(allSessions) mean(allSessions)],'k')
    sem = std(allSessions)./sqrt(length(allSessions));
    plot([1.1 1.1],[mean(allSessions)-sem mean(allSessions)+sem],'k')
    temp  = RespROIs.fractionResp.nResp_sh./RespROIs.fractionResp.nROIs;
    allSessions_sh = reshape(temp,size(temp,1)*size(temp,2),1); allSessions_sh(isnan(allSessions_sh)) = [];
    plot([0.8 1.2],[mean(allSessions_sh) mean(allSessions_sh)],'k:')
    
    allMice = sum(RespROIs.fractionResp.nResp,2,'omitnan')./sum(RespROIs.fractionResp.nROIs,2,'omitnan');
    scatter(1.8*ones(size(allMice)),allMice,'filled','markerfacecolor',[0 0 0],'markerfacealpha',0.2)
    plot([2 2.2],[mean(allMice) mean(allMice)],'k')
    sem = std(allMice)./sqrt(length(allMice));
    plot([2.1 2.1],[mean(allMice)-sem mean(allMice)+sem],'k')
    allMice_sh  = sum(RespROIs.fractionResp.nResp_sh,2,'omitnan')./sum(RespROIs.fractionResp.nROIs,2,'omitnan');
    plot([1.8 2.2],[mean(allMice_sh) mean(allMice_sh)],'k:')
    
    temp = sum(RespROIs.fractionResp.nResp,'all','omitnan')/sum(RespROIs.fractionResp.nROIs,'all','omitnan');
    plot([2.8 3.2],[temp temp],'k')
    temp = sum(RespROIs.fractionResp.nResp_sh,'all','omitnan')/sum(RespROIs.fractionResp.nROIs,'all','omitnan');
    plot([2.8 3.2],[temp temp],'k:')
    
    xticks([1:3]);xticklabels({'across sessions','across mice','across all'});xlim([0 4]);xtickangle(45)
    
    SaveFolder_FrResp = [SaveFolder filesep 'FrResp' filesep];
    if ~exist(SaveFolder_FrResp,'dir');mkdir(SaveFolder_FrResp);end
    saveas(gcf,[SaveFolder_FrResp filesep 'FrResp.tif'])
    savefig([SaveFolder_FrResp filesep 'FrResp'])
    
    
    % nBoutons
    nBoutons = 0; nROIs = 0;
    for i = 1:n_animals
        for ii = 1:n_pos(i)
            temp = sum(RespROIs.nBoutonsPerROI{i,ii});
            nBoutons = nBoutons+temp;
            temp = length(RespROIs.nBoutonsPerROI{i,ii});
            nROIs = nROIs+temp;
        end
    end
    disp([num2str(nBoutons) ' boutons from ' num2str(nROIs) ' ROIs'])
end

%% n az sensitive, ele sensitive boutons
% break ROIs down into boutons
if doStimPosSensitive
    fprintf('analyzing the space-dependency in the response...')
    SaveFolderSpaceSensitive = [SaveFolder filesep 'SpaceSensitive' filesep];
    if ~exist(SaveFolderSpaceSensitive,'dir')
        mkdir(SaveFolderSpaceSensitive)
    end
    IsMod = StimPosSensitive(RespROIs,SaveFolderSpaceSensitive);
    fprintf('done\n')
end

if do_elevation
    SaveFolderEleSpaceSensitive = [SaveFolder filesep 'SpaceSensitive/Elevation' filesep];
    if ~exist(SaveFolderEleSpaceSensitive,'dir')
        mkdir(SaveFolderEleSpaceSensitive)
    end
    ElePosSensitive(RespROIs,SaveFolderEleSpaceSensitive);
    
end

%% PCA
if do_PCA
    SaveFolderPCA = [SaveFolder 'PCA' filesep];
    if ~exist(SaveFolderPCA,'dir')
        mkdir(SaveFolderPCA)
    end
    
    cumvar_all = performPCA(RespROIs,SaveFolderPCA);
    
end

%% -- SMI --
if do_SMI
    % is working with new data structure and splits ROIs (1Sep2021)
    SaveFolderSMI = [SaveFolder filesep 'SMI' filesep];
    if ~exist(SaveFolderSMI,'dir');mkdir(SaveFolderSMI);end
    SMI = SMI_v2(RespROIs,SaveFolderSMI);
end



%% (normalized) grand average
if plotGrandAverage
    SaveFolder_GA = [SaveFolder filesep 'GrandAverage_not_x-val' filesep];
    if ~exist(SaveFolder_GA,'dir');mkdir(SaveFolder_GA);end
    
    dff_all = cell(n_animals,max(n_pos));
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            data = RespROIs.data{animal,pos};
            nROIs = size(data,1);
            
            % % - analysis time windows
            try
                resp_window = RespROIs.info.resp_window{animal,pos};
            catch
                resp_window = RespROIs.info.resp_window;
            end
            base_window = RespROIs.info.base_window;
            tAna = false(size(timeVect));
            tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
            tAna = tAna(1:size(data,2));
            tBase = false(size(timeVect));
            tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
            tBase = tBase(1:size(data,2));
            
            % % - do a dFF
            F0 = nanmean(data(:,tBase,:,:),2);
            ROI_wNegBase = 0;
            for ROI = 1:nROIs
                Fo=squeeze(F0(ROI,:,:,:));
                if any(Fo<0,'all')
                    ROI_wNegBase= ROI_wNegBase+1;
                    minF = min(data(ROI,:,:,:),[],'all');
                    data(ROI,:,:,:) = data(ROI,:,:,:)-minF;
                    F0(ROI,1,:,:) = nanmean(data(ROI,tBase,:,:),2);
                end
            end
            dFF = (data-F0)./F0;
            
            % - transform ROIs into boutons
            nBoutons_max = max(RespROIs.nBoutonsPerROI{animal,pos});
            for BoutonPerROI = 2:nBoutons_max
                multipleBoutons = RespROIs.nBoutonsPerROI{animal,pos}==BoutonPerROI;
                dFF = cat(1,dFF,repmat(dFF(multipleBoutons,:,:,:),BoutonPerROI-1,1,1,1));
            end
            
            
            dff_all{animal,pos} = mean(dFF(:,1:36,:,:),4);
        end
    end
    
    if normalize
        trialAvgData = cat(1,dff_all{:});
        respAvgData = squeeze(mean(trialAvgData(:,tAna,:),2));
        normTrialAvgData = trialAvgData./max(respAvgData,[],2);
        grandAverage = NaN(36,39,n_resample);
        for Rand = 1:n_resample
            selected = logical(cat(2,corrBoutons{:,:,Rand}));
            nROIs(Rand) = sum(selected);
            grandAverage(:,:,Rand) = squeeze(mean(normTrialAvgData(selected,:,:)));
        end
        % % - plot the normalized grand average
        grandAverage_mean = mean(grandAverage,3);
        grandAverage_CI = prctile(grandAverage,[5 95],3);
        ymin = min(grandAverage_CI(:,:,2),[],'all');
        ymax = max(grandAverage_CI(:,:,1),[],'all');
        ymax = ymax+0.1;
        TraceFig = figure;
        for i = 1:39
            figure(TraceFig)
            subtightplot(3,13,i); hold on
            plot([1 36],[0 0],'k:')
            fill([6 18 18 6 6],[ymin ymin ymax ymax ymin]',...
                'g','edgecolor','none','facealpha',.2)
            plot(grandAverage_mean(:,i),'k')
            fill([1:36 36:-1:1],[grandAverage_CI(:,i,1);flip(grandAverage_CI(:,i,2))]',...
                'k','edgecolor','none','facealpha',.2)
            ylim([ymin ymax]);
            if i~=27
                axis off
            end
        end
        HeatMapFig = figure;
        grandAverage_RespMean = mean(grandAverage_mean(tAna,:));
        grandAverage_RespMean_2D = reshape(grandAverage_RespMean,13,3);
        %     imagesc(grandAverage_RespMean_2D',[0 max(grandAverage_RespMean_2D,[],'all')])
        imagesc(grandAverage_RespMean_2D',[0 .3]);
        colormap gray
        c=colorbar;
        c.Ticks = [0 .3];
        
    else
        % % - not normalized
        
        % % % - all boutons - % % %
        %         grandAverage_all = squeeze(mean(trialAvgData));
        %         grandAverage_all_sem = squeeze(std(trialAvgData)./sqrt(size(trialAvgData,1)));
        %              ymin = min(grandAverage_all-grandAverage_all_sem,[],'all');
        %             ymax = max(grandAverage_all+grandAverage_all_sem,[],'all');
        %             TraceFig2 = figure;
        %             for i = 1:39
        %                 figure(TraceFig2)
        %                 subtightplot(3,13,i); hold on
        %                 plot([1 36],[0 0],'k:')
        %                 fill([6 18 18 6 6],[ymin ymin ymax ymax ymin]',...
        %                     'g','edgecolor','none','facealpha',.2)
        %                 plot(grandAverage_all(:,i),'k')
        %                 fill([1:36 36:-1:1],[grandAverage_all(:,i)-grandAverage_all_sem(:,i);...
        %                     flip(grandAverage_all(:,i)+grandAverage_all_sem(:,i))]',...
        %                     'k','edgecolor','none','facealpha',.2)
        %                 ylim([ymin ymax]);
        %                 if i~=27
        %                     axis off
        %                 end
        %             end
        %             HeatMapFig2 = figure;
        %             grandAverage_RespMean = mean(grandAverage_all(tAna,:));
        %             grandAverage_RespMean_2D = reshape(grandAverage_RespMean,13,3);
        %             %         imagesc(grandAverage_RespMean_2D',[0 max(grandAverage_RespMean_2D,[],'all')])
        %             imagesc(grandAverage_RespMean_2D',[0 .5]);
        %             colormap gray
        %             c=colorbar;
        %             c.Ticks = [0 .5];
        % % % ------------- % % %
        
        % % % - boutons per session - % % %
        for animal = 1:n_animals
            for pos = 1:n_pos(animal)
                grandAverage_all = squeeze(mean(dff_all{animal,pos}));
                grandAverage_all_sem = squeeze(std(dff_all{animal,pos})./sqrt(size(dff_all,1)));
                ymin = min(grandAverage_all-grandAverage_all_sem,[],'all');
                ymax = max(grandAverage_all+grandAverage_all_sem,[],'all');
                ymax = floor(ymax*10); ymax = ymax/10;
                
                TraceFig2 = figure;
                toPlot = movmean(grandAverage_all,3,1);
                for i = 1:39
                    figure(TraceFig2)
                    subtightplot(3,13,i); hold on
                    plot([1 36],[0 0],'k:')
                    %                         if strcmpi(stim_toAnalyze,'LED')
                    fill([24 30 30 24 24],[ymin ymin ymax ymax ymin]',...
                        'r','edgecolor','none','facealpha',.2)
                    fill([6 12 12 6 6],[ymin ymin ymax ymax ymin]',...
                        'g','edgecolor','none','facealpha',.2)
                    %                         elseif strcmpi(stim_toAnalyze,'SPKwn')
                    %                            fill([6 12 12 6 6],[ymin ymin ymax ymax ymin]',...
                    %                             'g','edgecolor','none','facealpha',.2)
                    %                         end
                    plot(toPlot(:,i),'k')
                    fill([1:36 36:-1:1],[toPlot(:,i)-grandAverage_all_sem(:,i);...
                        flip(toPlot(:,i)+grandAverage_all_sem(:,i))]',...
                        'k','edgecolor','none','facealpha',.2)
                    ylim([ymin ymax]);
                    if i~=27
                        axis off
                    end
                end
                ze_title = [RespROIs.info.animalID{animal} ', pos' num2str(pos)];
                sgtitle(ze_title)
                saveas(gcf,[SaveFolder_GA ze_title '_NonNorm_NonX-Val_Trace.tif'])
                saveas(gcf,[SaveFolder_GA ze_title '_NonNorm_NonX-Val_Trace.svg'])
                savefig([SaveFolder_GA ze_title '_NonNorm_NonX-Val_Trace'])
                close
                
                HeatMapFig2 = figure;
                grandAverage_RespMean = mean(grandAverage_all(tAna,:));
                grandAverage_RespMean_2D = reshape(grandAverage_RespMean,13,3);
                ymax = max(grandAverage_RespMean);
                ymax = floor(ymax*10); ymax = ymax/10;
                %         imagesc(grandAverage_RespMean_2D',[0 max(grandAverage_RespMean_2D,[],'all')])
                imagesc(grandAverage_RespMean_2D',[0 ymax]);
                colormap gray
                c=colorbar;
                c.Ticks = [0 ymax];
                sgtitle(ze_title)
                saveas(gcf,[SaveFolder_GA ze_title '_NonNorm_NonX-Val_HM.tif'])
                saveas(gcf,[SaveFolder_GA ze_title '_NonNorm_NonX-Val_HM.svg'])
                savefig([SaveFolder_GA ze_title '_NonNorm_NonX-Val_HM'])
                close
            end
        end
        % % % ------------- % % %
        
    end
end % end if plot grand averages

%% cross-validate ROIs
% % - Peak az and ele of cross-validated ROIs
if do_XvalBoutons
    n_animals = size(RespROIs.info.animalID ,2);
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    timeVect = 0:1/6.0962:7;
    
    corrBoutons = cell(n_animals,max(n_pos),n_resample);
    
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            
            data = RespROIs.data{animal,pos};
            nROIs = size(data,1);
            
            % % - analysis time windows
            try
                resp_window = RespROIs.info.resp_window{animal,pos};
            catch
                resp_window = RespROIs.info.resp_window;
            end
            base_window = RespROIs.info.base_window;
            tAna = false(size(timeVect));
            tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
            tAna = tAna(1:size(data,2));
            tBase = false(size(timeVect));
            tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
            tBase = tBase(1:size(data,2));
            
            % % - do a dFF
            F0 = nanmean(data(:,tBase,:,:),2);
            ROI_wNegBase = 0;
            for ROI = 1:nROIs
                Fo=squeeze(F0(ROI,:,:,:));
                if any(Fo<0,'all')
                    ROI_wNegBase= ROI_wNegBase+1;
                    minF = min(data(ROI,:,:,:),[],'all');
                    data(ROI,:,:,:) = data(ROI,:,:,:)-minF;
                    F0(ROI,1,:,:) = nanmean(data(ROI,tBase,:,:),2);
                end
            end
            dFF = (data-F0)./F0;
            
            % - transform ROIs into boutons
            nBoutons_max = max(RespROIs.nBoutonsPerROI{animal,pos});
            for BoutonPerROI = 2:nBoutons_max
                multipleBoutons = RespROIs.nBoutonsPerROI{animal,pos}==BoutonPerROI;
                dFF = cat(1,dFF,repmat(dFF(multipleBoutons,:,:,:),BoutonPerROI-1,1,1,1));
            end
            
            dFF_all{animal,pos} = dFF;
            
            % % - cross-correlation
            nRep = size(dFF,4);
            for Rand = 1:n_resample
                % % - template
                [template_data,idx] = datasample(dFF,nRep/2,4,'Replace',false);
                template_mean = squeeze(mean(mean(template_data(:,tAna,:,:),2),4));
                [max_template_mean,idx_template] = max(template_mean,[],2);
                
                % % - test
                test_idx = true(nRep,1);
                test_idx(idx) = false;
                test_data = dFF(:,:,:,test_idx);
                test_mean = squeeze(mean(mean(test_data(:,tAna,:,:),2),4));
                [~,idx_test] = max(test_mean,[],2);
                
                % % - correlation
                nROIs = size(template_mean,1);
                R = NaN(1,nROIs);
                for i = 1:nROIs
                    Rtemp = corrcoef(template_mean(i,:),test_mean(i,:));
                    R(i) = Rtemp(1,2);
                end
                correlatedBoutons = false(1,nROIs);
                correlatedBoutons(R>pearsonR_th) = true;
                
                corrBoutons{animal,pos,Rand} = correlatedBoutons;
            end % end loop through random subsampling
        end % end loop through positions
    end % end loop through animals
    
    %% (normalized) grand average
    if plotGrandAverage
        temp = cell(n_animals,max(n_pos));
        for animal = 1:n_animals
            for pos = 1:n_pos(animal)
                temp{animal,pos} = mean(dFF_all{animal,pos}(:,1:36,:,:),4);
            end
        end
        trialAvgData = cat(1,temp{:});
        %        data_all = cat(1,dFF_all{:});
        if normalize
            respAvgData = squeeze(mean(trialAvgData(:,tAna,:),2));
            normTrialAvgData = trialAvgData./max(respAvgData,[],2);
            grandAverage = NaN(36,39,n_resample);
            for Rand = 1:n_resample
                selected = logical(cat(2,corrBoutons{:,:,Rand}));
                nROIs(Rand) = sum(selected);
                grandAverage(:,:,Rand) = squeeze(mean(normTrialAvgData(selected,:,:)));
            end
            % % - plot the normalized grand average
            grandAverage_mean = mean(grandAverage,3);
            grandAverage_CI = prctile(grandAverage,[5 95],3);
            ymin = min(grandAverage_CI(:,:,2),[],'all');
            ymax = max(grandAverage_CI(:,:,1),[],'all');
            ymax = ymax+0.1;
            TraceFig = figure;
            for i = 1:39
                figure(TraceFig)
                subtightplot(3,13,i); hold on
                plot([1 36],[0 0],'k:')
                fill([6 18 18 6 6],[ymin ymin ymax ymax ymin]',...
                    'g','edgecolor','none','facealpha',.2)
                plot(grandAverage_mean(:,i),'k')
                fill([1:36 36:-1:1],[grandAverage_CI(:,i,1);flip(grandAverage_CI(:,i,2))]',...
                    'k','edgecolor','none','facealpha',.2)
                ylim([ymin ymax]);
                if i~=27
                    axis off
                end
            end
            HeatMapFig = figure;
            grandAverage_RespMean = mean(grandAverage_mean(tAna,:));
            grandAverage_RespMean_2D = reshape(grandAverage_RespMean,13,3);
            %     imagesc(grandAverage_RespMean_2D',[0 max(grandAverage_RespMean_2D,[],'all')])
            imagesc(grandAverage_RespMean_2D',[0 .3]);
            colormap gray
            c=colorbar;
            c.Ticks = [0 .3];
            
        else
            % % - not normalized
            
            % % % - all boutons - % % %
            %         grandAverage_all = squeeze(mean(trialAvgData));
            %         grandAverage_all_sem = squeeze(std(trialAvgData)./sqrt(size(trialAvgData,1)));
            %              ymin = min(grandAverage_all-grandAverage_all_sem,[],'all');
            %             ymax = max(grandAverage_all+grandAverage_all_sem,[],'all');
            %             TraceFig2 = figure;
            %             for i = 1:39
            %                 figure(TraceFig2)
            %                 subtightplot(3,13,i); hold on
            %                 plot([1 36],[0 0],'k:')
            %                 fill([6 18 18 6 6],[ymin ymin ymax ymax ymin]',...
            %                     'g','edgecolor','none','facealpha',.2)
            %                 plot(grandAverage_all(:,i),'k')
            %                 fill([1:36 36:-1:1],[grandAverage_all(:,i)-grandAverage_all_sem(:,i);...
            %                     flip(grandAverage_all(:,i)+grandAverage_all_sem(:,i))]',...
            %                     'k','edgecolor','none','facealpha',.2)
            %                 ylim([ymin ymax]);
            %                 if i~=27
            %                     axis off
            %                 end
            %             end
            %             HeatMapFig2 = figure;
            %             grandAverage_RespMean = mean(grandAverage_all(tAna,:));
            %             grandAverage_RespMean_2D = reshape(grandAverage_RespMean,13,3);
            %             %         imagesc(grandAverage_RespMean_2D',[0 max(grandAverage_RespMean_2D,[],'all')])
            %             imagesc(grandAverage_RespMean_2D',[0 .5]);
            %             colormap gray
            %             c=colorbar;
            %             c.Ticks = [0 .5];
            % % % ------------- % % %
            
            % % % - boutons per session - % % %
            for animal = 1:n_animals
                for pos = 1:n_pos(animal)
                    grandAverage_all = squeeze(mean(temp{animal,pos}));
                    grandAverage_all_sem = squeeze(std(trialAvgData)./sqrt(size(trialAvgData,1)));
                    ymin = min(grandAverage_all-grandAverage_all_sem,[],'all');
                    ymax = max(grandAverage_all+grandAverage_all_sem,[],'all');
                    floor(ymax*10)
                    TraceFig2 = figure;
                    for i = 1:39
                        figure(TraceFig2)
                        subtightplot(3,13,i); hold on
                        plot([1 36],[0 0],'k:')
                        fill([6 18 18 6 6],[ymin ymin ymax ymax ymin]',...
                            'g','edgecolor','none','facealpha',.2)
                        plot(grandAverage_all(:,i),'k')
                        fill([1:36 36:-1:1],[grandAverage_all(:,i)-grandAverage_all_sem(:,i);...
                            flip(grandAverage_all(:,i)+grandAverage_all_sem(:,i))]',...
                            'k','edgecolor','none','facealpha',.2)
                        ylim([ymin ymax]);
                        if i~=27
                            axis off
                        end
                    end
                    HeatMapFig2 = figure;
                    grandAverage_RespMean = mean(grandAverage_all(tAna,:));
                    grandAverage_RespMean_2D = reshape(grandAverage_RespMean,13,3);
                    %         imagesc(grandAverage_RespMean_2D',[0 max(grandAverage_RespMean_2D,[],'all')])
                    imagesc(grandAverage_RespMean_2D',[0 .5]);
                    colormap gray
                    c=colorbar;
                    c.Ticks = [0 .5];
                    pause
                    title([mouse{animal}])
                end
            end
            % % % ------------- % % %
            
        end
    end % end if plot grand averages
    
    %% peakAz = f(V1RF)
    scatterFig = figure; hold on
    frContraFig = figure; hold on
    meanPeakAz_mean = NaN(n_animals,max(n_pos));
    meanPeakAz_CI = NaN(n_animals,max(n_pos),2);
    meanPeakAz = NaN(n_animals,max(n_pos),n_resample);
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            N = NaN(13,n_resample);
            for Rand = 1:n_resample
                data = dFF_all{animal,pos};
                selected = corrBoutons{animal,pos,Rand};
                %                 smi = SMI{animal,pos};
                %                 smi_selec = false(1,size(smi,1));
                %                 smi_selec(smi>0.3) = true;
                selected2 = (selected);%&smi_selec);
                if sum(selected2)>10
                    AvgResp = squeeze(mean(mean(data(selected2,tAna,:,:),4),2));
                    AvgAzResp = reshape(AvgResp,size(AvgResp,1),13,3);
                    AvgAzResp = mean(AvgAzResp,3);
                    [~,peakResp] = max(AvgAzResp,[],2);
                    meanPeakAz(animal,pos,Rand) = median(peakResp);
                    
                    N(:,Rand) = histcounts(peakResp,[1:14],'Normalization','probability');
                    %                     if Rand==1
                    %                         jitter = rand(length(peakResp),1)*5-2.5;
                    %                         xData = repmat(RespROIs.V1az(animal,pos),length(peakResp),1);
                    %                         scatter(xData+jitter,peakResp,'k','filled','markerfacealpha',.2)
                    %                     end
                end
                %                 if animal == 8 && pos == 2
                %                     figure;histogram(peakResp);
                %                     keyboard
                %                 end
            end
            if animal == 8 && pos == 2
                keyboard
            end
            if  isnan(squeeze(meanPeakAz(animal,pos,:)))<n_resample
                figure(scatterFig)
                meanPeakAz_mean(animal,pos) = median(meanPeakAz(animal,pos,:),3,'omitnan');
                meanPeakAz_CI(animal,pos,:) = prctile(meanPeakAz(animal,pos,:),[5 95],3);
                if RespROIs.V1az(animal,pos)>100
                    RespROIs.V1az(animal,pos) = 100;
                end
                fill([RespROIs.V1az(animal,pos)-1 RespROIs.V1az(animal,pos)+1 RespROIs.V1az(animal,pos)+1 RespROIs.V1az(animal,pos)-1],...
                    [squeeze(meanPeakAz_CI(animal,pos,1)) squeeze(meanPeakAz_CI(animal,pos,1)) squeeze(meanPeakAz_CI(animal,pos,2)) squeeze(meanPeakAz_CI(animal,pos,2))],...
                    'k','facealpha',.5,'edgecolor','none')
                plot([RespROIs.V1az(animal,pos)-1 RespROIs.V1az(animal,pos)+1],[meanPeakAz_mean(animal,pos) meanPeakAz_mean(animal,pos)],'k')
                
                figure(frContraFig)
                median_N = median(N,2,'omitnan');
                frContra(animal,pos) = sum(median_N(4:end))./sum(median_N);
                scatter([RespROIs.V1az(animal,pos)],frContra(animal,pos),'k','filled')
            else
                disp(['mouse' num2str(animal) ', pos' num2str(pos) ' skipped'])
            end
        end
    end
    % - plot
    X = reshape(RespROIs.V1az,n_animals*max(n_pos),1);
    X = [ones(length(X),1) X];
    
    figure(scatterFig)
    Y = reshape(meanPeakAz_mean,n_animals*max(n_pos),1);
    [b,bint,~,~,stats]=regress(Y,X);
    yCalc = b(2).*[0 100]+b(1);
    yCalc1 = bint(2,1).*[0 100]+bint(1,1);
    yCalc2 = bint(2,2).*[0 100]+bint(1,2);
    plot([0 100],yCalc,'r')
    plot([-20 100],[1 13],'k:')
    fill([0 100 100 0],[yCalc1 flip(yCalc2)],'r','facealpha',.2,'edgecolor','none')
    ylim([0 13]);yticks([1 3 13]);yticklabels({'-20','0','100'})
    xlim([-20 110]); xticks([-20 0 100])
    title(['r2=' num2str(stats(1),4) ', p=' num2str(stats(3),4)])
    axis square
    
    figure(frContraFig)
    Y = reshape(frContra,n_animals*max(n_pos),1);
    [b,bint,~,~,stats]=regress(Y,X);
    yCalc = b(2).*[0 100]+b(1);
    yCalc1 = bint(2,1).*[0 100]+bint(1,1);
    yCalc2 = bint(2,2).*[0 100]+bint(1,2);
    plot([0 100],yCalc,'r')
    fill([0 100 100 0],[yCalc1 flip(yCalc2)],'r','facealpha',.2,'edgecolor','none')
    xlim([-20 110]); xticks([-20 0 100])
    ylim([0 1.2]);yticks([0 .5 1]);plot([0 110],[.5 .5],'k:')
    title(['r2=' num2str(stats(1)) ', p=' num2str(stats(3),4)])
    
end

%% -- Ipsi vs contra / Mono vs. Bino --
% % % - some definitions
% n_animals = length(animalID);
% ispos = ~cellfun(@isempty,RespROIs.info.data_details);
% n_pos = sum(ispos,2);
%
% azPos_toUse = false(16,1);
% azPos_toUse([1:5,7:2:end]) = true;
%
% base_window = RespROIs.info.base_window;
% resp_window = RespROIs.info.resp_window{1,1};
% timeVect = 0:1/6.0962:10;
% tBase = false(size(timeVect));
% tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
% tAna = false(size(timeVect));
% tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
%
% % % - do the calculation, per mice
% IpsiTuned = NaN(n_animals,max(n_pos));
% ContraTuned = NaN(n_animals,max(n_pos));
% idx_perMouse = cell(n_animals,max(n_pos));
% for mouse = 1:n_animals
%     for pos = 1:n_pos(mouse)
%         data = cat(1,RespROIs.data{mouse,pos});
%         nRep = size(data,4); nROIs = size(data,1);
%         Xvalidated = false(nROIs,1);
%         pearson_R = NaN(1,nROIs);
%         tAna_temp = tAna(1:size(data,2));
%         tBase_temp = tBase(1:size(data,2));
%         F0 = nanmean(data(:,tBase_temp,:,:),2);
%         dFF = (data-F0)./F0;
%         balanced_Resp = squeeze(nanmean(dFF(:,tAna_temp,:,:),2));
%         temp = randperm(nRep,nRep/2);
%         halfTrials = false(nRep,1); halfTrials(temp) = true;
%         A = balanced_Resp(:,:,halfTrials); B = balanced_Resp(:,:,~halfTrials);
%         for i = 1:nROIs
%             A2 = squeeze(A(i,:,:)); A2=A2(:);
%             B2 = squeeze(B(i,:,:)); B2=B2(:);
%             R_temp = corrcoef(A2,B2);
%             pearson_R(i) = R_temp(1,2);
%             if pearson_R(i)>0.2
%                 Xvalidated(i) = true;
%             end
%         end
%         balanced_meanResp = median(balanced_Resp(Xvalidated,:,:),3);
%         [~,temp]=max(balanced_meanResp,[],2);
%         idxMax = mod(temp,13); idxMax(idxMax==0)=13;
%         IpsiTuned(mouse,pos) = length(find(idxMax<4));
%         ContraTuned(mouse,pos) = length(find(idxMax>=4));
%         idx_perMouse{mouse,pos} = idxMax;
%     end
% end
% %
% frIpsi = sum(IpsiTuned,2,'omitnan')./(sum(IpsiTuned,2,'omitnan')+sum(ContraTuned,2,'omitnan'));
% k = sum(IpsiTuned,2,'omitnan')+sum(ContraTuned,2,'omitnan');
% nXval = find(k>10);
% frContra = 1-frIpsi;
%
% % % - plot
% figure; subplot(1,3,1);hold on
% N = NaN(n_animals,13);
% for mouse = 1:n_animals
%     if ismember(mouse,nXval)
%     IdxToPlot = cat(1,idx_perMouse{mouse,:});
%     N(mouse,:)= histcounts(IdxToPlot,[.5:1:13.5],'Normalization','probability');
%     plot(N(mouse,:))
%     end
% end
% plot(nanmean(N),'k','linewidth',3)
% fill([1:13,13:-1:1],[nanmean(N)-std(N,'omitnan')./sqrt(7) flip(nanmean(N)+std(N,'omitnan')./sqrt(7))],...
%     [0 0 0],'edgecolor','none','facealpha',.2)
% xlim([0 14]); xticks([1 3 8 13]); xticklabels({'-20','0','50','100'})
% ylim([0 .5]); yticks([0 .25 .5])
% set(gca,'tickdir','out')
% title('Peak azimuth, Xval')
%
% subplot(1,3,2); hold on
% scatter([ones(length(nXval),1); 2*ones(length(nXval),1)],[frIpsi(nXval);frContra(nXval)],'filled')
% plot(repmat([1;2],1,n_animals),[frIpsi,frContra]','color',[.5 .5 .5])
% xlim([.5 2.5]);xticks([1 2]); xticklabels({'Fraction ipsi','Fraction contra'})
% ylim([0 1]); yticks([0 .5 1])
% title(['mean/sem = ' num2str(nanmean(frContra(nXval)),3) '+-' num2str(std(frContra(nXval))./sqrt(length(nXval)),3)])
% set(gca,'tickdir','out')
%
% subplot(1,3,3); hold on
% V1az_all = RespROIs.V1az(:);
% % k=isnan(V1az_all);
% % V1az_all(k)=[];
% nXVal = IpsiTuned+ContraTuned;
% nXVal = nXVal(:);
% k = find(nXVal>10);
% frIpsi = IpsiTuned./(IpsiTuned+ContraTuned);
% frIpsi = frIpsi(:);
% V1az_all = RespROIs.V1az(:);
% scatter(V1az_all(k),frIpsi(k),'filled')
% ylim([0 1]); yticks([0 .5 1]);ylabel('Fraction ipsi')
% xlabel('V1 az ')
%
% X = [V1az_all(k) ones(length(k),1)];
% [b,~,~,~,stats]=regress(frIpsi(k),X);
% title(['regression,r2=' num2str(stats(1),3) ', p=' num2str(stats(3),3)])
%
% % % - save
% saveas(gcf,[SaveFolder filesep 'IpsiVsContra.tif']);
% saveas(gcf,[SaveFolder filesep 'IpsiVsContra.svg']);
% savefig([SaveFolder filesep 'IpsiVsContra']);
% pause(0.1);


% - as a function of V1 pos
% V1az_all = reshape(RespROIs.V1az',size(RespROIs.V1az,1)*size(RespROIs.V1az,2),1);
% V1az_all(isnan(V1az_all))=[];
% nROIs_all = cat(1,nROIs{:});
% nROIs_cumsum = cumsum(nROIs_all);
% nROIs_cumsum = [0;nROIs_cumsum];
% for animal = 1:length(nROIs_cumsum)-1
%     n_mono = sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),6:13)))+sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),19:26)))+sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),33:39)));
%     n_bino = sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),1:5))) +sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),14:18)))+sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),27:32)));
%     index1(animal) = (n_mono-n_bino)/(n_mono+n_bino);
%
%     n_ipsi = sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),1:3)))+sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),14:16)))+sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),27:29)));
%     n_contra = sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),4:13)))+sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),17:26)))+sum(sum(idex(nROIs_cumsum(animal)+1:nROIs_cumsum(animal+1),30:39)));
%     index2(animal) = (n_contra-n_ipsi)/(n_ipsi+n_contra);
% end
% figure; hold on
% yyaxis left
% scatter(V1az_all,index1,'filled','markerfacealpha',0.5)
% ylim([-1 1]); ylabel('Index #mono-#bino'); yticks([-1 -.5 .5 1])
% xlim([-30 110]); xlabel('V1 RF')
% plot([-30 110],[0 0],'k--')
% text(-30,+0.5,'mono','horizontalalignment','left','verticalalignment','top','Rotation',90)
% text(-30,-0.5,'bino','horizontalalignment','left','verticalalignment','top','Rotation',90)
% yyaxis right
% scatter(V1az_all,index2,'filled','markerfacealpha',0.5)
% ylim([-1 1]); ylabel('Index #contra-#ipsi'); yticks([-1 -.5 .5 1])
% text(110,+0.5,'contra','horizontalalignment','left','verticalalignment','bottom','Rotation',90)
% text(110,-0.5,'ipsi'  ,'horizontalalignment','left','verticalalignment','bottom','Rotation',90)
%
% x=V1az_all;
% y=index1';
% X = [ones(length(x),1) x];
% b = X\y;
% yCalc2 = X*b;
% plot(x,yCalc2,'b--')
% Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)
%
% y=index2';
% b = X\y;
% yCalc2 = X*b;
% plot(x,yCalc2,'r--')
% Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)

%% Axon bouton RF azimuth as function of that of V1
if do_RF
    % is working with new data structure and splits ROIs (1Sep2021)
    SaveFolder_RF = [SaveFolder filesep 'RF'];
    if ~exist(SaveFolder_RF,'dir')
        mkdir(SaveFolder_RF)
    end
    
    [RF,GOF,hasRF] = AxonRF_vs_V1RF(RespROIs,rsquare_th,Boutons,corrBoutons,SaveFolder_RF,[dataset '_' stim_toAnalyze]);
    
    % % --- RF distribution
    n_animals = length(animalID);
    RF_permice = cell(n_animals,1);
    for animal = 1:n_animals
        RF_permice{animal} = cat(2,RF{animal,:});
    end
    figure; hold on
    for animal = 1:n_animals % loop through the animals
        if length(RF_permice{animal})>10
            [N,~]  = histcounts(RF_permice{animal},-20:10:100,'Normalization','probability');
            plot(N,'color',[.5 .5 .5])
        end
    end
    
    RF_all = cat(2,RF_permice{:});
    [temp,~] = histcounts(RF_all,-20:10:100,'Normalization','probability');
    plot(temp,'k-','linewidth',3)
    xlabel('RF azimuth'); xticks(1:2:12);xlim([0 13]);xticklabels({'-20','0','20','40','60','80'})
    ylabel('Probability')
    nBoutons=length(cat(2,RF_permice{:}));
    if Boutons == 1
        ROItype='boutons';
    elseif Boutons == 0
        ROItype='axons';
    else
        ROItype='ROIs';
    end
    title(['Bouton RF center across mice' newline 'n=' num2str(nBoutons) ROItype])
    
    saveas(gcf,[SaveFolder_RF filesep ROItype '_RFcenterProba.tif'])
    savefig([SaveFolder_RF filesep ROItype '_RFcenterProba'])
end

%% SMI = f(RF)
if do_SMI && do_RF
    n_animals = size(RespROIs.info.animalID ,2);
    for animal= 1:n_animals
        RF_permice{animal} = cat(2,RF{animal,:});
        SMI_permice{animal} = cat(1,SMI{animal,:});
        hasRF_permice{animal} = cat(2,hasRF{animal,:});
    end
    
    clear GOF_all SMI_all
    RF_all = cat(2,RF_permice{:});
    hasRF_all = logical(cat(2,hasRF_permice{:}));
    SMI_all = cat(1,SMI_permice{:});
    
    figure;
    subplot(1,2,1)
    scatter(SMI_all(hasRF_all),RF_all,'filled','markerfacealpha',.2)
    xlabel('SMI'); ylabel('RF')
    subplot(1,2,2)
    scatter(RF_all,SMI_all(hasRF_all),'filled','markerfacealpha',.2)
    ylabel('SMI'); xlabel('RF')
end

%% Mean RF across boutons in the same V1 az position bin
if do_RF
    clear V1az_all is_member
    V1azbins = -20:20:120;
    % V1azbins(1) = -19;
    for animal = 1:n_animals
        V1az=[];
        for pos=1:length(RespROIs.V1az(animal,:))
            temp = RespROIs.V1az(animal,pos)*ones(length(RF{animal,pos}),1);
            V1az = [V1az;temp];
            V1az_all{animal} =  V1az;
        end
    end
    V1az_all = cat(1,V1az_all{:});
    
    [C,~,bin]=histcounts(V1az_all,V1azbins);
    
    for i = 1:length(V1azbins)
        tf = ismember(bin',i);
        is_member(i,:) = tf;
        RFavg_binned(i) = mean(RF_all(tf));
    end
    
    % -- plot figure
    figure;hold on
    scatter(V1az_all,RF_all,...
        'MarkerEdgeColor','none','MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.2);
    for i = 1:length(V1azbins)
        plot([V1azbins(i)+5 V1azbins(i)+15],[RFavg_binned(i) RFavg_binned(i)],'k-','linewidth',5)
    end
    ylim([-30 110]);xlim([-10 130]);
    plot([0 100],[0 100],'k--')
    axis equal
    xlabel('V1 neuron RF');ylabel(['AC ' ROItype ' RF'])
    
    saveas(gcf,[SaveFolder_RF filesep ROItype '_ACRFvsV1RF_binned.tif'])
    savefig([SaveFolder_RF filesep ROItype '_ACRFvsV1RF_binned'])
end

%% detect "notches" in the RF - don't do that if not needed, it take ages (with bootstrap)
if do_notches
    is_notche = 0;
    notches = cell(n_animals,7);
    f=waitbar(0);
    for animal = 1:n_animals
        waitbar(animal/n_animals,f)
        parfor pos=1:length(RespROIs.V1az(animal,:))
            if ~isempty(RF{animal,pos})
                %             waitbar(pos/length(RespROIs.V1az(animal,:)),f)
                %                 if RF{animal,pos}(i)<40
                %                     keyboard
                data = squeeze(RespROIs.data{animal,pos});
                dataID = 1:size(data,1);
                x_axis = [1:size(data,2)]/6-1;
                %                 nBoutons_max = max(RespROIs.nBoutonsPerROI{animal,pos});
                %                 for BoutonPerROI = 2:nBoutons_max
                %                     multipleBoutons = find(RespROIs.nBoutonsPerROI{animal,pos}==BoutonPerROI);
                %                     dataID = [dataID repmat(dataID(multipleBoutons),1,BoutonPerROI-1)];
                %                 end
                %                 temp=find(hasRF{animal,pos});
                %                 index=temp(i);
                %                 bouton_index = dataID(index);
                F0 = mean(data(:,1:6,:,:),2);
                F0 = repmat(F0,1,length(x_axis),1,1);
                nROIs = size(data,1);
                
                for ROI = 1:nROIs
                    Fo=squeeze(F0(ROI,:,:,:));
                    if any(Fo<=0,'all')
                        minF = min(data(ROI,:,:,:),[],'all');
                        data(ROI,:,:,:) = data(ROI,:,:,:)-minF;
                        F0(ROI,1,:,:) = nanmean(data(ROI,1:6,:,:),2);
                    end
                end
                dFF = (data-F0)./F0;
                
                %                     ymax=max(mean(dFF(bouton_index,:,:,:),4),[],'all');
                %                     ymin=min(mean(dFF(bouton_index,:,:,:),4),[],'all');
                %                     figure;
                %                 keyboard
                nootches= NaN(size(dFF,1),39);
                for bouton_index = 1:size(dFF,1)
                    %                 CI = NaN(39,2);
                    %                 for k =1:39 % all the speaker positions
                    %                     %                         subplot(3,13,k); hold on
                    %                     %                         toPlot = median(dFF(bouton_index,:,k,:),4);
                    %                     %                         ci = bootci(1000,@median,squeeze(dFF(bouton_index,:,k,:))');
                    %                     temp = bootci(1000,@median,squeeze(mean(dFF(bouton_index,7:15,k,:),2))');
                    %                     CI(k,:) = temp;
                    %                     %                         fill([x_axis flip(x_axis)]',[ci(1,:),flip(ci(2,:))],[0.7 0.7 0.7],'LineStyle','none')
                    %                     %                         plot(x_axis,toPlot,'k')
                    %                     %                         ylim([ymin-0.1 ymax+0.1])
                    %                 end
                    %                 medianResp=squeeze(mean(median(dFF(bouton_index,7:15,:,:),4),2));
                    %                 %                     figure; hold on
                    %                 c=0; temp=zeros(39,1);
                    %                 for kk=1:3
                    %                     %                         plot(medianResp(1+(kk-1)*13:13+(kk-1)*13)+10*(kk-1),'k')
                    %                     for k =1:13
                    %                         c=c+1;
                    %                         %                             plot([k k],[CI(k+(kk-1)*13,1) CI(k+(kk-1)*13,2)]+10*(kk-1),'k')
                    %                         if mod(c,13)~=0
                    %                             if medianResp(k+(kk-1)*13)>medianResp(k+1+(kk-1)*13)
                    %                                 temp(c)=CI(k+(kk-1)*13,1)-CI(k+1+(kk-1)*13,2);
                    %                             else
                    %                                 temp(c)=-CI(k+(kk-1)*13,2)+CI(k+1+(kk-1)*13,1);
                    %                             end
                    %                         end
                    %                     end
                    %                 end
                    %                 isnotch = find(temp>0);
                    %                 if ~isempty(isnotch)
                    %                     nootches(bouton_index,1:length(isnotch)) = isnotch;
                    %                 end
                    for k = 1:39
                        if mod(k,13)~=0 && mod(k,13)~=1
                            a=squeeze(mean(dFF(bouton_index,7:15,k-1,:),2));
                            c=squeeze(mean(dFF(bouton_index,7:15,k+1,:),2));
                            [~,h]=signrank(a,c);
                            if h
                                b=squeeze(mean(dFF(bouton_index,7:15,k,:),2));
                                [~,h1]=signrank(a,b);
                                [~,h2]=signrank(b,c);
                                if h1 && h2
                                    is_notche = is_notche+1;
                                end
                            end
                        end
                    end
                end % end loop through boutons
                notches{animal,pos} = nootches;
            end
        end % end loop through positions
    end % end loop through mice
    close(f)
end

%
% for animal = 1:n_animals
% notches_permice{animal} = cat(1,notches{animal,:});
% end
% notches_all = cat(1,notches_permice{:});
% n_notches = ~isnan(notches_all);
% n_notches = sum(n_notches,'all');
% fr_notches = n_notches/(size(notches_all,1)*36);

%% -- tuning curves --
if do_TuningCurves
    SaveFolder_TuningCurve = [SaveFolder filesep 'TuningCurves' filesep];
    if ~exist(SaveFolder_TuningCurve,'dir');mkdir(SaveFolder_TuningCurve);end
    %     TuningCurves(RespROIs,XVal,nShuffles,magn_th,smi_range,SMI,hasRF,SaveFolder_TuningCurve,[dataset '_' stim_toAnalyze])
    % until Feb 2023, following line:
    %     TuningCurves(RespROIs,XVal,nShuffles,magn_th,smi_range,SMI,[],SaveFolder_TuningCurve,[dataset '_' stim_toAnalyze])
    peakAz = TuningCurves_v2(RespROIs,XVal,nShuffles,magn_th,[],[],[],pearsonR_th,SaveFolder_TuningCurve,[dataset '_' stim_toAnalyze]);
    
end

%% Bayesian decoder
if do_decoder
    SaveFolder_Decoder = [SaveFolder filesep 'BayesianDecoder' filesep];
    if ~exist(SaveFolder_Decoder,'dir')
        mkdir(SaveFolder_Decoder)
    end
    Decoder_Output = BayesianDecoder_batch_v3(RespROIs,[dataset '_' stim_toAnalyze],SaveFolder_Decoder,doSingleSession,TwentykHz,plotSingleSession,doMouseagregate,doAccVsV1);
    
    % use v2 for the reviewer's figure or use the plot saved (LMaxons -
    % LED- 21 Aug)
    % Decoder_Output = BayesianDecoder_batch_v2(RespROIs,[dataset '_' stim_toAnalyze],SaveFolder_Decoder,1);
end

%% Use peak instead of fit
if do_Peak
    BoutonPeakResp_vs_V1RF;
end

%% grand averages
if do_GrandAverages
    n_animals = length(animalID);
    timeVect = 0:1/6.0962:7;
    base_window = RespROIs.info.base_window;
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    
    GrandAvg = NaN(n_animals,max(n_pos),39);
    GrandAvg_norm = NaN(n_animals,max(n_pos),39);
    for mouse = 1:n_animals
        for pos = 1:n_pos(mouse)
            data = RespROIs.data{mouse,pos};
            nROIs = size(data,1);
            
            resp_window = RespROIs.info.resp_window{mouse,pos};
            tAna = false(size(timeVect));
            tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
            tAna = tAna(1:size(data,2));
            tBase = false(size(timeVect));
            tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
            tBase = tBase(1:size(data,2));
            
            % - do a dFF
            F0 = nanmean(data(:,tBase,:,:),2);
            for ROI = 1:nROIs
                Fo=squeeze(F0(ROI,:,:,:));
                if any(Fo<0,'all')
                    minF = min(data(ROI,:,:,:),[],'all');
                    data(ROI,:,:,:) = data(ROI,:,:,:)-minF;
                    F0(ROI,1,:,:) = nanmean(data(ROI,tBase,:,:),2);
                end
            end
            dFF = (data-F0)./F0;
            
            % - transform ROIs into boutons
            nBoutons_max = max(RespROIs.nBoutonsPerROI{mouse,pos});
            for BoutonPerROI = 2:nBoutons_max
                multipleBoutons = RespROIs.nBoutonsPerROI{mouse,pos}==BoutonPerROI;
                dFF = cat(1,dFF,repmat(dFF(multipleBoutons,:,:,:),BoutonPerROI-1,1,1,1));
            end
            
            % - stimulus response
            temp = squeeze(nanmean(nanmean(dFF(:,tAna,:,:),2),4));
            GrandAvg(mouse,pos,:) = squeeze(nanmean(temp,1));
            GrandAvg_norm(mouse,pos,:) = squeeze(nanmean(temp./max(temp,2),1));
            
            % --- to plot indiv session data
%             SessionFig = figure;
%             traceAVgs = nanmean(dFF,4);
%             traceAVgs_norm = traceAVgs./max(nanmean(traceAVgs(:,tAna,:),2),[],3);
%             traceGrandAvg = squeeze(nanmean(traceAVgs_norm,1));
%             traceSEM = squeeze(std(traceAVgs_norm,[],1,'omitnan'))./sqrt(size(traceAVgs,1));
%             ymax = max(traceGrandAvg+traceSEM,[],'all'); ymin = min(traceGrandAvg-traceSEM,[],'all'); ymin = min([0 ymin]);
%             
%             for i = 1:39
%                 subtightplot(3,13,i); hold on
%                 toPlot = movmean(traceGrandAvg(:,i),3);
%                 SEM = movmean(traceSEM(:,i),3);
%                 plot([1 size(traceGrandAvg,1)],[0 0],'k:')
%                 fill([6 12 12 6 6],[-0.05 -0.05 0.5 0.5 -0.05],'k','edgecolor','g','facecolor','none')
%                 fill([24 30 30 24 24],[-0.05 -0.05 0.5 0.5 -0.05],'r','edgecolor','none','facealpha',0.2)
%                 fill([1:size(traceGrandAvg,1) size(traceGrandAvg,1):-1:1],[toPlot+SEM;flip(toPlot-SEM)]','k','edgecolor','none','facealpha',0.2)
%                 plot(toPlot,'k','linewidth',1)
%                 ylim([ymin ymax])
%                 xlim([1 36])
%                 if i ~= 27
%                     axis off
%                 end
%             end
            % ------------------------------
            
            clear data temp
        end % end loop thru positions
    end % end loop thru mice
    
    % % - plot the heatmap - optionala
    %     figure; c=1;
    %      for mouse = 1:n_animals
    %         for pos = 1:max(n_pos)
    %             if ispos(mouse,pos)
    %             subplot(n_animals,max(n_pos),c)
    %
    % %                 cmax = max(GrandAvg_norm(mouse,pos,:),[],'all');
    % %             reshaped = reshape(squeeze(GrandAvg_norm(mouse,pos,:)),3,13);
    %
    %             cmax = max(GrandAvg(mouse,pos,:),[],'all');
    %             reshaped = reshape(squeeze(GrandAvg(mouse,pos,:)),3,13);
    %             imagesc(reshaped,[0 cmax])
    %             end
    %             c=c+1;
    %         end
    %      end
    
    % % - smi
    SqrResp = GrandAvg.^2;
    GrandAvg_meanResp = nanmean(GrandAvg,3);
    respDiff = (GrandAvg-GrandAvg_meanResp).^2;
    smi = sum(respDiff,3)./sum(SqrResp,3);
    
    % - % figure
    figure;
    subplot(1,2,1);
    violinplot(smi(:));
    ylim([0 1]); yticks([0 .5 1])
    title('SMI Grand average - all sessions')
    subplot(1,2,2);
    scatter(ones(n_animals,1),nanmean(smi,2),'filled')
    ylim([0 1]); yticks([0 .5 1])
    set(gca,'TickDir','out')
    title('SMI Grand average - per mouse')
    
    SaveFolder_GrandAverage = [SaveFolder filesep 'GrandAverage' filesep];
    if ~exist(SaveFolder_GrandAverage,'dir')
        mkdir(SaveFolder_GrandAverage)
    end
    saveas(gcf,[SaveFolder_GrandAverage dataset '_avgAcrossSessions.tif']);
    saveas(gcf,[SaveFolder_GrandAverage dataset '_avgAcrossSessions.svg']);
    savefig([SaveFolder_GrandAverage dataset '_avgAcrossSessions']);
    
    % % - Best Az as a function of V1RF
    az_vector = -20:10:100;
    bestAz_pos = NaN(n_animals,max(n_pos));
    for mouse = 1:n_animals
        for pos = 1:n_pos(mouse)
            data = squeeze(GrandAvg(mouse,pos,:));
            temp = reshape(data,13,3);
            az_avgd = mean(temp,2);
            [~,bestAz_Idx] = max(az_avgd);
            bestAz_pos(mouse,pos) = az_vector(bestAz_Idx);
        end
    end
    V1az_all = RespROIs.V1az(:);
    k = isnan(V1az_all);
    V1az_all(k) = [];
    bestAz_pos_all = bestAz_pos(:);
    bestAz_pos_all(k) = [];
    
    figure;
    subplot(1,3,1); hold on
    scatter(V1az_all,bestAz_pos_all,'filled')
    xlim([-25 105]); xlabel('V1 RF')
    ylim([-25 105]); ylabel('GrdAvg best Az')
    axis square
    plot([-25 105],[-25 105],'k:')
        [b,bint,~,~,stats] = regress(bestAz_pos_all,[ones(size(V1az_all)),V1az_all]);
     yCalc = b(2).*[min(V1az_all) max(V1az_all)]+b(1);
    yCalc1 = bint(2,1).*[min(V1az_all) max(V1az_all)]+bint(1,1);
    yCalc2 = bint(2,2).*[min(V1az_all) max(V1az_all)]+bint(1,2);
    plot([0 100],yCalc,'r')
    fill([0 100 100 0],[yCalc1 flip(yCalc2)],'r','facealpha',.2,'edgecolor','none')
    title('95% CI, wald technique')
    
    subplot(1,3,2); hold on
    scatter(V1az_all,bestAz_pos_all,'filled')
    xlim([-25 105]); xlabel('V1 RF')
    ylim([-25 105]); ylabel('GrdAvg best Az')
    axis square
    plot([-25 105],[-25 105],'k:')
    ci = bootci(1000,{@regress,bestAz_pos_all,[ones(size(V1az_all)),V1az_all]}, ...
    'Type','normal');
    yCalc = b(2).*[min(V1az_all) max(V1az_all)]+b(1);
    yCalc_boot1 = ci(1,2).*[min(V1az_all) max(V1az_all)]+ci(1,1);
    yCalc_boot2 = ci(2,2).*[min(V1az_all) max(V1az_all)]+ci(2,1);
    plot([0 100],yCalc,'r')
    fill([0 100 100 0],[yCalc_boot1 flip(yCalc_boot2)],'r','facealpha',.2,'edgecolor','none')
    title('95% CI, bootstrap')

    subplot(1,3,3); hold on
    scatter(V1az_all,bestAz_pos_all,'filled')
    xlim([-25 105]); xlabel('V1 RF')
    ylim([-25 105]); ylabel('GrdAvg best Az')
    axis square
    plot([-25 105],[-25 105],'k:')
    yfit = [ones(size(V1az_all)),V1az_all]*b;
    resid = bestAz_pos_all - yfit;
    ci = bootci(1000,{@(bootr)regress(yfit+bootr,[ones(size(V1az_all)),V1az_all]),resid}, ...
    'Type','normal');
    yCalc = b(2).*[min(V1az_all) max(V1az_all)]+b(1);
    yCalc_boot1 = ci(1,2).*[min(V1az_all) max(V1az_all)]+ci(1,1);
    yCalc_boot2 = ci(2,2).*[min(V1az_all) max(V1az_all)]+ci(2,1);
    plot([0 100],yCalc,'r')
    fill([0 100 100 0],[yCalc_boot1 flip(yCalc_boot2)],'r','facealpha',.2,'edgecolor','none')
    title('95% CI, bootstrap residuals')

%     text(100,-20,['R2=' num2str(stats(1),4) ', p=' num2str(stats(3),4)],...
%         'horizontalalignment','right','verticalalignment','bottom')
        
    saveas(gcf,[SaveFolder_GrandAverage dataset 'BoutonVsV1RF.tif']);
    saveas(gcf,[SaveFolder_GrandAverage dataset 'BoutonVsV1RF.svg']);
    savefig([SaveFolder_GrandAverage dataset 'BoutonVsV1RF']);
    
    GrandAverageSourceData.bestAz_pos_all = bestAz_pos_all;
    GrandAverageSourceData.V1az_all = V1az_all;
    save([SaveFolder_GrandAverage 'GrandAverageSourceData'],'GrandAverageSourceData')
end

%% noise corraltions
if do_noiseCorrelations
    n_animals = length(animalID);
    timeVect = 0:1/6.0962:7;
    base_window = RespROIs.info.base_window;
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    
    rsq = NaN(n_animals,max(n_pos),13);
    rsqCtrl = NaN(n_animals,max(n_pos),13);
    %%
    % rsq_all = cell(n_animals,max(n_pos),100);
    for mouse = 1:n_animals
        for pos =1:n_pos(mouse)
            % - % do a dFF
            data = RespROIs.data{mouse,pos};
            nROIs = size(data,1);
            
            resp_window = RespROIs.info.resp_window{mouse,pos};
            tAna = false(size(timeVect));
            tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
            tAna = tAna(1:size(data,2));
            tBase = false(size(timeVect));
            tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
            tBase = tBase(1:size(data,2));
            
            F0 = nanmean(data(:,tBase,:,:),2);
            ROI_wNegBase = 0;
            for ROI = 1:nROIs
                Fo=squeeze(F0(ROI,:,:,:));
                if any(Fo<0,'all')
                    ROI_wNegBase= ROI_wNegBase+1;
                    minF = min(data(ROI,:,:,:),[],'all');
                    data(ROI,:,:,:) = data(ROI,:,:,:)-minF;
                    F0(ROI,1,:,:) = nanmean(data(ROI,tBase,:,:),2);
                end
            end
            dFF = (data-F0)./F0;
            
            % % - noise in the responses
            allResp = nanmean(dFF(:,tAna,:,:),2);
            meanResp = nanmean(nanmean(dFF(:,tAna,:,:),2),4);
            noise = squeeze(abs(allResp-meanResp));
            
            % % - decorrelate for a control
            noise_decorr = NaN(nROIs,39,size(allResp,4));
            for ROI = 1:nROIs
                for trialtype = 1:39
                    noise_decorr(ROI,trialtype,:) = noise(ROI,trialtype,randperm(size(allResp,4)));
                end
            end
            
            %% whole population
            noise_reshaped = reshape(noise,size(noise,1),39*size(allResp,4));
            noise_corr = corr(noise_reshaped','rows','pairwise');
            for k =1:length(noise_corr)
                noise_corr(k,k) = NaN;
            end
            rsq_perpos(mouse,pos) = mean(noise_corr,'all','omitnan');
            rsq_all{mouse,pos,:} = noise_corr(:);
            
            
            %% general correlations
            %         clear rsq
            [~,temp] = max(meanResp,[],3);
            %         Idx_maxResp = mod(temp,13);
            %         Idx_maxResp(Idx_maxResp==0)=13;
            Idx_maxResp=  temp;
            for trialtype = 1:39%13
                temp = noise(Idx_maxResp==trialtype,:,:);
                temp2 = noise(randperm(size(temp,1)),:,:);
                if ~isempty(temp)
                    noise_reshaped = reshape(temp,size(temp,1),39*size(allResp,4));
                    noise_corr = corr(noise_reshaped','rows','pairwise');
                    for k =1:length(noise_corr)
                        noise_corr(k,k) = NaN;
                    end
                    rsq(mouse,pos,trialtype) = mean(noise_corr,'all','omitnan');
                    
                    noiseCtrl_reshaped = reshape(temp2,size(temp,1),39*size(allResp,4));
                    noiseCtrl_corr = corr(noiseCtrl_reshaped','rows','pairwise');
                    for k =1:length(noise_corr)
                        noiseCtrl_corr(k,k) = NaN;
                    end
                    rsqCtrl(mouse,pos,trialtype) = mean(noiseCtrl_corr,'all','omitnan');
                end
            end
        end
    end
    %%
    rsq_perpos(rsq_perpos==0)=NaN;
    rsq_perpos_permice = nanmean(rsq_perpos,2);
    figure; hold
    scatter(3*ones(n_animals,1),rsq_perpos_permice,'r','filled')
    xlim([.5 3.5]);xticks([1:3]);xticklabels({'ACaud','LMaud','LMvis'})
    ylabel('r2')
    figure;hold on
    aaaah = cat(1,rsq_all{:});
    histogram(aaaah,[0:0.01:1],'facecolor',[1 1 1],'edgecolor','r')
    %%
    rSq_perMice = nanmean(nanmean(rsq,3),2);
    rSqCtrl_perMice = nanmean(nanmean(rsqCtrl,3),2);
    %
    figure; hold on
    plot(ones(2,n_animals).*[1;2],[rSq_perMice rSqCtrl_perMice]','color',[.5 .5 .5])
    scatter(ones(n_animals,1),rSq_perMice,'r','filled')
    scatter(2*ones(n_animals,1),rSqCtrl_perMice,'k','filled')
    [~,p]=ttest(rSq_perMice,rSqCtrl_perMice);
    title(['p=' num2str(p,3)])
    xlim([.5 2.5]);xticks([1 2]);xticklabels({'same az','random'})
    ylim([0 0.1]);yticks([0 .05 .1]);ylabel('r2')
    
    %% correlation as a function of stimulus azimtuh
    %         [~,temp] = max(meanResp,[],3);
    %         Idx_maxResp = mod(temp,13);
    %         Idx_maxResp(Idx_maxResp==0)=13;
    %         for trialtype = 1:13
    %             temp = noise(Idx_maxResp==trialtype,:,:);
    %             if ~isempty(temp)
    %                 for trialtype2 = 1:13
    %                 az2selec = trialtype2+13*[0;1;2];
    %                 noise_reshaped = reshape(temp(:,az2selec,:),size(temp,1),3*size(allResp,4));
    %                 noise_corr = corr(noise_reshaped','rows','pairwise');
    %                 for k =1:length(noise_corr)
    %                     noise_corr(k,k) = NaN;
    %                 end
    %                 rsq(mouse,pos,trialtype,trialtype2) = mean(noise_corr,'all','omitnan');
    %                 end
    %             end
    %         end
    %     end
    % end
    % temp = squeeze(nanmean(rsq,2));
    % % figure
    % % for i = 1:n_animals
    % %     for ii = 1:13
    % %         subplot(1,13,ii); hold on
    % %     dataToplot = squeeze(temp(i,ii,:));
    % %     plot(dataToplot)
    % %     end
    % %
    % % end
    % %    for ii = 1:13
    % %  subplot(1,13,ii); hold on
    % %  plot(squeeze(nanmean(temp(:,ii,:),1)),'k','linewidth',3)
    % %    end
    %
    % % sort according to distance to best az
    %    rsq_Sorted = NaN(n_animals,13,27);
    %    for j = 1:13
    %        rsq_Sorted(:,j,16-j:28-j) = temp(:,j,:);
    %        %           idx(j,:) = [1:13]-j;
    %    end
    %    temp2 = squeeze(nanmean(rsq_Sorted,2));
    %    delta_az = -120:10:120;
    %    %
    %    figure; hold on
    %    plot(delta_az,temp2./max(temp2,[],2))
    %    plot(delta_az,mean(temp2./max(temp2,[],2)),'k','linewidth',3)
    %    % average the same distance to best az
    %    for i = 1:12
    %    temp3(:,i) = mean([temp2(:,i),temp2(:,25-i+1)],2);
    %    end
    %    temp3(:,13) = temp2(:,13);
    %    % normalize
    %    figure; hold on
    %    plot([temp3./max(temp3,[],2)]')
    %    plot(mean(temp3./max(temp3,[],2)),'k','linewidth',3)
    %    anova1([temp3./max(temp3,[],2)])
    %% this is doing the same analysis as Boffi et al (Figure 3)
    % %         corrMatricesFig = figure;
    %         for j = 1:13
    %             noise_reshaped = reshape(noise(:,j,:),size(noise,1),20);
    %             %             noise_reshaped(noise_reshaped<0)=NaN;
    %             noise_corr = corr(noise_reshaped','rows','pairwise');
    %             %             T = clusterdata(noise_corr,'Linkage','complete','Cutoff',1);
    %                         %             k=T;
    %             cgObj = clustergram(noise_corr,'Colormap',redbluecmap);
    %             k = str2num(cell2mat(cgObj.RowLabels));
    %
    %             % % - shuffle
    %             noise_reshaped_decorr = reshape(noise_decorr(:,j,:),size(noise,1),20);
    %             noise_corr_decorr = corr(noise_reshaped_decorr','rows','pairwise');
    %             cgObj = clustergram(noise_corr_decorr,'Colormap',redbluecmap);
    %             k2 = str2num(cell2mat(cgObj.RowLabels));
    %
    %             for jj = 1:length(k)
    %                 maaatrix1(jj,:) = noise_corr(k(jj),k);
    %             end
    %             for jjj =1:13
    %                 noise_reshaped = reshape(noise(:,jjj,:),size(noise,1),20);
    %                 %                              noise_reshaped(noise_reshaped<0)=NaN;
    %                 noise_corr = corr(noise_reshaped','rows','pairwise');
    %                 for jj = 1:length(k)
    %                     maaatrix(jj,:) = noise_corr(k(jj),k);
    %                 end
    %                 temp = pdist2(maaatrix1,maaatrix,'cosine');
    % %                 for jjjj = 1:size(data,1)
    % %                     temp2 = temp(jjjj,jjjj);
    % %                 end
    % %                 D(j,jjj) = nanmean(temp2);
    %                 D(j,jjj) = mean(temp,'all','omitnan');
    %
    % %                 if j==1 && jjj == 1
    % %                     figure(corrMatricesFig)
    % %                     subplot(3,3,1); imagesc(maaatrix);
    % %                 elseif j==1 && jjj == 6
    % %                     figure(corrMatricesFig)
    % %                     subplot(3,3,2); imagesc(maaatrix);
    % %                 elseif j==1 && jjj == 13
    % %                     figure(corrMatricesFig)
    % %                     subplot(3,3,3); imagesc(maaatrix);
    % %                 end
    %             end
    %
    %              for jj = 1:length(k2)
    %                 maaatrix_decorr1(jj,:) = noise_corr_decorr(k2(jj),k2);
    %             end
    %             for jjj =1:13
    %                 noise_reshaped = reshape(noise_decorr(:,jjj,:),size(noise,1),20);
    %                 noise_corr = corr(noise_reshaped','rows','pairwise');
    %                 for jj = 1:length(k2)
    %                     maaatrix_decorr(jj,:) = noise_corr(k2(jj),k2);
    %                 end
    %                 temp = pdist2(maaatrix_decorr1,maaatrix_decorr,'cosine');
    % %                 for jjjj = 1:size(data,1)
    % %                     temp2 = temp(jjjj,jjjj);
    % %                 end
    % %                 D_decorr(j,jjj) = nanmean(temp2);
    %                 D_decorr(j,jjj) = mean(temp,'all','omitnan');
    %
    %             end
    %
    %         end
    %         %%
    %        figure; hold on
    %        histogram(D_decorr(:),[0.7:0.01:1],'facecolor',[1 1 1],'edgecolor','k')
    %        histogram(D(:),[0.7:0.01:1],'facecolor',[1 1 1],'edgecolor','r')
    %
    %        %%
    %        D_Sorted = NaN(13,27);
    %        for j = 1:13
    %            D_Sorted(j,16-j:28-j) = D(j,:);
    % %           idx(j,:) = [1:13]-j;
    %        end
    %
    %        %
    %        figure; hold on
    %        plot(nanmean(D_Sorted(:,15:end)))
    %     end
    % end
end
%%
diary off;
%% end of the function
% end
% function Batch_dataResponsive_v2(dataset,stim_toAnalyze)
%% Preambule
clearvars -except SMI_LMvis SMI_LMaud SMI_ACvis SMI_ACaud smi_LMvis smi_LMaud smi_ACaud Decoder_Output_AC Decoder_OutpuLMspkt Decoder_OutputLMled Decoder_Output_Bl6_80kHz Decoder_Output_CBA_80kHz Decoder_Output_RotSPK%RF_bl6_20kHz RF_bl6_80kHz RF_CBA RF_LMad RF_LMvis RF_bl6_vis
MainFolder = 'D:\AVspace_final';

dataset = 'ACaxons';           % 'ACaxons', 'LMaxons'
stim_toAnalyze = 'SPKwn';      % 'SPKwn', 'LED'
soundStim = '2-20';            % '2-20', '2-80' max fqcy (in kHz). Only if ACaxons
bckgd = 'bl6';                 % 'bl6', 'cba'. Only if 2-80 kHz sound


% ------------------------------------
% --------------- /!\ ----------------
removeFaceSession = false; % true means analysis on sessions w/o detectable face movement (see motion energy analysis)
% --------------- /!\ ----------------
% ------------------------------------


% -- Parameters used for Mazo et al., 2024 --
% - this is to reload the correct dataset ---
metric = 'dFF';                 % 'dFF'; 'spike'
selection_method = 'wilcoxon';  % 'ttest'; 'wilcoxon'; 'bootstrap'
alpha_threshold = 0.01;         % stat test alpha
amp_threshold = 0;              % amp above sd of baseline
dFF_threshold = 0.15;           %
% ------------------------------------------

%% select what to do
% - nROIs responsives and nBoutons
basicQuantif = false;

% - ROI sensitive to az/ele?
doStimPosSensitive = false; % -- Fig2c
do_elevation = false; % resample azimtuh to create a 3x3 grid, 20 deg steps

% - SMI
% for Fig 2d,f you need to save the 'SMI' variable and run 'SMI_plotALL.m' for plotting all data on the same graph
do_SMI = false;

% - PCA
do_PCA = false;  % -- not presented in Mazo et al., 2024

% - Cross validate Boutons
do_peakAzVsV1RF = true; % -- Fig 3g, Fig 4b
n_resample = 100;     % 100 in Mazo et al., 2024
pearsonR_th = 0.3;    % also used for tuning curves

% - ROI RF: fit a gaussian to the boutons' responses and compare with V1 RF (similar to Marques et al., Nat Neuro 2018)
% Alternative to Fig4b
% Not shown in Mazo et al., 2024
do_RF = false;
Boutons = 1;          % Boutons = 1, boutons; 2, axons; else ROIs
rsquare_th = 0.3;     % what is a good fit?

% - bayesian decoder
% for Fig 3c, you need to save the 'Decoder_Output' variable and run 'Decoder_v2.m' for plotting all data on the same graph
do_decoder = false; % -- Fig3, 4c, Supp Fig 4, Supp Fig 8
doSingleSession = true;    % Decoding error per session -- Fig 3b, Fig 4c, Supp Fig 4a,c. Some quirks with plotting
plotSingleSession = false; 
doMouseagregate = false;   % Mean decoding error as a function of number of axons -- Fig 3c
doAccVsV1 = false;         % Mean decoding error as a function of V1 RF -- Supp Fig 4b

% - Grand Average
do_GrandAverages = false;  % -- Fig2f and SuppFig7
plotIndivSessions = false;

% - tuning curves. uses 'pearsonR_th' and 'do_XvalBoutons'
do_TuningCurves = false; % -- Supp Fig 2e,f
nShuffles = 100;        % H0: there is no tuning; 100 used in Mazo et al., 2024
XVal = 1;               % whether to cross validate to determine best position or not

% - noise correlations
do_noiseCorrelations = false; % -- not presented in Mazo et al., 2024

%% Correct dependencies
if doAccVsV1
   doSingleSession = true; 
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
        
end

DataFolder = [MainFolder filesep '01_Data'  filesep '03_ResponsiveData' filesep];

dt = datestr(now,'yyyymmdd_HHMM');
if strcmp(stim_toAnalyze,'SPKwn')
    stimType = [stim_toAnalyze soundStim 'kHz' '_' bckgd];
else
    stimType = stim_toAnalyze;
end
temp = num2str(alpha_threshold); temp2 = strfind(temp,'.'); alpha = temp(temp2+1:end);
temp = num2str(dFF_threshold); temp2 = strfind(temp,'.'); dFF = temp(temp2+1:end);
SaveFolder = [MainFolder filesep '03_Plots' filesep dataset filesep stimType filesep ...
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
load([DataFolder filesep dataname '.mat'])
diaryName = 'AnalysisParameters';
diary([SaveFolder diaryName])
fprintf(1,'done \n')
disp(RespROIs.info);

RespROIs.info.stim_toAnalyze = stim_toAnalyze;

if removeFaceSession
    animalID = {'CMad103','CMad104','CMad97'};
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
    
    SaveFolder = [MainFolder filesep '03_Plots' filesep 'ExludeBehavior' filesep dataset filesep stimType filesep ...
        dt '_' selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF filesep];
    
    if ~exist(SaveFolder,'dir')
        mkdir(SaveFolder)
    end
end

%% Definition of a few useful variables
n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
timeVect = 0:1/6.0962:7;

%% quantif nBoutons and nROIs
if basicQuantif
    colors = distinguishable_colors(n_animals);
    
    % - fraction responsive
    figure; subplot(1,3,[1 2]);hold on
    for i = 1:n_animals
        frResp = RespROIs.fractionResp.nResp(i,~isnan(RespROIs.fractionResp.nResp(i,:)))./RespROIs.fractionResp.nROIs(i,~isnan(RespROIs.fractionResp.nROIs(i,:)));
        frResp(isnan(frResp))=[];
        scatter(i*ones(n_pos(i),1),frResp,...
            'filled','markerfacecolor',colors(i,:),'markerfacealpha',0.5)
        sh_mean = sum(RespROIs.fractionResp.nResp_sh(i,~isnan(RespROIs.fractionResp.nResp_sh(i,:))))./sum(RespROIs.fractionResp.nROIs(i,~isnan(RespROIs.fractionResp.nROIs(i,:))));
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
    SMI = SMI_v2(RespROIs,SaveFolderSMI); % save SMI to reuse it with 'SMI_olotALL'
end

%% cross-validate ROIs
% % - Peak az and ele of cross-validated ROIs
if do_peakAzVsV1RF
    n_animals = size(RespROIs.info.animalID ,2);
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    timeVect = 0:1/6.0962:7;
    
    % % - Cross-validate the boutons
    corrBoutons = cell(n_animals,max(n_pos),n_resample);
    idx_perMouse_az = cell(n_animals,max(n_pos),n_resample);
    idx_perMouse_ele = cell(n_animals,max(n_pos),n_resample);
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
                
                Resp = squeeze(nanmean(dFF(:,tAna,:,:),2));
                meanResp = mean(Resp(correlatedBoutons,:,:),3);

            temp = reshape(meanResp,size(meanResp,1),13,3);
            az_avged = mean(temp,3);
            ele_avged = squeeze(mean(temp,2));
            
            [~,idxMax_az]=max(az_avged,[],2);
            [~,idxMax_ele]=max(ele_avged,[],2);
            
            idx_perMouse_az{animal,pos,Rand} = idxMax_az;
            idx_perMouse_ele{animal,pos,Rand} = idxMax_ele;
            
                corrBoutons{animal,pos,Rand} = correlatedBoutons;
            end % end loop through random subsampling
        end % end loop through positions
    end % end loop through animals
    
% % - distribution - X-val boutons
N_az = NaN(13,n_resample);
N_ele = NaN(3,n_resample);
for Rand = 1:n_resample
    distrib_az = cat(1,idx_perMouse_az{:,:,Rand});
    N_az(:,Rand) = histcounts(distrib_az,[.5:1:13.5],'Normalization','probability');
    
    distrib_ele = cat(1,idx_perMouse_ele{:,:,Rand});
    N_ele(:,Rand) = histcounts(distrib_ele,[.5:1:3.5],'Normalization','probability');
    
    nROIs(Rand) = length(distrib_ele);
end
% - azimtuh
mean_N_az = median(N_az,2);
CI_N_az = prctile(N_az,[5 95],2);

% - elevation
mean_N_ele = median(N_ele,2);
CI_N_ele = prctile(N_ele,[5 95],2);

% - plot
az_vector = -20:10:100;
ele_vector = 20:-20:-20;
figure;
subplot(1,4,[1 3]);hold on
plot([1 13],[1/13 1/13],'k:')
fill([[1:13] [13:-1:1]],[CI_N_az(:,1);flip(CI_N_az(:,2))]','r','edgecolor','none','facealpha',0.5)
plot(mean_N_az,'k')
xlabel('Peak azimtuh');xlim([1,13]);xticks([1:13]);xticklabels(az_vector)
ylabel('Fraction');ylim([0 .25])

subplot(1,4,4);hold on
plot([1 3],[1/3 1/3],'k:')
fill([[1:3] [3:-1:1]],[CI_N_ele(:,1);flip(CI_N_ele(:,2))]','r','edgecolor','none','facealpha',0.5)
plot(mean_N_ele,'k')
xlabel('Peak elevation');xlim([1,3]);xticks([1:3]);xticklabels(ele_vector)
ylim([0 1])

    
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
                selected2 = (selected);%&smi_selec); % if you further want to select based on SMI...
                if sum(selected2)>10
                    AvgResp = squeeze(mean(mean(data(selected2,tAna,:,:),4),2));
                    AvgAzResp = reshape(AvgResp,size(AvgResp,1),13,3);
                    AvgAzResp = mean(AvgAzResp,3);
                    [~,peakResp] = max(AvgAzResp,[],2);
                    meanPeakAz(animal,pos,Rand) = median(peakResp);
                    
                    N(:,Rand) = histcounts(peakResp,[1:14],'Normalization','probability');
                    
                end
                
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
                
                %  --- not shown in Mazo et al., 2024 ----
                figure(frContraFig)
                median_N = median(N,2,'omitnan');
                frContra(animal,pos) = sum(median_N(4:end))./sum(median_N);
                scatter([RespROIs.V1az(animal,pos)],frContra(animal,pos),'k','filled')
                %  ---------------------------------------
                
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
    ylabel('Bouton best azimuth'); ylim([0 13]);yticks([1 3 13]);yticklabels({'-20','0','100'})
    xlabel('V1 RF az'); xlim([-20 110]); xticks([-20 0 100])
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
    xlabel('V1 RF az'); xlim([-20 110]); xticks([-20 0 100])
    ylabel('Fraction contra-preferring'); ylim([0 1.2]);yticks([0 .5 1]);plot([0 110],[.5 .5],'k:')
    title(['r2=' num2str(stats(1)) ', p=' num2str(stats(3),4)])
    
end % end X-validate

%% Axon bouton RF azimuth as function of that of V1
if do_RF
    SaveFolder_RF = [SaveFolder filesep 'RF'];
    if ~exist(SaveFolder_RF,'dir')
        mkdir(SaveFolder_RF)
    end
    
    [RF,GOF,hasRF] = AxonRF_vs_V1RF(RespROIs,rsquare_th,Boutons,SaveFolder_RF,[dataset '_' stim_toAnalyze]);
    
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


%% -- tuning curves --
if do_TuningCurves
    SaveFolder_TuningCurve = [SaveFolder filesep 'TuningCurves' filesep];
    if ~exist(SaveFolder_TuningCurve,'dir');mkdir(SaveFolder_TuningCurve);end
    TuningCurves_v2(RespROIs,XVal,nShuffles,pearsonR_th,SaveFolder_TuningCurve,[dataset '_' stim_toAnalyze]);
    
end

%% Bayesian decoder
if do_decoder
    SaveFolder_Decoder = [SaveFolder filesep 'BayesianDecoder' filesep];
    if ~exist(SaveFolder_Decoder,'dir')
        mkdir(SaveFolder_Decoder)
    end
    Decoder_Output = BayesianDecoder_batch_v3(RespROIs,[dataset '_' stim_toAnalyze],SaveFolder_Decoder,doSingleSession,TwentykHz,plotSingleSession,doMouseagregate,doAccVsV1);

end

%% grand averages
if do_GrandAverages
    SaveFolder_GrandAverage = [SaveFolder filesep 'GrandAverage' filesep];
    if ~exist(SaveFolder_GrandAverage,'dir')
        mkdir(SaveFolder_GrandAverage)
    end
    
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
            
            if plotIndivSessions
                SaveFolder_GrandAverageIndiv = [SaveFolder filesep 'GrandAverage' filesep 'IndivSessions'];
                mkdir(SaveFolder_GrandAverageIndiv)
                
                SessionFig = figure;
                traceAVgs = nanmean(dFF,4);
                traceAVgs_norm = traceAVgs./max(nanmean(traceAVgs(:,tAna,:),2),[],3);
                traceGrandAvg = squeeze(nanmean(traceAVgs_norm,1));
                traceSEM = squeeze(std(traceAVgs_norm,[],1,'omitnan'))./sqrt(size(traceAVgs,1));
                ymax = max(traceGrandAvg+traceSEM,[],'all'); ymin = min(traceGrandAvg-traceSEM,[],'all'); ymin = min([0 ymin]);
                
                for i = 1:39
                    subtightplot(3,13,i); hold on
                    toPlot = movmean(traceGrandAvg(:,i),3);
                    SEM = movmean(traceSEM(:,i),3);
                    plot([1 size(traceGrandAvg,1)],[0 0],'k:')
                    if strcmpi(stim_toAnalyze,'LED')
                        fill([24 30 30 24 24],[ymin ymin ymax ymax ymin]',...
                            'r','edgecolor','none','facealpha',.2)
                        fill([6 12 12 6 6],[ymin ymin ymax ymax ymin]',...
                            'k','facecolor','none','edgecolor','g','edgealpha',.2)
                    elseif strcmpi(stim_toAnalyze,'SPKwn')
                        fill([6 12 12 6 6],[ymin ymin ymax ymax ymin]',...
                            'g','edgecolor','none','facealpha',.2)
                        fill([24 30 30 24 24],[ymin ymin ymax ymax ymin]',...
                            'r','facecolor','none','edgecolor','r','edgealpha',.2)
                    end
                    fill([1:size(traceGrandAvg,1) size(traceGrandAvg,1):-1:1],[toPlot+SEM;flip(toPlot-SEM)]','k','edgecolor','none','facealpha',0.2)
                    plot(toPlot,'k','linewidth',1)
                    ylim([ymin ymax])
                    xlim([1 36])
                    if i ~= 27
                        axis off
                    end
                end
                ze_title = [RespROIs.info.animalID{animal} '_pos' num2str(pos)];
                sgtitle(ze_title,'interpreter','none')
                saveas(gcf,[SaveFolder_GrandAverageIndiv filesep dataset '_avgAcrossSessions.tif']);
                savefig([SaveFolder_GrandAverageIndiv filesep dataset '_avgAcrossSessions']);
            end
            
            clear data temp
        end % end loop thru positions
    end % end loop thru mice
    
    
    % % - SMI
    SqrResp = GrandAvg.^2;
    GrandAvg_meanResp = nanmean(GrandAvg,3);
    respDiff = (GrandAvg-GrandAvg_meanResp).^2;
    smi = sum(respDiff,3)./sum(SqrResp,3);
    
    % -- figure
    figure;
    subplot(1,2,1);
    violinplot(smi(:));
    ylabel('SMI'); ylim([0 1]); yticks([0 .5 1])
    title('SMI Grand average - all sessions')
    subplot(1,2,2);
    scatter(ones(n_animals,1),nanmean(smi,2),'filled')
    ylim([0 1]); yticks([0 .5 1])
    title('SMI Grand average - per mouse')
    
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
            [~,temp] = max(meanResp,[],3);
            
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
    
    rSq_perMice = nanmean(nanmean(rsq,3),2);
    rSqCtrl_perMice = nanmean(nanmean(rsqCtrl,3),2);
    
    % - plot
    figure; hold on
    plot(ones(2,n_animals).*[1;2],[rSq_perMice rSqCtrl_perMice]','color',[.5 .5 .5])
    scatter(ones(n_animals,1),rSq_perMice,'r','filled')
    scatter(2*ones(n_animals,1),rSqCtrl_perMice,'k','filled')
    [~,p]=ttest(rSq_perMice,rSqCtrl_perMice);
    title(['p=' num2str(p,3)])
    xlim([.5 2.5]);xticks([1 2]);xticklabels({'same az','random'})
    ylim([0 0.1]);yticks([0 .05 .1]);ylabel('r2 per mouse')
end

%%
diary off;
% end
%% Preambule
clearvars -except RespROIs
MainFolder = 'D:\AVspace_final';

[indx,tf] = listdlg('PromptString','choose the analysis','ListString',{'Ipsi vs Contra','Array vs Rotating'});
if indx == 1
    az_vector = -90:20:90;   % Analyze Ipsi vs Contra
else
    az_vector = -10:10:100;  % Analyze LED and speaker array vs Rotating speaker
end

% az_vector = -90:20:90;     % Analyze Ipsi vs Contra
% az_vector = -10:10:100;  % Analyze LED and speaker array vs Rotating speaker

load_the_data = true;

do_SMI     = false; % For Fig 9b
do_Peak    = false; % For Fig 9d

do_XvalBoutons = false; % for ipsi- vs contra-tuned - Fig 5c
rsq_th = 0.3;
n_subsampling = 100;
plotGrandAverage = false;

do_decoder = true;         % For Fig 5d
doSingleSession = true;   % 
plotIndivSession = false; 
doAnimalData = false;      % Agregate data per animal agregate -- Not presented in Mazo et al., 2024
PlotAnimalData = false;
do_MickeyMouse = false;     % Decoding error as a function of number of axons.
% For Supp Fig9c, save the 'Decoder_Output' and run 'Decoder_SingleRow.m' in parallel. Then run 'DecoderRotVsSingleRow.m' to combine the 2

%% load the data
if isequal(az_vector,-90:20:90)         % ipsi vs contra -- Fig 5
    animalID ={'CMad97','CMad98',...    %,'CMad103'};%,'CMad99'};          % Jan 2021
        'CMad122','CMad123','CMad124'};
elseif isequal(az_vector,-10:10:100)     % control for LED and Speaker array -- Supp Fig 9
    animalID ={'CMad97','CMad98'};
end
n_animals = length(animalID);

if load_the_data
    fprintf(1,'\nloading the data...')
    if isequal(az_vector,-90:20:90)
        loadName = [MainFolder filesep '01_Data\03_ResponsiveData\respData_RotSPK_dFF_wilcoxon_alpha01_magn15_amp0.mat']; % loads data from CMad97, 98, 122, 123 and 124
    elseif isequal(az_vector,-10:10:100)
        loadName = [MainFolder filesep '01_Data\03_ResponsiveData\respData_RotSPK_dFF_wilcoxon_alpha01_magn15_amp0_' num2str(length(animalID)) 'mice.mat'];  % loads data from CMad97, 98
    end
    load(loadName)
    fprintf(1,'done \n')
end

%% save folder
dt = datestr(now,'yyyymmdd_HHMM');
SaveFolder = [MainFolder filesep '03_Plots' filesep 'RotSPK' filesep dt filesep];
if length(az_vector)==10
    SaveFolder = [SaveFolder 'IpsiVsContra\'];
else
    SaveFolder = [SaveFolder 'SuppFig\'];
end
if ~exist(SaveFolder,'dir');mkdir(SaveFolder);end

%% -- SMI --
if do_SMI
    SaveFolderSMI = [SaveFolder filesep 'SMI' filesep];
    if ~exist(SaveFolderSMI,'dir');mkdir(SaveFolderSMI);end
    SMI = SMI_v2(RespROIs,SaveFolderSMI);
end

%% -- Ipsi vs contra --
if do_XvalBoutons
    % % - some definitions
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    
    azPos_toUse = false(15,1);
    azPos_toUse([1:6,8:2:end]) = true;

    for mouse = 1:n_animals
        if strcmp(animalID{mouse},'CMad122') || strcmp(animalID{mouse},'CMad123') || strcmp(animalID{mouse},'CMad124')
            timeVect = 0:1/6.0962:12; % frame rate, 6Hz; trial length, 12s.
            base_window = [0 1];
            resp_window = [1.2  2.2];
        else
            timeVect = 0:1/6.0962:10;
            base_window = [2 3];
            resp_window = [3.2  4.2];
        end
        
        data = cat(1,RespROIs.data{mouse,:});
        nRep = size(data,4); nROIs = size(data,1);
        Xvalidated = false(nROIs,1);
        
        tBase = false(size(timeVect));
        tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
        tBase = tBase(1:size(data,2));
        tAna = false(size(timeVect));
        tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
        tAna = tAna(1:size(data,2));
        
        % - do a dFF
        F0 = nanmean(data(:,tBase,:,:),2);
        dFF = (data-F0)./F0;
        
        balanced_Resp = squeeze(nanmean(dFF(:,tAna,azPos_toUse,:),2));
        
        for Rand = 1:n_subsampling
            Xvalidated = false(1,nROIs);
            temp = randperm(nRep,nRep/2);
            halfTrials = false(nRep,1); halfTrials(temp) = true;
            A = mean(balanced_Resp(:,:,halfTrials),3);
            B = mean(balanced_Resp(:,:,~halfTrials),3);
            for i = 1:nROIs
                A2=A(i,:);
                B2=B(i,:);
                R_temp = corrcoef(A2,B2);
                pearson_R(i) = R_temp(1,2);
                if pearson_R(i)>rsq_th
                    Xvalidated(i) = true;
                end
            end
            balanced_meanResp = median(balanced_Resp(Xvalidated,:,:),3);
            [~,idxMax]=max(balanced_meanResp,[],2);
            
            idx_perMouse{mouse,Rand} = idxMax;
        end
    end
    
    % % - cross-validated distribution of peak Azimuths
    idxPerSubsampling = cell(1,n_subsampling);
    for Rand = 1:n_subsampling
        idxPerSubsampling{Rand} = cat(1,idx_perMouse{:,Rand});
    end
    
    N = NaN(10,n_subsampling);
    nROIs = NaN(n_subsampling,1);
    for Rand = 1:n_subsampling
        distrib = cat(1,idx_perMouse{:,Rand});
        N(:,Rand) = histcounts(distrib,[1:1:11],'Normalization','probability');
        nROIs(Rand) = length(distrib);
    end
    
    mean_N = mean(N,2);
    CI_N = prctile(N,[5 95],2);
    
    figure;
    yl=[0 0.2];
    subplot(2,5,[1 2]);hold on
    plot([1 10],[.1 .1],'k:')
    fill([[1:10] [10:-1:1]],[CI_N(:,1);flip(CI_N(:,2))]','r','edgecolor','none','facealpha',0.5)
    plot(mean_N,'k')
    ylim(yl); ylabel('Fraction of boutons')
    xticks(1:10);xticklabels(az_vector); xlim([0 11]);xlabel('Peak azimuth')
    text(1,0,['n=' num2str(floor(mean(nROIs))) 'boutons'],...
        'horizontalalignment','left','verticalalignment','bottom')
    
    subplot(2,5,6:7); hold on
    rel_mean(1) = mean(mean_N([1,10]));
    rel_mean(2) = mean(mean_N([2,9]));
    rel_mean(3) = mean(mean_N([3,8]));
    rel_mean(4) = mean(mean_N([4,7]));
    rel_mean(5) = mean(mean_N([5,6]));
    
    N_rel(1,:) = squeeze(mean(N([1,10],:)));
    N_rel(2,:) = squeeze(mean(N([2,9],:)));
    N_rel(3,:) = squeeze(mean(N([3,8],:)));
    N_rel(4,:) = squeeze(mean(N([4,7],:)));
    N_rel(5,:) = squeeze(mean(N([5,6],:)));
    CI_N = prctile(N_rel,[5 95],2);
    
    plot([1 5],[.1 .1],'k:')
    fill([[1:5] [5:-1:1]],[CI_N(:,1);flip(CI_N(:,2))]','r','edgecolor','none','facealpha',0.5)
    plot(rel_mean,'k')
    ylim(yl); ylabel('Fraction of boutons')
    xticks(1:5);set(gca,'xdir','reverse');xticklabels(abs(az_vector(1:5))); xlim([0 6]);xlabel('Peak azimuth')
       
    % % - per mice
    N = NaN(n_animals,10,100);
    nROIs = NaN(n_animals,100);
    mean_N = NaN(n_animals,10);
    for i = 1:n_animals
        %     N = NaN(10,n_subsampling);
        for Rand = 1:n_subsampling
            distrib = cat(1,idx_perMouse{i,Rand});
            N(i,:,Rand) = histcounts(distrib,[1:1:11],'Normalization','probability');
            nROIs(i,Rand) = length(distrib);
        end
        mean_N(i,:) = mean(squeeze(N(i,:,:)),2);
    end
    
    subplot(2,5,3:4);hold on
    plot([1 10],[.1 .1],'k:')
    for i = 1:n_animals
        plot(mean_N(i,:),'color',[.5 .5 .5])
    end
    plot(mean(mean_N),'k','linewidth',2)
    xticks(1:10);xticklabels(az_vector); xlim([0 11]);xlabel('speaker pos.')
    t = table(mean_N(:,1),mean_N(:,2),mean_N(:,3),mean_N(:,4),mean_N(:,5),mean_N(:,6),mean_N(:,7),mean_N(:,8),mean_N(:,9),mean_N(:,10),...
        'VariableNames',{'meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10'});
    Meas = table(1:10,'VariableNames',{'Measurements'});
    rm = fitrm(t,'meas1-meas10~1');
    ranovatbl = ranova(rm);
    p = table2array(ranovatbl(1,5));
    titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
    c = multcompare(rm,'Time');
    title({'Per mice',titleText})
    
    
    rel_perMice(:,1) = mean(mean_N(:,[1,10]),2);
    rel_perMice(:,2) = mean(mean_N(:,[2,9]),2);
    rel_perMice(:,3) = mean(mean_N(:,[3,8]),2);
    rel_perMice(:,4) = mean(mean_N(:,[4,7]),2);
    rel_perMice(:,5) = mean(mean_N(:,[5,6]),2);
    subplot(2,5,8:9);hold on
    for i = 1:n_animals
        plot([mean(mean_N(i,[1,10]),2) mean(mean_N(i,[2,9]),2),...
            mean(mean_N(i,[3,8]),2) mean(mean_N(i,[4,7]),2) mean(mean_N(i,[5,6]),2)],...
            'color',[.5 .5 .5])
    end
    plot(mean(rel_perMice),'k','linewidth',2)
    t = table(rel_perMice(:,1),rel_perMice(:,2),rel_perMice(:,3),rel_perMice(:,4),rel_perMice(:,5),...
        'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
    rm = fitrm(t,'meas1-meas5~1');
    ranovatbl = ranova(rm);
    p = table2array(ranovatbl(1,5));
    xticks(1:5);set(gca,'xdir','reverse');xticklabels(abs(az_vector(1:5))); xlim([0 6]);xlabel('Peak azimuth')
    title(['ANOVA p=' num2str(p,4)])
    
    subplot(2,5,5);hold on
    FractionContra = squeeze(sum(N(:,6:10,:),2)./sum(N,2));
    CI_FrContra = prctile(FractionContra,[5 95],2);
    plot([0.5 5.5],[.5 .5],'k:')
    scatter([1:5]'.*ones(n_animals,1),mean(FractionContra,2),'k','filled')
    for i = 1:n_animals
        fill([.9 1.1 1.1 0.9]+(i-1),[CI_FrContra(i,1) CI_FrContra(i,1) CI_FrContra(i,2) CI_FrContra(i,2)],'k','facealpha',0.2)
    end
    xlim([.5 5.5]);xticks([1:n_animals]);xlabel('mouse #')
    ylim([0 1]); yticks([0 .5 1]); ylabel('fraction contra')
    [~,chi,p_chi] = crosstab(mean(FractionContra,2),[1-mean(FractionContra,2)]);
    title(['chi2, p =' num2str(p_chi,4)]);
    
    subplot(2,5,10);hold on
    FractionIpsi = squeeze(sum(N(:,1:3,:),2)./(sum(N,2).*3));
    scatter(ones(2,1),mean(FractionIpsi(1:2,:),2),'k','filled')
    scatter(ones(3,1),mean(FractionIpsi(3:5,:),2),'r','filled')
    FractionMid = squeeze(sum(N(:,4:7,:),2)./(sum(N,2).*4));
    scatter(2*ones(2,1),mean(FractionMid(1:2,:),2),'k','filled')
    scatter(2*ones(3,1),mean(FractionMid(3:5,:),2),'r','filled')
    FractionContra = squeeze(sum(N(:,8:10,:),2)./(sum(N,2).*3));
    scatter(3*ones(2,1),mean(FractionContra(1:2,:),2),'k','filled')
    scatter(3*ones(3,1),mean(FractionContra(3:5,:),2),'r','filled')
    
    for i = 1:n_animals
        plot([1:3],[mean(FractionIpsi(i,:),2) mean(FractionMid(i,:),2) mean(FractionContra(i,:),2)],'k')
    end
    
    plot([.9 1.1;1.9 2.1;2.9 3.1]',[mean(mean(FractionIpsi,2)) mean(mean(FractionMid,2)) mean(mean(FractionContra,2)); mean(mean(FractionIpsi,2)) mean(mean(FractionMid,2)) mean(mean(FractionContra,2))],'k','linewidth',3)
    t = table(mean(FractionIpsi,2),mean(FractionMid,2),mean(FractionContra,2),...
        'VariableNames',{'meas1','meas2','meas3'});
    rm = fitrm(t,'meas1-meas3~1');
    ranovatbl = ranova(rm);
    p = table2array(ranovatbl(1,5));

    c = multcompare(rm,'Time');
    mc = multcompare(rm,'Time', 'ComparisonType', 'tukey');
    plot([0.5 3.5],[.1 .1],'k:')
    xlim([.5 3.5]);xticks([1:3]);xticklabels({'Ipsi','Mid','Contra'})
    ylim([0 .3]); yticks([0 .1 .2 .3]); ylabel('fraction contra')
    title({['RM 1W ANOVA, p =' num2str(p,4)],['1vs2, p=' num2str(table2array(c(1,5)),4) '|2vs3,p=' num2str(table2array(c(4,5)),4)]});
    
    % % - save
    set(gcf,'units','normalized','position',[.1 .1 .8 .6])
    saveas(gcf,[SaveFolder filesep 'FractionIpsiVsContra.tif']);
    saveas(gcf,[SaveFolder filesep 'FractionIpsiVsContra.svg']);
    savefig([SaveFolder filesep 'FractionIpsiVsContra']);
    pause(0.1);
end

%% Bayesian decoder
if do_decoder
    SaveFolder_Decoder = [SaveFolder '\BayesianDecoder\'];
    
    if ~exist(SaveFolder_Decoder,'dir')
        mkdir(SaveFolder_Decoder)
    end
    Decoder_Output = BayesianDecoder_batch_RotSPK_v3(RespROIs,SaveFolder_Decoder,az_vector,doSingleSession,doAnimalData,do_MickeyMouse,PlotAnimalData,plotIndivSession);
end

%% PCA
if do_PCA
    SaveFolderPCA = [SaveFolder filesep 'PCA' filesep];
    if ~exist(SaveFolderPCA,'dir')
        mkdir(SaveFolderPCA)
    end
    
    performPCA_RotSPK(RespROIs,az_vector,SaveFolderPCA)
    
end
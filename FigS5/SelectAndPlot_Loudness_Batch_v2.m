%% CM Aug. 2022
% edits July2023
% v2: Nov2023. Uses responsive boutons in a loudness-dependent manner

% data is already dFF. F0 is trial baseline

% Call like RespROIs = SelectAndPlot_Loudness_Batch_v1(RespROIs) to reuse
% the data in the workspace
% SelectAndPlot_Loudness_Batch_v1 will either reselect the responsive ROIs
% or just reload the 'RespData' file, depending on "do_ROIselection"

function RespROIs = SelectAndPlot_Loudness_Batch_v2(varargin)
%% Paramters
% % - to select responsive ROIs or load appropriate data
do_ROIselection = false; % do it when new data added for instance
timeVect = 0:1/6.0962:7;
base_window = [0 1];
resp_window = [1.2  2.8];
metric = 'dFF';
selection_method = 'wilcoxon'; % 'ttest'; 'wilcoxon'; 'bootstrap'
alpha_threshold = 0.01;
magn_th = 0.15;

% % - do SMI, dFF and fraction of resp ROIs as a function of loudness
doFirstOrderMetrics = false;

% % - for the decoder
do_decoder = false;
doSingleSession = false;

% % - for bouton RF = f(V1RF)
do_boutonVsV1RF = true;
n_resample = 100;
pearsonR_th = 0.3;

MainFolder = 'C:\Users\camil\AudSpace\Data\ACaxonsLoudness';
subfix = ['_' metric '_' selection_method];
saveFolder = 'C:\Users\camil\AudSpace\Plots\ACaxonsLoudness\new';

% % - dataset
animalID ={'CMad73','CMad74',... % Apr 2021
    'CMad124'};                  % Aug 2023

V1az =[
    78.6 52.1 13.2 60.9 34   NaN;  % CMad73  - pos 1,5,8,9,12
    79.1 11.9 90.4 27.6 56.2 NaN;  % CMad74  - pos 2,3,7,8,11
    43.1 82.3 20.9 100 25.9 73.2]; % CMad124 - pos 06-11

%% get the analysis and base window as logicals
tAna=false(1,size(timeVect,2));
tAna(timeVect>resp_window(1)&timeVect<resp_window(2)) = true;
tBase=false(1,size(timeVect,2));
tBase(timeVect>base_window(1)&timeVect<base_window(2)) = true;

%% create ops with the ROI selection details - fed to 'SelectAndPlot' function
ops.timeVect = timeVect;
ops.tBase = tBase;
ops.tAna = tAna;
ops.metric = metric;
ops.selection_method = selection_method;
ops.alpha_th = alpha_threshold;
ops.magn_th = magn_th;
ops.nAttenuation = 4;
ops.loud_levels = 85:-10:55;
ops.nRep = 15;
ops.nAzPos = 13;
ops.az_vector = -20:10:100;
ops.n_animals = length(animalID);

%% load each dataset and select responsive ROIs
if do_ROIselection
    disp('reselect the ROIs')
    RespROIs.data = cell(ops.n_animals,size(V1az,2));
    for animal = 1:ops.n_animals
        fprintf(2,'\n %s \n', animalID{animal})
        files=subdir([MainFolder filesep animalID{animal}]);
        for pos = 1:length(files)
            if contains(files(pos).name,'dataSorted','IgnoreCase',true)
                
                str_num = regexp(files(pos).name,'pos');
                fprintf('%s',files(pos).name(str_num:str_num+4))
                load(files(pos).name);
                RespROIs.info.data_details{animal,pos}=Ftraces_merged.header.DataID;
                
                input_data = Ftraces_merged.data.tracesDist;
                
                [data,SMI,output_metrics] = SelectAndPlot_Loudness_v2(input_data,Ftraces_merged.nBoutonsPerROI,ops);
                RespROIs.data{animal,pos} = data;
                RespROIs.SMI{animal,pos} = SMI;
                RespROIs.output_metrics{animal,pos} = output_metrics;
                
                % pearson's corr hasn't been calculated, so need to do it here
                nROIs = size(Ftraces_merged.Fvalues,1);
                Fval = reshape(Ftraces_merged.Fvalues,nROIs,size(Ftraces_merged.Fvalues,2)*size(Ftraces_merged.Fvalues,3));
                R_fval = corrcoef(Fval');
                RespROIs.pearsonR{animal,pos} = R_fval;
                
            end
            fprintf('\n')
        end % end loop thru positions
    end % end loop thru mice
    
    RespROIs.info.animalID = animalID;
    RespROIs.info.selection_method = selection_method;
    RespROIs.info.alpha_threshold  = alpha_threshold;
    
    RespROIs.info.tBase = tBase;
    RespROIs.info.tAna  = tAna;
    
    RespROIs.V1az  = V1az;
    
    % % - save
    disp('saving data...')
    saveFolder = [MainFolder filesep 'RespData' filesep];
    if ~exist(saveFolder,'dir')
        mkdir(saveFolder)
    end
    dFF_th_Str = num2str(magn_th); k = strfind(dFF_th_Str,'.'); dFF_th_toSave = dFF_th_Str(k+1:end);
    temp = num2str(alpha_threshold); temp2 = strfind(temp,'.'); alpha = temp(temp2+1:end);
    saveName = [saveFolder metric '_' selection_method '_dot' alpha '_dot' dFF_th_toSave];
    save(saveName,'RespROIs', '-v7.3')
else
    if nargin==1
        disp('reuse RespData from the workspace')
        RespROIs = varargin{1};
    else
        dFF_th_Str = num2str(magn_th); k = strfind(dFF_th_Str,'.'); dFF_th_toSave = dFF_th_Str(k+1:end);
        temp = num2str(alpha_threshold); temp2 = strfind(temp,'.'); alpha = temp(temp2+1:end);
        dataname = [metric '_' selection_method '_dot' alpha '_dot' dFF_th_toSave];
        disp(['simply reload RespData file: ' dataname])
        load([MainFolder filesep 'RespData' filesep dataname '.mat'])
    end
end % do ROI selection

%% plots 'first-order' metrics
if doFirstOrderMetrics
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    % norm_mean = cell(ops.n_animals,1);
    colors = [1 0 0; 1 0 1;.5 0 1;0 0 1];
    SMIfig = figure;
    RepMagnFig = figure;
    for animal = 1:ops.n_animals
        % % - SMI distribution for different loudness
        figure(SMIfig)
        subplot(1,ops.n_animals+2,animal); hold on
        smi = cat(1,RespROIs.SMI{animal,:});
        for i = 1:ops.nAttenuation
            [f,x] = ecdf(smi(:,i));
            plot(f,x,'color',colors(i,:))
        end
        xlabel('SMI');
        if animal == 1
            ylabel('Cum. frequency')
        else
            ylabel('')
        end
        legend({'85dB','75dB','65dB','55dB'},'location','northwest')
        title(animalID{animal})
        
        subplot(1,ops.n_animals+2,ops.n_animals+1); hold on
        for i = 1:4
            [f,x] = ecdf(smi(:,i));
            p1 = plot(f,x,'color',colors(i,:));
            p1.Color = [colors(i,:) .5];
        end
        if animal == 1
            xlabel('SMI');
        end
        
        smi_meanPerMice(animal,:) = nanmean(smi);
        subplot(1,ops.n_animals+2,ops.n_animals+2); hold on
        scatter([1:4],nanmean(smi),'filled','markeredgecolor','none')
        
        % % - Resp magnitude per loudness
        figure(RepMagnFig)
        subplot(2,ops.n_animals+1,animal); hold on
        temp = cell(1,n_pos(animal));
        for i = 1:n_pos(animal)
            selec = RespROIs.output_metrics{animal,i}.selection;
            temp{i} = RespROIs.output_metrics{animal,i}.GlobalResp_acrossLoud(selec,:);
        end
        meandFF_acrossLoud = cat(1,temp{:});
        plot(nanmean(meandFF_acrossLoud,1),'o-')
        xlabel('dB');xticks(1:4); xticklabels({'85','75','65','55'});set(gca,'xdir','reverse')
        ax = gca; ylim([0 ax.YLim(2)]);
        if animal == 1
            ylabel('mean dF/F')
        end
        title(animalID{animal})
        
        % % - Fraction of resp ROIs per loudness -- not presented in Mazo et al., 2024
        subplot(2,ops.n_animals+1,ops.n_animals+1); hold on
        m = max(meandFF_acrossLoud,[],2,'omitnan');
        norm_mean = meandFF_acrossLoud./m;
        norm_mean(m<0,:) = [];
        plot(nanmean(norm_mean),'o-','color',[.5 .5 .5])
        norm_mean_perMouse(animal,:) = nanmean(norm_mean);
        
        subplot(2,ops.n_animals+1,ops.n_animals+1+animal); hold on
        temp = NaN(max(n_pos),ops.nAttenuation);
        for i = 1:n_pos(animal)
            temp(i,:) = sum(RespROIs.output_metrics{animal,i}.sel_acrossLoud)/size(RespROIs.output_metrics{animal,i}.sel_acrossLoud,1);
        end
        allResp_acrossLoud = nanmean(temp);
        plot(allResp_acrossLoud,'o-')
        xlabel('dB'); xticks(1:4); xticklabels({'85','75','65','55'});set(gca,'xdir','reverse')
        ax = gca; ylim([0 ax.YLim(2)]);
        if animal == 1
            ylabel('fr RespROIs')
        end
        RespAcrossLoud_perMice(animal,:) = allResp_acrossLoud;
        subplot(2,ops.n_animals+1,2*(ops.n_animals+1)); hold on
        plot(allResp_acrossLoud,'o-','color',[.5 .5 .5])
    end % end loop thru animals
    
    % % - SMI
    figure(SMIfig)
    subplot(1,ops.n_animals+2,ops.n_animals+1); hold on
    smi_all = cat(1,RespROIs.SMI{:});
    for i = 1:4
        [f,x] = ecdf(smi_all(:,i));
        plot(f,x,'color',colors(i,:),'linewidth',5)
    end
    Y = [smi_all(:,1),ones(size(smi_all,1),1)];
    Y= [Y;smi_all(:,2),2*ones(size(smi_all,1),1)];
    Y= [Y;smi_all(:,3),3*ones(size(smi_all,1),1)];
    Y= [Y;smi_all(:,4),4*ones(size(smi_all,1),1)];
    Y(isnan(Y),:)=[];
    ADtest = AnDarksamtest(Y);
    title(['AD test,p=' num2str(ADtest,4)])
    
    subplot(1,ops.n_animals+2,ops.n_animals+2); hold on
    plot([.9 1.9 2.9 3.9;1.1 2.1 3.1 4.1],[mean(smi_meanPerMice); mean(smi_meanPerMice)] ,'k')
    p = anova1(smi_meanPerMice,[],'off');
    title(['1W RM ANOVA, p=' num2str(p,4)])
    ylim([0 1]); ylabel('mean SMI')
    xlabel('Loudness (dB)');xticks(1:4);xticklabels(ops.loud_levels)
    set(gca,'xdir','reverse')
    
    % % - Resp magnitude
    figure(RepMagnFig)
    subplot(2,ops.n_animals+1,ops.n_animals+1); hold on
    plot(mean(norm_mean_perMouse),'k','linewidth',2)
    ylim([.2 .7]);ylabel('norm. dF/F')
    xlabel('dB');xticks(1:4); xticklabels({'85','75','65','55'});set(gca,'xdir','reverse')
    title('All mice')
    
    % % - Fraction responsive
    subplot(2,ops.n_animals+1,2*(ops.n_animals+1)); hold on
    plot(mean(RespAcrossLoud_perMice),'k-','linewidth',2)
    ylim([0 .2])
    xlabel('dB'); xticks(1:4); xticklabels({'85','75','65','55'});set(gca,'xdir','reverse')
    
    % % - save
    figure(SMIfig)
    set(gcf,'units','normalized','position',[.1 .1 .7 .2])
    savefig([saveFolder filesep 'SMI'])
    saveas(gcf,[saveFolder filesep 'SMI.tif']);
    saveas(gcf,[saveFolder filesep 'SMI.svg']);
    
    figure(RepMagnFig)
    savefig([saveFolder filesep 'RespVsLoudness'])
    saveas(gcf,[saveFolder filesep 'RespVsLoudness.tif']);
    saveas(gcf,[saveFolder filesep 'RespVsLoudness.svg']);
end

%% Decoder
if do_decoder
    Decoder_Output = BayesianDecoder_loudness_v1(RespROIs,saveFolder,ops,doSingleSession);
end

%% RF as a function of V1RF
if do_boutonVsV1RF
    
    % % - cross-validate first
    ispos = ~cellfun(@isempty,RespROIs.info.data_details);
    n_pos = sum(ispos,2);
    corrBoutons = cell(ops.n_animals,max(n_pos),n_resample,ops.nAttenuation);
    for att = 1:ops.nAttenuation
        for animal = 1:ops.n_animals
            for pos = 1:n_pos(animal)
                
                data = RespROIs.data{animal,pos};
                nROIs = size(data,1);
                
                selection = logical(RespROIs.output_metrics{animal,pos}.sel_acrossLoud(:,att));
                data = data(selection,:,:,:,:);
                % % - cross-correlation
                nRep = ops.nRep;
                for Rand = 1:n_resample
                    % % -- template
                    [template_data,idx] = datasample(data,floor(nRep/2),5,'Replace',false);
                    template_mean = squeeze(mean(mean(template_data(:,tAna,:,att,:),2),5));
                    
                    % % -- test
                    test_idx = true(nRep,1);
                    test_idx(idx) = false;
                    test_data = data(:,:,:,:,test_idx);
                    test_mean = squeeze(mean(mean(test_data(:,tAna,:,att,:),2),5));
                    
                    % % -- correlation
                    nROIs = size(template_mean,1);
                    R = NaN(1,nROIs);
                    for i = 1:nROIs
                        Rtemp = corrcoef(template_mean(i,:),test_mean(i,:));
                        R(i) = Rtemp(1,2);
                    end
                    correlatedBoutons = false(1,nROIs);
                    correlatedBoutons(R>pearsonR_th) = true;
                    
                    corrBoutons{animal,pos,Rand,att} = correlatedBoutons;
                end % end loop through random subsampling
            end % end loop through positions
        end % end loop through animals
    end % end loop thru attenuations
    
    scatterFig = figure; hold on
    meanPeakAz_mean = NaN(ops.n_animals,max(n_pos),ops.nAttenuation);
    meanPeakAz_CI = NaN(ops.n_animals,max(n_pos),2,ops.nAttenuation);
    meanPeakAz = NaN(ops.n_animals,max(n_pos),n_resample,ops.nAttenuation);
    for att = 1:ops.nAttenuation
        for animal = 1:ops.n_animals
            for pos = 1:n_pos(animal)
                for Rand = 1:n_resample
                    data = RespROIs.data{animal,pos};
                    selection = logical(RespROIs.output_metrics{animal,pos}.sel_acrossLoud(:,att));
                    data = data(selection,:,:,:,:);
                    selected = corrBoutons{animal,pos,Rand,att};
                    selected2 = selected; % that was to select based on SMI
                    
                    AvgAzResp = squeeze(mean(mean(data(selected2,tAna,:,att,:),2),5));
                    [~,peakResp] = max(AvgAzResp,[],2);
                    
                    % use median of the distribution...
                    meanPeakAz(animal,pos,Rand,att) = median(peakResp);
                    
                end % end loop thru randomizations
                
                if  sum(isnan(squeeze(meanPeakAz(animal,pos,:,att))))<n_resample
                    figure(scatterFig); subplot(1,ops.nAttenuation,att); hold on
                    meanPeakAz_mean(animal,pos,att) = median(meanPeakAz(animal,pos,:,att),3,'omitnan');
                    meanPeakAz_CI(animal,pos,:,att) = prctile(meanPeakAz(animal,pos,:,att),[5 95],3);
                    if RespROIs.V1az(animal,pos)>100
                        RespROIs.V1az(animal,pos) = 100;
                    end
                    fill([RespROIs.V1az(animal,pos)-1 RespROIs.V1az(animal,pos)+1 RespROIs.V1az(animal,pos)+1 RespROIs.V1az(animal,pos)-1],...
                        [squeeze(meanPeakAz_CI(animal,pos,1,att)) squeeze(meanPeakAz_CI(animal,pos,1,att)) squeeze(meanPeakAz_CI(animal,pos,2,att)) squeeze(meanPeakAz_CI(animal,pos,2,att))],...
                        'k','facealpha',.5,'edgecolor','none')
                    plot([RespROIs.V1az(animal,pos)-1 RespROIs.V1az(animal,pos)+1],[meanPeakAz_mean(animal,pos,att) meanPeakAz_mean(animal,pos,att)],'k')
                    
                else
                    disp(['mouse' num2str(animal) ', pos' num2str(pos) ' skipped'])
                end
            end % end loop thru pos
        end % end loop thru animals
    end % end loop thru attenuations
    
    % % - plot
    X = reshape(RespROIs.V1az,ops.n_animals*max(n_pos),1);
    X = [ones(length(X),1) X];
    
    for att = 1:ops.nAttenuation
        figure(scatterFig)
        subplot(1,ops.nAttenuation,att); hold on
        Y = reshape(meanPeakAz_mean(:,:,att),ops.n_animals*max(n_pos),1);
        [b,bint,~,~,stats]=regress(Y,X);
        
        % % - regress, wald technique
        %         yCalc = b(2).*[0 100]+b(1);
        %         yCalc1 = bint(2,1).*[0 100]+bint(1,1);
        %         yCalc2 = bint(2,2).*[0 100]+bint(1,2);
        
        % % - bootstrap residuals
        yfit = X*b;
        resid = Y - yfit;
        ci = bootci(1000,{@(bootr)regress(yfit+bootr,X),resid}, ...
            'Type','normal');
        yCalc = b(2).*[0 100]+b(1);
        yCalc1 = ci(1,2).*[0 100]+ci(1,1);
        yCalc2 = ci(2,2).*[0 100]+ci(2,1);
        
        plot([0 100],yCalc,'r')
        plot([-20 100],[1 13],'k:')
        fill([0 100 100 0],[yCalc1 flip(yCalc2)],'r','facealpha',.2,'edgecolor','none')
        ylim([0 13]);yticks([1 3 13]);yticklabels({'-20','0','100'})
        xlim([-20 110]); xticks([-20 0 100])
        title({['loudness: ' num2str(ops.loud_levels(att))],['r2=' num2str(stats(1),4) ', p=' num2str(stats(3),4)]})
        axis square
    end
    
    % % - save
    set(gcf,'units','normalized','position',[.1 .1 .8 .3])
    savefig([saveFolder filesep 'BountonVsV1RF'])
    saveas(gcf,[saveFolder filesep 'BountonVsV1RF.tif']);
    saveas(gcf,[saveFolder filesep 'BountonVsV1RF.svg']);
end
end
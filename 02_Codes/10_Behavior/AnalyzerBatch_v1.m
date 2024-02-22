%% November 2022
clear all;
SoundType = '2-20kHz'; %'2-20kHz','2-80kHz','all'

do_SPK = true;
do_LED = false; % not presented in Mazo et al., 2024
PlotIndivMouse = false;
PlotIndivSession = false;
doDecoder = true;

SaveFolder = ['D:\AVspace_final\03_Plots\10_BehaviorVideos\Pupil' filesep SoundType filesep '2DSpace'];
Dir = ['D:\AVspace_final\01_Data\10_Behavior' filesep];

%%
equation = 'y=%.2fx+%.2fb';
az_labels = -20:10:100;

if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end

%% get path and animal names
switch SoundType
    case '2-20kHz'
        directories{1} = [Dir '2-20kHz\ACaxons'];
        directories{2} = [Dir '2-20kHz\LMaxons'];
        ACanimals = {'CMad58','CMad62',...
            'CMad65','CMad67','CMad68'};
        LManimals = {'CMad85','CMad86','RD10278'};
        
    case '2-80kHz'
        directories{1} = [Dir '2-80kHz\Bl6'];
        %         directories{2} = [Dir '2-80kHz\CBA'];
        ACanimals = {'CMad97','CMad98'};
        %         CBAanimals = {'CMad103','CMad104','CMad106','CMad109','CMad110','CMad111'};
    case 'all'
end

%% load and store data in ACaxons and LMaxons structure
switch SoundType
    case '2-20kHz'
        for j = 1:length(directories)
            dirlist = dir(directories{j});
            mouse_counter = 0;
            for i = 1:length(dirlist)
                if any(strcmp(dirlist(i).name,ACanimals)) || any(strcmp(dirlist(i).name,LManimals)) % if a mouse folder
                    mouse_counter = mouse_counter+1;
                    pos_counter = 0;
                    poslist = dir([dirlist(i).folder filesep dirlist(i).name]);
                    for ii = 1:length(poslist)
                        if any(contains(poslist(ii).name,'eye_data_sorted'))
                            pos_counter = pos_counter+1;
                            temp = load([poslist(ii).folder filesep poslist(ii).name]);
                            if j == 1
                                ACaxons.eyedata{mouse_counter,pos_counter} = temp.eye_data_sorted;
                                ACaxons.info{mouse_counter,pos_counter} = [dirlist(i).name filesep poslist(ii).name];
                            else
                                LMaxons.eyedata{mouse_counter,pos_counter} = temp.eye_data_sorted;
                                LMaxons.info{mouse_counter,pos_counter} = [dirlist(i).name filesep poslist(ii).name];
                            end
                        end
                    end
                end
            end
        end
        eyedata = [ACaxons.eyedata;LMaxons.eyedata];
        animals = [ACanimals LManimals];
        
    case '2-80kHz'
        for j = 1:length(directories)
            dirlist = dir(directories{j});
            mouse_counter = 0;
            for i = 1:length(dirlist)
                if any(strcmp(dirlist(i).name,ACanimals))% || any(strcmp(dirlist(i).name,CBAanimals))% if a mouse folder
                    mouse_counter = mouse_counter+1;
                    pos_counter = 0;
                    poslist = dir([dirlist(i).folder filesep dirlist(i).name]);
                    for ii = 1:length(poslist)
                        if any(contains(poslist(ii).name,'eye_data_sorted'))
                            pos_counter = pos_counter+1;
                            temp = load([poslist(ii).folder filesep poslist(ii).name]);
                            if j == 1
                                Bl6mice.eyedata{mouse_counter,pos_counter} = temp.eye_data_sorted;
                                Bl6mice.info{mouse_counter,pos_counter} = [dirlist(i).name filesep poslist(ii).name];
                            else
                                CBAmice.eyedata{mouse_counter,pos_counter} = temp.eye_data_sorted;
                                CBAmice.info{mouse_counter,pos_counter} = [dirlist(i).name filesep poslist(ii).name];
                            end
                        end
                    end
                end
            end
        end
        eyedata = [Bl6mice.eyedata];%;CBAanimals.eyedata];
        animals = [ACanimals];% CBAanimals];
end

%% analyze stimulus responses
if do_SPK
    Resp = [1.5 3];
else
    Resp = [4.5 6];
end
if PlotIndivSession; indivDataFig = figure; end% this is the figure wiht the p-values per session
% pupilAreaVsAz = figure;
counter=0;
smi = NaN(length(animals),size(eyedata,2));
test_reshape_all = NaN(length(animals),size(eyedata,2),13,3);
avgTrace_perMouse = NaN(length(animals),39,120);

for mice = 1:length(animals)
    if PlotIndivMouse; MouseFig = figure; end
    counter = counter+1;
    clear session_avg
    n_pos =0;
    avgTrace = NaN(size(eyedata,2),39,120);
    for pos = 1:size(eyedata,2)
        if ~isempty(eyedata{mice,pos})
            n_pos = n_pos+1;
            if do_SPK
                area_SPK = eyedata{mice,pos}.SPK.area;
            else
                area_SPK = eyedata{mice,pos}.LED.area;
            end
            %             area_LED = eyedata{mice,pos}.LED.area;
            fps = size(area_SPK,3)/6;
            x_vect = 0:1/fps:6; x_vect=x_vect(1:size(area_SPK,3));
            
            tAna = false(size(area_SPK,3),1);
            tAna(x_vect>=Resp(1) & x_vect<Resp(2)) = true;
            tAna(end) = false; % this is for LED (coz one session has 121 frame, not 120), doesn't matter for spk
            session_avg(pos,:) = nanmean(nanmean(area_SPK(:,:,tAna),3),2);
            
            % -- reminder of quantif for each session ---------------------
            nRep = size(area_SPK,2);
            pupil_res = nanmean(area_SPK(:,:,tAna),3)';
            pupil_res_anova = reshape(pupil_res,nRep,13,3);
            pupil_res_anova2 = [pupil_res_anova(:,:,1);pupil_res_anova(:,:,2);pupil_res_anova(:,:,3)];
            % % have to use anovan because of NaNs
            pupil_res_anova3 = reshape(pupil_res_anova2,nRep*13*3,1);
            g1=[ones(nRep,1); 2*ones(nRep,1); 3*ones(nRep,1)];
            g1=repmat(g1,13,1);
            g2=[1:13].*ones(nRep*3,1);
            g2=reshape(g2,nRep*13*3,1);
            [p,~,~]=anovan(pupil_res_anova3,{g2 g1},'model','interaction','varnames',{'g1','g2'},'display','off');
            
            avgTrace(pos,:,1:size(area_SPK,3)) = squeeze(median(area_SPK,2,'omitnan'));
            if PlotIndivMouse
                figure(MouseFig);
                for i = 1:39
                    subplot(4,13,i); hold on
                    toPlot = squeeze(movmean(avgTrace(pos,i,:),10));
                    plot(toPlot,'color',[.5 .5 .5])
                    ylim([-1 2])
                end
%                 end % this was for plotting a specific mouse
            end
            
            if PlotIndivSession
                figure(indivDataFig)
                for i = 1:3
                    subplot(1,3,i);hold on;
                    scatter(counter,p(i),'r','filled');
                end
                k=find(p<0.05);
                for j = 1:length(k)
                    subplot(1,3,k(j));
                    text(counter,p(k(j)),num2str(pos),...
                        'horizontalalignment','center','verticalalignment','top')
                end
            end
            % -------------------------------------------------------------
                                    
        end
    end % end loop through imaging positions
    
    avgTrace_perMouse(mice,:,:) = squeeze(nanmean(avgTrace(:,:,1:120),1));
    if PlotIndivMouse
%         if mice == 4
%             keyboard
%         end
        figure(MouseFig);
        for i = 1:39
            subplot(4,13,i); hold on
            fill([20 40 40 20 20],[-1 -1 2 2 -1],'r','edgecolor','none','facealpha',0.2)
            toPlot = squeeze(movmean(avgTrace_perMouse(mice,i,:),10));
            plot(toPlot,'k')
%             if i ~= 27 
%                 axis off
%             end
        end
        
        X = [az_labels; ones(1,13)]';
        temp = reshape(mean(avgTrace_perMouse(mice,:,tAna),3),13,3);
        for i = 1:3
            subplot(4,13,[13*3+2*(i-1)+1 13*3+2*(i-1)+2]);hold on
            plot(temp(:,i),'k','linewidth',2)
            ylim([-0.5 1.5])
            
            y = temp(:,i);
            [b,~,~,~,stats] = regress(y,X);
            title([sprintf(equation , b(1)*100, b(2)) ' | p=' num2str(stats(3),3)])
            xlim([0 14]); xticks(1:13);xticklabels(az_labels);
            if i == 1
                ylabel('pupil area (s.d)'); xlabel('azimuth');
            end
        end
        
        subplot(4,13,13*3+7);hold on
        plot(nanmean(temp),'k','linewidth',2)
        y = nanmean(temp)';   X = [1:3;ones(1,3)]';
        [b,~,~,~,stats] = regress(y,X);
        title([sprintf(equation , b(1)*100, b(2)) ' | p=' num2str(stats(3),3)])
        xlim([0 4]); xticks(1:3);xticklabels({'+20','0','-20'});xlabel('elevation')
        
        pupil_res_anova = reshape(session_avg,n_pos,13,3);
        pupil_res_anova2 = [pupil_res_anova(:,:,1);pupil_res_anova(:,:,2);pupil_res_anova(:,:,3)];
        [p,tbl,stats]=anova2(pupil_res_anova2,n_pos,'off');
        
        %     pupil_res_anova3=reshape(pupil_res_anova2,n_pos*13*3,1);
        %     g1=[ones(n_pos,1); 2*ones(n_pos,1); 3*ones(n_pos,1)];
        %     g1=repmat(g1,13,1);
        %     g2=[1:13].*ones(n_pos*3,1);
        %     g2=reshape(g2,n_pos*13*3,1);
        %     [p,tbl,stats]=anovan(pupil_res_anova3,{g2 g1},'model','interaction','varnames',{'g1','g2'},'display','off');
        
        p_val(mice,:) = p;
        sgtitle([animals{mice} newline '2W ANOVA, p=' num2str(p,3)])
        
        set(gcf,'units','normalized','position',[.1 .3 .8 .6])
%         saveas(gcf,[SaveFolder filesep animals{mice} '.png'])
%         savefig([SaveFolder filesep animals{mice} '.fig'])
    end
end % end loop through mice

%% Decoder
if doDecoder
        BayesianDecoder_Pupil_batch(Resp,animals,eyedata,'allMice',SaveFolder)
end

%%  pupil area at the mouse level, across mice
if PlotIndivMouse
    n_color = 3;
    left_color = [1 0 0];
    right_color = [.5 0 .5];
    colors_r = [linspace(left_color(1),right_color(1),n_color)', linspace(left_color(2),right_color(2),n_color)', linspace(left_color(3),right_color(3),n_color)'];
    colors_r=[colors_r];
    
    time_labels = 0:5;
    
    AllMiceFig = figure;
    figure(AllMiceFig);
    for i = 1:39
        subplot(4,13,i);hold on
        fill([20,40,40,20,20],[-0.5,-0.5,1.5,1.5,-0.5 ],...
            [0.4660 0.6740 0.1880],'edgecolor','none','facealpha',0.2)
        for ii = 1:length(animals)
            p1 = plot(squeeze(avgTrace_perMouse(ii,i,:)),'color',[.5 .5 .5]);
        end
        p2 = plot(squeeze(nanmean(avgTrace_perMouse(:,i,:),1)));
        if i <14
            p2.Color = colors_r(1,:);
        elseif i < 27
            p2.Color = colors_r(2,:);
        else
            p2.Color = colors_r(3,:);
        end
        if i == 27
            xlabel('time (s)'); ylabel('pupil area (std)')
        end
        xticks([20:20:120]);xticklabels(time_labels)
        ylim([-0.5 1.5])
    end
    
    X = [az_labels; ones(1,13)]';
    
    subplot(4,13,40:41); hold on
    temp1 = nanmean(avgTrace_perMouse(:,1:13,tAna),3);
    toPlot = nanmean(temp1);
    sem = std(toPlot,'omitnan')./sqrt(5);
    fill([1:13,13:-1:1],[toPlot-sem flip(toPlot+sem)],...
        colors_r(1,:),'edgecolor','none','facealpha',0.5)
    plot(toPlot,'color',colors_r(1,:))
    ylim([0 1])
    y = nanmean(temp1)';
    [b,~,~,~,stats] = regress(y,X);
    p_regress(1) = stats(3);
    title([sprintf(equation , b(1)*100, b(2)) ' | p=' num2str(stats(3),3)])
    ylabel('pupil area (s.d)'); xlabel('azimuth');
    yCalc1 = X*b;
    plot(yCalc1,'k--');
    
    subplot(4,13,49:51); hold on
    fill([1:13,13:-1:1],[toPlot-sem flip(toPlot+sem)],...
        colors_r(1,:),'edgecolor','none','facealpha',0.5)
    plot(toPlot,'color',colors_r(1,:))
    
    
    subplot(4,13,43:44); hold on
    temp2 = nanmean(avgTrace_perMouse(:,14:26,tAna),3);
    toPlot = nanmean(temp2);
    sem = std(toPlot,'omitnan')./sqrt(5);
    fill([1:13,13:-1:1],[toPlot-sem flip(toPlot+sem)],...
        colors_r(2,:),'edgecolor','none','facealpha',0.5)
    plot(toPlot,'color',colors_r(2,:))
    ylim([0 1])
    y = nanmean(temp2)';
    [b,~,~,~,stats] = regress(y,X);
    p_regress(2) = stats(3);
    title([sprintf(equation , b(1)*100, b(2)) ' | p=' num2str(stats(3),3)])
    yCalc1 = X*b;
    plot(yCalc1,'k--');
    
    subplot(4,13,49:51); hold on
    fill([1:13,13:-1:1],[toPlot-sem flip(toPlot+sem)],...
        colors_r(2,:),'edgecolor','none','facealpha',0.5)
    plot(toPlot,'color',colors_r(2,:))
    
    
    subplot(4,13,46:47); hold on
    temp3 = nanmean(avgTrace_perMouse(:,27:39,tAna),3);
    toPlot = nanmean(temp3);
    sem = std(toPlot,'omitnan')./sqrt(5);
    fill([1:13,13:-1:1],[toPlot-sem flip(toPlot+sem)],...
        colors_r(2,:),'edgecolor','none','facealpha',0.5)
    plot(toPlot,'color',colors_r(2,:))
    ylim([0 1])
    y = nanmean(temp3)';
    [b,~,~,~,stats] = regress(y,X);
    p_regress(3) = stats(3);
    title([sprintf(equation , b(1)*100, b(2)) ' | p=' num2str(stats(3),3)])
    yCalc1 = X*b;
    plot(yCalc1,'k--');
    
    subplot(4,13,49:51); hold on
    fill([1:13,13:-1:1],[toPlot-sem flip(toPlot+sem)],...
        colors_r(3,:),'edgecolor','none','facealpha',0.5)
    plot(toPlot,'color',colors_r(3,:))
    ylim([0 1])
    xticks(1:13);xticklabels(az_labels);
    
    set(gcf,'units','normalized','position',[.1 .3 .8 .6])
    saveas(gcf,[SaveFolder filesep 'AllMice.svg'])
    saveas(gcf,[SaveFolder filesep 'AllMice.png'])
    savefig([SaveFolder filesep 'AllMice.fig'])
end

%% Plot across mice, more agregate quantification
temp = nanmean(avgTrace_perMouse(:,:,tAna),3);
az_tuning = reshape(temp,size(temp,1),13,3);
az_tuning = nanmean(az_tuning,3)./max(temp,[],2);
az_tuning_norm = az_tuning;%./max(az_tuning,[],2);
sd = std(az_tuning_norm)./sqrt(size(az_tuning_norm,1));

figure;
subplot(1,2,1); hold on
fill([1:13 13:-1:1],[mean(az_tuning_norm)+sd flip(mean(az_tuning_norm)-sd)],...
    [.5 .5 .5],'edgecolor','none','facealpha',.2)
plot(mean(az_tuning_norm), 'k')

ylim([0 1]); ylabel('norm dFF')
xticks(1:13);xticklabels(az_labels);

columnNames = {'az1'};
for i = 2:13
    columnNames = [columnNames, ['az' num2str(i)]];
end
t = array2table(az_tuning_norm,'VariableNames',columnNames);
Meas = table([1:13]','VariableNames',{'azimuth'});
rm = fitrm(t,'az1-az13 ~ 1','WithinDesign',Meas);
[rmTable] = ranova(rm);
p = table2array(rmTable(1,5));

%     [p,~,stats]=anova1(az_tuning_norm,[],'off');

X = [[-20:10:100]' ones(13,1)];
[b,~,~,~,stats] = regress(mean(az_tuning_norm)',X);
p_regress = stats(3);
yCalc1 = X*b;
plot(yCalc1,'k--');

title(['1w anova,p=' num2str(p,4) newline 'lin reg, p=' num2str(p_regress,4)])

subplot(1,2,2); hold on
ele_tuning = reshape(temp,size(temp,1),13,3);
ele_tuning = squeeze(nanmean(ele_tuning,2))./max(temp,[],2);
ele_tuning_norm = ele_tuning;%./max(az_tuning,[],2);
sd = std(ele_tuning_norm)./sqrt(size(ele_tuning_norm,1));

fill([1:3 3:-1:1],[mean(ele_tuning_norm)+sd flip(mean(ele_tuning_norm)-sd)],...
    [.5 .5 .5],'edgecolor','none','facealpha',.2)
plot(mean(ele_tuning_norm), 'k')
xticks(1:3);xticklabels([-20 0 20]);
ylim([0 1]); ylabel('norm dFF')

columnNames = {'ele1'};
for i = 2:3
    columnNames = [columnNames, ['ele' num2str(i)]];
end
t = array2table(ele_tuning_norm,'VariableNames',columnNames);
Meas = table([1:3]','VariableNames',{'elevation'});
rm = fitrm(t,'ele1-ele3 ~ 1','WithinDesign',Meas);
[rmTable] = ranova(rm);
p = table2array(rmTable(1,5));

Y = [[-20 0 20]' ones(3,1)];
[b,~,~,~,stats] = regress(mean(ele_tuning_norm)',Y);
p_regress2 = stats(3);
yCalc2 = Y*b;
plot(yCalc2,'k--');

text(20,1,['n=' num2str(length(animals)) 'mice'],...
    'horizontalalignment','left','verticalalignment','top')
title(['1w anova,p=' num2str(p,4) newline 'lin reg, p=' num2str(p_regress2,4)])

saveas(gcf,[SaveFolder filesep 'AllMiceQuantif.svg'])
saveas(gcf,[SaveFolder filesep 'AllMiceQuantif.png'])
savefig([SaveFolder filesep 'AllMiceQuantif.fig'])


%% compare LED and SPK responses -- not presented in Mazo et al., 2024
if do_LED
    Resp_SPK = [1.5 3];
    Resp_LED = [4.5 6];
    
    x_vect = 0:1/20:6;
    
    tAna_SPK = false(120,1);
    tAna_SPK(x_vect>=Resp_SPK(1) & x_vect<Resp_SPK(2)) = true;
    tAna_LED = false(120,1);
    tAna_LED(x_vect>=Resp_LED(1) & x_vect<Resp_LED(2)) = true;
    tAna_LED(end) = false;
    
    avgResp_SPK = NaN(length(animals),size(eyedata,2),39);
    avgResp_LED = NaN(length(animals),size(eyedata,2),39);
    for mice = 1:length(animals)
        n_pos =0;
        avgTrace = NaN(size(eyedata,2),39,120);
        for pos = 1:size(eyedata,2)
            if ~isempty(eyedata{mice,pos})
                n_pos = n_pos+1;
                area_SPK = eyedata{mice,pos}.SPK.area;
                area_LED = eyedata{mice,pos}.LED.area;
                avgResp_SPK(mice,pos,:) = squeeze(nanmean(nanmean(area_SPK(:,:,tAna_SPK),3),2));
                avgResp_LED(mice,pos,:) = squeeze(nanmean(nanmean(area_SPK(:,:,tAna_LED),3),2));
            end
        end % end loop through imaging positions
    end % end loop through mice
    
    avgResp_SPK_perMicePerPos = nanmean(nanmean(avgResp_SPK,3),2);
    avgResp_LED_perMicePerPos = nanmean(nanmean(avgResp_LED,3),2);
    %
    % jitter = randn(5,1)*0.2;
    figure; hold on
    plot([.5 2.5],[0 0],'k:')
    scatter(ones(length(animals),1),avgResp_SPK_perMicePerPos,[],[0.4660 0.6740 0.1880],'filled')
    scatter(2*ones(length(animals),1),avgResp_LED_perMicePerPos,[],[0 0.4470 0.7410],'filled')
    for i = 1:length(animals)
        plot([1 2],[avgResp_SPK_perMicePerPos(i) avgResp_LED_perMicePerPos(i)],'color',[.5 .5  .5])
    end
    plot([.8 1.2],[mean(avgResp_SPK_perMicePerPos) mean(avgResp_SPK_perMicePerPos)],'color',[0.4660 0.6740 0.1880],'linewidth',5)
    plot([1.8 2.2],[mean(avgResp_LED_perMicePerPos) mean(avgResp_LED_perMicePerPos)],'color',[0 0.4470 0.7410],'linewidth',5)
    xlim([0 3]);xticks([1,2]);xticklabels({'speaker','LED'})
    ylabel('pupil area (sd)')
    [~,p]=ttest(avgResp_SPK_perMicePerPos,avgResp_LED_perMicePerPos);
    title(['p ttest, p=' num2str(p,3)])
    ylim([-.5 1.5]) %jbtest or adtest (avgResp_SPK_perMicePerPos)
        
%     saveas(gcf,[SaveFolder filesep 'SPKvsLED.svg'])
%     saveas(gcf,[SaveFolder filesep 'SPKvsLED.png'])
%     savefig([SaveFolder filesep 'SPKvsLED.fig'])
end

%% Pupil-activity correlation
% -- load the data
% % - SMI
fprintf('load SMI...')
switch SoundType
    case '2-20kHz'
        load('C:\Users\camil\AudSpace\Data\SMI\Wilcoxon_dot01_0sd_0dot15ampTh\SMI_ACaud_LMaud.mat')
    case '2-80kHz'
        %         load('C:\Users\camil\AudSpace\Data\SMI\
end

% % - 2P data
for j = 1:length(directories)
    if j == 1
        switch SoundType
            case '2-20kHz'
                fprintf('load AC axons data...')
                eyedata = ACaxons.eyedata; animals = ACanimals;
                load('C:\Users\camil\AudSpace\Data\ACaxons\RespData\ACaxons_SPKwn_dFF_wilcoxon_dot01_0_dot15.mat')
                for i = 4:8
                    for ii = 1:7
                        if i == 5
                            if ii == 1
                            else
                                tempdata{i-3,ii-1} = RespROIs.data{i,ii};
                                temp_nBoutons{i-3,ii-1} = RespROIs.nBoutonsPerROI{i,ii};
                                temp_SMI{i-3,ii-1} = SMI_ACaud{i,ii};
                                temp_info{i-3,ii-1} = RespROIs.info.data_details{i,ii};
                                temp_respWindow{i-3,ii-1,:} = RespROIs.info.resp_window{i,ii};
                            end
                        elseif i == 8
                            if ii < 4
                                tempdata{i-3,ii} = RespROIs.data{i,ii};
                                temp_nBoutons{i-3,ii} = RespROIs.nBoutonsPerROI{i,ii};
                                temp_SMI{i-3,ii} = SMI_ACaud{i,ii};
                                temp_info{i-3,ii} = RespROIs.info.data_details{i,ii};
                                temp_respWindow{i-3,ii,:} = RespROIs.info.resp_window{i,ii};
                            elseif ii>4
                                tempdata{i-3,ii-1} = RespROIs.data{i,ii};
                                temp_nBoutons{i-3,ii-1} = RespROIs.nBoutonsPerROI{i,ii};
                                temp_SMI{i-3,ii-1} = SMI_ACaud{i,ii};
                                temp_info{i-3,ii-1} = RespROIs.info.data_details{i,ii};
                                temp_respWindow{i-3,ii-1,:} = RespROIs.info.resp_window{i,ii};
                            end
                        else
                            tempdata{i-3,ii} = RespROIs.data{i,ii};
                            temp_nBoutons{i-3,ii} = RespROIs.nBoutonsPerROI{i,ii};
                            temp_SMI{i-3,ii} = SMI_ACaud{i,ii};
                            temp_info{i-3,ii} = RespROIs.info.data_details{i,ii};
                            temp_respWindow{i-3,ii,:} = RespROIs.info.resp_window{i,ii};
                        end
                    end
                end % correct for missing pupil data
                RespROIs.data=tempdata;
                RespROIs.nBoutonsPerROI = temp_nBoutons;
                SMI_ACaud = temp_SMI;
                RespROIs.info.data_details = temp_info;
                RespROIs.info.resp_window = temp_respWindow;
                animalID = {'CMad58','CMad62','CMad65','CMad67','CMad68'};RespROIs.info.animalID = animalID;
                
            case '2-80kHz'
                fprintf('load Bl6 axons data...')
                eyedata = Bl6mice.eyedata; animals = ACanimals;
                load('C:\Users\camil\AudSpace\Data\ACaxons\RespData\ACaxons_SPKwn_2-80_bl6_dFF_wilcoxon_dot01_0_dot15.mat')
        end
    else
        switch SoundType
            case '2-20kHz'
                eyedata = LMaxons.eyedata; animals = LManimals;
                fprintf('load LM axons data...')
                load('C:\Users\camil\AudSpace\Data\LMaxons\RespData\LMaxons_SPKwn_dFF_wilcoxon_dot01_0_dot15.mat')
            case '2-80kHz'
        end
        fprintf('done!\n')
    end
    
    for mice = 1:length(animals)
        for pos = 1:size(eyedata,2)
            if ~isempty(RespROIs.data{mice,pos})
                %% ----------- draft
                data = RespROIs.data{mice,pos}  ;
                % % --- do a dFF
                timeVect = 0:1/6.0962:7;
                base_window = [0 1]; resp_window = [1.2  2.8];
                
                tAna_dFF = false(size(timeVect));
                tAna_dFF(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
                tAna_dFF = tAna_dFF(1:size(data,2));
                tBase_dFF = false(size(timeVect));
                tBase_dFF(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
                tBase_dFF = tBase_dFF(1:size(data,2));
                
                F0 = nanmean(data(:,tBase_dFF,:,:),2);
                ROI_wNegBase = 0;
                for ROI = 1:size(data,1)
                    Fo=squeeze(F0(ROI,:,:,:));
                    if any(Fo<0,'all')
                        ROI_wNegBase= ROI_wNegBase+1;
                        minF = min(data(ROI,:,:,:),[],'all');
                        data(ROI,:,:,:) = data(ROI,:,:,:)-minF;
                        F0(ROI,1,:,:) = nanmean(data(ROI,tBase_dFF,:,:),2);
                    end
                end
                dFF = (data-F0)./F0;
                
                % % -- pupil response
                fps = size(eyedata{mice,pos}.SPK.area,3)/6;
                x_vect = 0:1/fps:6; x_vect=x_vect(1:size(eyedata{mice,pos}.SPK.area,3));
            
                tAna = false(size(eyedata{mice,pos}.SPK.area,3),1);
                tAna(x_vect>=Resp(1) & x_vect<Resp(2)) = true;
                tAna(end) = false; % this is for LED (coz one session has 121 frame, not 120), doesn't matter for spk
                
                pupilresp = nanmean(eyedata{mice,pos}.SPK.area(:,:,tAna),3);
                pupilMeanresp = nanmean(pupilresp,2);
                pupil_noise = pupilresp-pupilMeanresp; % pupilresp
                %                 pupilresp_reshaped = reshape(pupilresp,39*20,1);
                pupilresp_reshaped = reshape(pupil_noise,39*20,1);
                k = isnan(pupilresp_reshaped);
                pupilresp_reshaped(k) = [];
                
                % % -- single ROI noise correlation
                dFFresp = squeeze(nanmean(dFF(:,tAna_dFF,:,:),2));
                dFFmeanresp = nanmean(dFFresp,3);
                dFF_noise = dFFresp-dFFmeanresp; % dFFresp
                
                rsq = NaN(1,size(dFF,1));
                rsq_shuffle = NaN(1,size(dFF,1));
                for ROI = 1:size(dFF,1)
                    indivROIdata = reshape(dFF_noise(ROI,:,:),39*20,1);
                    indivROIdata(k) = [];
                    temp = corrcoef(indivROIdata,pupilresp_reshaped);
                    rsq(ROI) = temp(1,2);
                                        
                    indivROIdata_shuffle = reshape(dFF_noise(ROI,:,:),39*20,1);
                    indivROIdata_shuffle = indivROIdata_shuffle(randperm(39*20),1);
                    indivROIdata_shuffle(k) = [];
                    temp = corrcoef(indivROIdata_shuffle,pupilresp_reshaped);
                    rsq_shuffle(ROI) = temp(1,2);
                end
                
                % % - ROI->boutons
                nBoutons_max = max(RespROIs.nBoutonsPerROI{mice,pos});
                for BoutonPerROI = 2:nBoutons_max
                    multipleBoutons = find(RespROIs.nBoutonsPerROI{mice,pos}==BoutonPerROI);
                    rsq = [rsq,repmat(rsq(multipleBoutons),1,BoutonPerROI-1)];
                    rsq_shuffle = [rsq_shuffle,repmat(rsq(multipleBoutons),1,BoutonPerROI-1)];
                end
                
                rsq_allROIs{j,mice,pos} = rsq;
                rsqSh_allROIs{j,mice,pos} = rsq_shuffle;
                
                
            end
        end
    end
    
end

% % - correlate all ROIs
max_nROIs = max(cellfun(@length,rsq_allROIs),[],'all');
LM_rsq = NaN(length(LManimals),size(rsq_allROIs,3),max_nROIs);
LM_rsqSh = NaN(length(LManimals),size(rsq_allROIs,3),max_nROIs);
for i = 1:length(LManimals)
    for ii = 1: size(rsq_allROIs,3)
        if ~isempty(rsq_allROIs{2,i,ii})
            LM_rsq(i,ii,1:size(rsq_allROIs{2,i,ii},2),1) = rsq_allROIs{2,i,ii};
            LM_rsqSh(i,ii,1:size(rsq_allROIs{2,i,ii},2),1) = rsqSh_allROIs{2,i,ii};
        end
    end
end
LM_rsq_perMice = nanmean(nanmean(LM_rsq,3),2);
LM_rsqSh_perMice = nanmean(nanmean(LM_rsqSh,3),2);

AC_rsq = NaN(length(ACanimals),size(rsq_allROIs,3));
AC_rsqSh = NaN(length(ACanimals),size(rsq_allROIs,3));
for i = 1:length(ACanimals)
    for ii = 1: size(rsq_allROIs,3)
        if ~isempty(rsq_allROIs{1,i,ii})
            AC_rsq(i,ii,1:size(rsq_allROIs{1,i,ii},2)) =  rsq_allROIs{1,i,ii};
            AC_rsqSh(i,ii,1:size(rsq_allROIs{1,i,ii},2)) =  rsqSh_allROIs{1,i,ii};
        end
    end
end
AC_rsq_perMice = nanmean(nanmean(AC_rsq,3),2);
AC_rsqSh_perMice = nanmean(nanmean(AC_rsqSh,3),2);

% plot average per mice
[~,pLM]=ttest(LM_rsq_perMice,LM_rsqSh_perMice);
[~,pAC]=ttest(AC_rsq_perMice,AC_rsqSh_perMice);

figure; hold on
scatter(ones(5,1)-0.2,AC_rsq_perMice,'b','filled')
scatter(ones(5,1)+0.2,AC_rsqSh_perMice,'k','filled')
plot(repmat([0.8;1.2],1,length(ACanimals)),[AC_rsq_perMice AC_rsqSh_perMice]','color',[.5 .5 .5])
scatter(2*ones(3,1)-0.2,LM_rsq_perMice,'r','filled')
scatter(2*ones(3,1)+0.2,LM_rsqSh_perMice,'k','filled')
plot(repmat([1.8;2.2],1,length(LManimals)),[LM_rsq_perMice LM_rsqSh_perMice]','color',[.5 .5 .5])
plot([0.5 0.7],[mean(AC_rsq_perMice) mean(AC_rsq_perMice)],'k','linewidth',3)
plot([.6 .6],[mean(AC_rsq_perMice)-std(AC_rsq_perMice)./sqrt(5) mean(AC_rsq_perMice)+std(AC_rsq_perMice)./sqrt(5)],...
    'k','linewidth',2)
plot([1.3 1.5],[mean(AC_rsqSh_perMice) mean(AC_rsqSh_perMice)],'k','linewidth',3)
plot([1.4 1.4],[mean(AC_rsqSh_perMice)-std(AC_rsqSh_perMice)./sqrt(5) mean(AC_rsqSh_perMice)+std(AC_rsqSh_perMice)./sqrt(5)],...
    'k','linewidth',2)
plot([1.55 1.75],[mean(LM_rsq_perMice) mean(LM_rsq_perMice)],'k','linewidth',3)
plot([1.65 1.65],[mean(LM_rsq_perMice)-std(LM_rsq_perMice)./sqrt(3) mean(LM_rsq_perMice)+std(LM_rsq_perMice)./sqrt(3)],...
    'k','linewidth',2)
plot([2.3 2.5],[mean(LM_rsqSh_perMice) mean(LM_rsqSh_perMice)],'k','linewidth',3)
plot([2.4 2.4],[mean(LM_rsqSh_perMice)-std(LM_rsqSh_perMice)./sqrt(3) mean(LM_rsqSh_perMice)+std(LM_rsqSh_perMice)./sqrt(3)],...
    'k','linewidth',2)
xlim([0 3]); xticks([1 2])
ylim([0 0.12]); yticks=([0 0.05 0.1]);
text(1,0.12,['p=' num2str(pAC,3)],...
    'horizontalalignment','center','verticalalignment','top')
text(2,0.12,['p=' num2str(pLM,3)],...
    'horizontalalignment','center','verticalalignment','top')

[~,p]=ttest2(AC_rsq_perMice,LM_rsq_perMice);
title(['indiv bouton vs pupil correlation' newline 'ttest, p=' num2str(p,4)])

savefig([SaveFolder filesep 'indivBoutons_AcrossMice'])
saveas(gcf,[SaveFolder filesep 'indivBoutons_AcrossMice.png'])
saveas(gcf,[SaveFolder filesep 'indivBoutons_AcrossMice.svg'])

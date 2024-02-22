%% 
% data is ROIs,frames,LED/SPKid,Brightness,SPKpos,loudness,nRep
% brightness (dim 4) = 1: SPK, no LEDs
% spk pos (dim 5) = 1 & loudness = 1 (dim 6): LED, no SPK
% Pos 7-8 are -20 and +20 in elevation -- not presented in Mazo et al., 2024

clearvars -except data_all

%% parametrize
MainFolder = 'D:\AVspace_final\01_Data\03_ResponsiveData\AV_Cells';
saveFolder = 'D:\AVspace_final\03_Plots\AVcells';

animalID = {'CMad76','CMad79','CMad80','CMad81',...
        'CAMP2','CAMP3','CMad133','CMad135'}; 

loudness_vec = [85 70 55];
brightness_vec = [1 3 5];
pos_labels = -40:20:40; % in azimuth, relative to V1RF (degrees)
pos_vec = [2 3 1 4 5];

selection_mthd = 'ttest'; % 'ttest','sd', 'snr', 'all'. ttest used in Mazo et al., 2024
alpha_th = 0.05;
magn_th = 0.05;

n_sh = 50;

tAna = 7:12;
tBase = 1:6;
nLEDs = 2;

doAllcells = false;

%% Load data
if ~exist('data_all')
    disp('loading data...');
    load([MainFolder filesep 'data_all_v5.mat']);
    
    clear file path
else
    disp('reusing data from workspace')
end

%% select responsive cells
n_animals = length(animalID);

visResp = cell(n_animals,1);
AVResp = cell(n_animals,1);
AVmod = cell(n_animals,1);
visTrace = cell(n_animals,1);
RespDist = cell(n_animals,1);
% audResp_traces = cell(n_animals,1);

disp(['selection: ' selection_mthd])

for mouse = 1:n_animals
    fprintf([animalID{mouse}])
    data = cat(1,data_all{mouse,:});
    nROIs = size(data,1);
    nRep = size(data,7);
    
    % % - initialize variables
    selection = false(1,nROIs);
    visResp{mouse} = NaN(nROIs,3);
    AVResp{mouse} = NaN(nROIs,3,7,3);
    AVmod{mouse} = NaN(nROIs,3,7,3);
    RespDist{mouse} = NaN(nROIs,3,8,4,7);    
    
    bestLED = squeeze(mean(mean(mean(mean(mean(data(:,tAna,:,2:4,2:6,2:4,:),2),4),5),6),7)); % use all LED trials (V and AV)
    [~,idx_bestLED] = max(bestLED,[],2);    
    
    testSD = NaN(1,nROIs);
    counter = 0; counter2 = 0;
    for i = 1:nROIs
        if strcmp(selection_mthd,'ttest') 
            % % - select using t-test...
            temp = squeeze(mean(data(i,tAna,idx_bestLED(i),2:4,2:6,2:4,:),2));
            avgResp_perTrial = mean(temp,4);
            maxResp = max(avgResp_perTrial,[],'all');
            VnAV_respDist = reshape(temp,3*5*3*nRep,1);
            temp = squeeze(mean(data(i,tBase,idx_bestLED(i),2:4,2:6,2:4,:),2));
            VnAV_baseDist = reshape(temp,3*5*3*nRep,1);
            h = ttest(VnAV_respDist,VnAV_baseDist,'alpha',alpha_th);
            
            if h && maxResp>magn_th %nanmean(VnAV_respDist)>magn_th
                temp = squeeze(data(i,tBase,idx_bestLED(i),2:4,2:6,2:4,:));
                temp = reshape(temp,numel(temp),1);
                sd = std(temp);
                if sd<2
                    selection(i) = true;
                end
            end
            
        elseif strcmp(selection_mthd,'sd')
            % % - ... or use sd instead
            temp = squeeze(mean(data(i,tAna,idx_bestLED(i),2:4,2:6,2:4,:),2));
            temp = reshape(temp,3*7*3*nRep,1);
            sd = std(temp);
            testSD(i) = sd; % for plotting the distribution of sd
            if sd<10%sd<5e-8 %sd<10
                selection(i) = true;
            end
            
        elseif strcmp(selection_mthd,'snr')
            % % - ... or use snr instead
            temp = squeeze(data(i,tBase,idx_bestLED(i),2:4,2:6,2:4,:));
            temp = reshape(temp,numel(temp),1);
            sd = std(temp);
            %             testSD(i) = sd;
            temp = squeeze(mean(data(i,tAna,idx_bestLED(i),2:4,2:6,2:4,:),2));
            m = max(mean(temp,4),[],'all');
            SNR = m/sd;
            testSD(i) = SNR;
            if SNR>1 
                selection(i) = true;
            end
            
        elseif strcmp(selection_mthd,'all')
            selection(i) = true;          
        end
        % same for both
        if selection(i)
            counter = counter+1;
            visResp{mouse}(counter,:) = squeeze(mean(mean(data(i,tAna,idx_bestLED(i),2:4,1,1,:),2),7));
            visTrace{mouse}(counter,:,:,:,:) = squeeze(mean(data(i,:,idx_bestLED(i),:,:,:,:),7));
            AVResp{mouse}(counter,:,:,:)  = squeeze(mean(mean(data(i,tAna,idx_bestLED(i),2:4,2:8,2:4,:),2),7));
            AVmod{mouse}(counter,:,:,:) = AVResp{mouse}(counter,:,:,:)-visResp{mouse}(counter,:);
            RespDist{mouse}(counter,:,:,:,:) = squeeze(mean(data(i,tAna,idx_bestLED(i),2:4,:,:,:),2));
        end
        
    end % end loop thru ROIs  
    fprintf([' | ' num2str(size(data,1)) ' cells, ' num2str(counter) ' resp ROIs\n'])
end % end loop thru mice  

%% traces with all cells -- Fig 6b,e
allCell = cat(1,visTrace{:});
nROIs = size(allCell,1);

figure; 
% % - V and AV traces, averaged across BR and LOUD
subplot(1,2,1); hold on
fill([5 11 11 5],[0 0 0.2 0.2],...
    'g','edgecolor','none','facealpha',.2)
vistrace = squeeze(nanmean(mean(allCell(:,:,2:4,1,1),3)));
avtrace1 = squeeze(nanmean(mean(mean(allCell(:,:,2:4,2,2:4),3),5)));
avtrace2 = squeeze(nanmean(mean(mean(allCell(:,:,2:4,3,2:4),3),5)));
avtrace3 = squeeze(nanmean(mean(mean(allCell(:,:,2:4,6,2:4),3),5)));
vistrace = movmean(vistrace,3);
avtrace1 = movmean(avtrace1,3);
avtrace2 = movmean(avtrace2,3);
avtrace3 = movmean(avtrace3,3);
sem = squeeze(std(mean(allCell(:,:,2:4,1,1),3)))./sqrt(nROIs);
fill([1:18 18:-1:1],[vistrace-sem flip(vistrace+sem)],...
    'r','edgecolor','none','facealpha',.2)
sem = squeeze(std(mean(mean(allCell(:,:,2:4,2,2:4),3),5)))./sqrt(nROIs);
fill([1:18 18:-1:1],[avtrace1-sem flip(avtrace1+sem)],...
    'b','facecolor',[0 .5 1],'edgecolor','none','facealpha',.2)
sem = squeeze(std(mean(mean(allCell(:,:,2:4,3,2:4),3),5)))./sqrt(nROIs);
fill([1:18 18:-1:1],[avtrace2-sem flip(avtrace2+sem)],...
    'b','facecolor',[0 .1 1],'edgecolor','none','facealpha',.2)
sem = squeeze(std(mean(mean(allCell(:,:,2:4,6,2:4),3),5)))./sqrt(nROIs);
fill([1:18 18:-1:1],[avtrace3-sem flip(avtrace3+sem)],...
    'b','edgecolor','none','facealpha',.2)

p1 = plot(vistrace,'r');
p2 = plot(avtrace1,'color',[0 .5 1]);
p3 = plot(avtrace2,'color',[0 .5 1],'linestyle','--');
p4 = plot(avtrace3,'b');
p6 = plot([1 18],[0 0],'k:');
legend([p1,p2,p3,p4],{'V','0','-40','+40'})
text(0,.2,['n=' num2str(nROIs) ' cells'],...
    'horizontalalignment','left','verticalalignment','top')

visResp = squeeze(mean(mean(allCell(:,tAna,2:4,1,1),2),3));
AVResp = squeeze(mean(mean(mean(mean(allCell(:,tAna,2:4,2:8,2:4),2),3),4),5));
[~,p]=ttest(visResp,AVResp);
title(['all cells | p. ttest,p=' num2str(p,4)])

% % - A only traces
subplot(1,2,2); hold on
fill([5 11 11 5],[0 0 0.2 0.2],...
    'g','edgecolor','none','facealpha',.2)
atrace1 = squeeze(nanmean(mean(allCell(:,:,1,2,2:4),5)));
atrace2 = squeeze(nanmean(mean(allCell(:,:,1,3,2:4),5)));
atrace3 = squeeze(nanmean(mean(allCell(:,:,1,6,2:4),5)));
atrace1 = movmean(atrace1,3);
atrace2 = movmean(atrace2,3);
atrace3 = movmean(atrace3,3);

sem = squeeze(std(mean(allCell(:,:,1,2,2:4),5)))./sqrt(nROIs);
fill([1:18 18:-1:1],[atrace1-sem flip(atrace1+sem)],...
    'b','facecolor',[0 .5 1],'edgecolor','none','facealpha',.2)
sem = squeeze(std(mean(allCell(:,:,1,3,2:4),5)))./sqrt(nROIs);
fill([1:18 18:-1:1],[atrace2-sem flip(atrace2+sem)],...
    'b','facecolor',[0 .1 1],'edgecolor','none','facealpha',.2)
sem = squeeze(std(mean(allCell(:,:,1,6,2:4),5)))./sqrt(nROIs);
fill([1:18 18:-1:1],[atrace3-sem flip(atrace3+sem)],...
    'b','edgecolor','none','facealpha',.2)

p2 = plot(atrace1,'color',[0 .5 1]);
p3 = plot(atrace2,'color',[0 .5 1],'linestyle','--');
p4 = plot(atrace3,'b');
p6 = plot([1 18],[0 0],'k:');
legend([p2,p3,p4],{'0','-40','+40'})
text(0,.2,['n=' num2str(nROIs) ' cells'],...
    'horizontalalignment','left','verticalalignment','top')

title('all cells | A resp')

set(gcf,'units','normalized','position',[.1 .1 .4 .3])
saveas(gcf,[saveFolder filesep 'Traces_allCells.svg'])
saveas(gcf,[saveFolder filesep 'Traces_allCells.jpg'])
savefig([saveFolder filesep 'Traces_allCells.fig'])

%% traces per mice - response to unimodal stim, w/ stats -- Fig6c, Supp Fig 10a,b
yl_mouse = [-.2,.4]; yl_all = [0 .3];
figure;
color_toUse = [0.5 .5 .5;.33 .33 .33; 0 0 0];
visResp = NaN(n_animals,3); audResp = NaN(n_animals,3);
for mouse = 1:n_animals
    allCell = cat(1,visTrace{mouse});
    nROIs = size(allCell,1);

    % % - V
    subplot(3,n_animals+1,mouse); hold on
    for BR = 2:4
        vistrace_mean = squeeze(nanmean(allCell(:,:,BR,1,1)));
        vistrace_mean = movmean(vistrace_mean,3);
        vistrace_sem = squeeze(std(allCell(:,:,BR,1,1)))./sqrt(nROIs);
        fill([1:18 18:-1:1],[vistrace_mean-vistrace_sem flip(vistrace_mean+vistrace_sem)],...
            color_toUse(BR-1,:),'edgecolor','none','facealpha',.2);
        plot(1:18,vistrace_mean,'color',color_toUse(BR-1,:));
        visResp(mouse,BR-1) = squeeze(mean(mean(allCell(:,tAna,BR,1,1),2),1));
    end
    if mouse == 1
        ylabel({'V','dF/F'})
    end
    title(['mouse' num2str(mouse) ' | n=' num2str(nROIs) ' neurons'])
    ylim(yl_mouse);
    
    % % - A
    subplot(3,n_animals+1,n_animals+1+mouse); hold on
    for LOUD = 2:4
        audtrace_mean = squeeze(nanmean(mean(allCell(:,:,1,2:6,LOUD),4)));
        audtrace_mean = movmean(audtrace_mean,3);
        audtrace_sem = squeeze(std(mean(allCell(:,:,1,2:6,LOUD),4)))./sqrt(nROIs);
        fill([1:18 18:-1:1],[audtrace_mean-audtrace_sem flip(audtrace_mean+audtrace_sem)],...
            color_toUse(3-LOUD+2,:),'edgecolor','none','facealpha',.2);
        plot(1:18,audtrace_mean,'color',color_toUse(3-LOUD+2,:));
        audResp(mouse,LOUD-1) = squeeze(mean(mean(mean(allCell(:,tAna,1,2:6,LOUD),2),4),1));
    end
    xlabel('Time (frames)')
    if mouse == 1
        ylabel({'A - avg across speakers','dF/F'})
    end
    ylim(yl_mouse)
    
    %     for SPKpos = 2:8
    subplot(3,n_animals+1,(n_animals+1)*2+mouse); hold on
    for LOUD = 2:4
        audRes_pos(mouse,LOUD-1,:) = squeeze(mean(mean(mean(allCell(:,tAna,1,2:6,LOUD),2),5),1));
        plot(squeeze(audRes_pos(mouse,LOUD-1,pos_vec)),'color',color_toUse(3-LOUD+2,:))
    end
    plot(squeeze(mean(audRes_pos(mouse,:,pos_vec),2)),'k','linewidth',2)
    xticks(1:5);xlim([.5 5.5]);xticklabels(pos_labels);xlabel('pos')
    ylim(yl_all)
    %     end
end

% % - V
subplot(3,n_animals+1,n_animals+1); hold on
plot(1:3,visResp,'color',[.5 .5 .5])
plot(1:3,mean(visResp),'k','linewidth',2)
xticks(1:3);xlim([.5 3.5]);xlabel('brightness')
ylim(yl_all)

t = table(visResp(:,1),visResp(:,2),visResp(:,3),...
    'VariableNames',{'meas1','meas2','meas3'});
Meas = table(1:3,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas3~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title({'Per mice',titleText})

% % - A, LOUD
subplot(3,n_animals+1,(n_animals+1)*2); hold on
plot(1:3,audResp,'color',[.5 .5 .5])
plot(1:3,mean(audResp),'k','linewidth',2)
set ( gca, 'xdir', 'reverse' );xticks([1:3]);xlim([.5 3.5]);xlabel('loudness')
ylim([0 .12])

t = table(audResp(:,1),audResp(:,2),audResp(:,3),...
    'VariableNames',{'meas1','meas2','meas3'});
Meas = table(1:3,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas3~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title(titleText)

% % - A, POS
subplot(3,n_animals+1,(n_animals+1)*3); hold on
audResp_pos_mean = squeeze(mean(audRes_pos(:,:,pos_vec),2));
plot(1:5,audResp_pos_mean,'color',[.5 .5 .5])
plot(1:5,mean(audResp_pos_mean,1),'k','linewidth',2)
xticks(1:5);xlim([.5 5.5]);xticklabels(pos_labels);xlabel('pos')
ylim([0 .12])

t = table(audResp_pos_mean(:,1),audResp_pos_mean(:,2),audResp_pos_mean(:,3),audResp_pos_mean(:,4),audResp_pos_mean(:,5),...
    'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
Meas = table(1:5,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas5~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title(titleText)

% % - finish plotting, (save)
sgtitle(['Pop resp to unimodal responses'])
set(gcf, 'units', 'normalized','position',[.05 .4 .9 .5])
saveas(gcf,[saveFolder filesep 'Traces_perMice.svg'])
saveas(gcf,[saveFolder filesep 'Traces_perMice.jpg'])
savefig([saveFolder filesep 'Traces_perMice.fig'])


%% unimodal responses to A @ 55dB -- Supp Fig 10g
figure; hold on
audResp_pos_mean = squeeze(mean(audRes_pos(:,3,pos_vec),2));
plot(1:5,audResp_pos_mean,'color',[.5 .5 .5])
plot(1:5,mean(audResp_pos_mean,1),'k','linewidth',2)
xticks(1:5);xlim([.5 5.5]);xticklabels(pos_labels);xlabel('pos')
ylim([0 0.12])

t = table(audResp_pos_mean(:,1),audResp_pos_mean(:,2),audResp_pos_mean(:,3),audResp_pos_mean(:,4),audResp_pos_mean(:,5),...
    'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
Meas = table(1:5,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas5~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title(titleText)
saveas(gcf,[saveFolder filesep 'A_55dB.svg'])
saveas(gcf,[saveFolder filesep 'A_55dB.jpg'])
savefig([saveFolder filesep 'A_55dB.fig'])


%% plot the effect per mouse, avged over brightness, loudness or pos -- Fig 6f and Supp Fig c,d
yl = [-.05 .05];

dataPos_perMice = NaN(n_animals,5);
dataLOUD_perMice = NaN(n_animals,3);
dataBR_perMice = NaN(n_animals,3);
figure;
for mouse = 1:n_animals
    dataPos_perMice(mouse,:) = squeeze(nanmean(mean(mean(AVmod{mouse}(:,:,1:5,:),2),4),1));
    dataLOUD_perMice(mouse,:) = squeeze(nanmean(mean(mean(AVmod{mouse}(:,:,1:5,:),2),3),1));
    subplot(1,4,3); hold on
    plot(dataPos_perMice(mouse,pos_vec),'color',[.5 .5 .5])
    subplot(1,4,2); hold on
    plot(dataLOUD_perMice(mouse,:),'color',[.5 .5 .5])
    dataBR_perMice(mouse,:) = squeeze(nanmean(mean(mean(AVmod{mouse}(:,:,1:5,:),3),4),1));
    subplot(1,4,1); hold on
    plot(dataBR_perMice(mouse,:),'color',[.5 .5 .5])
end

% % - BR
subplot(1,4,1);
plot([1 3],[0 0],'k:')
plot(mean(dataBR_perMice,1),'k','linewidth',2)
xlim([.5 3.5]);xticks([1:3]);xticklabels(brightness_vec);
ylim(yl); ylabel('AV-V (dFF)')

t = table(dataBR_perMice(:,1),dataBR_perMice(:,2),dataBR_perMice(:,3),...
    'VariableNames',{'meas1','meas2','meas3'});
Meas = table(1:3,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas3~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title({'brightness',titleText})

% % - LOUD
subplot(1,4,2);
plot([1 3],[0 0],'k:')
plot(mean(dataLOUD_perMice,1),'k','linewidth',2)
set(gca,'xdir','reverse')
xlim([.5 3.5]);xticks([1:3]);xticklabels(loudness_vec);
ylim(yl);

t = table(dataLOUD_perMice(:,1),dataLOUD_perMice(:,2),dataLOUD_perMice(:,3),...
    'VariableNames',{'meas1','meas2','meas3'});
Meas = table(1:3,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas3~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title({'loudness',titleText})

% % - Relative Az Position
subplot(1,4,3);
plot([1 7],[0 0],'k:')
plot(mean(dataPos_perMice(:,pos_vec),1),'k','linewidth',2)
xlim([.5 5.5]);xticks([1:5]);xticklabels(pos_labels);
ylim(yl); 
title('spk az pos')

t = table(dataPos_perMice(:,2),dataPos_perMice(:,3),dataPos_perMice(:,1),dataPos_perMice(:,4),dataPos_perMice(:,5),...
    'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
Meas = table(1:5,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas5~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title({'spk az pos',titleText})

% % - Absolute Az Position
subplot(1,4,4); hold on
clear temp
plot([1 3],[0 0],'k:')
temp(:,1) = mean(dataPos_perMice(:,[2,5]),2);
temp(:,2) = mean(dataPos_perMice(:,[3,4]),2);
temp(:,3) = dataPos_perMice(:,1);
plot(temp','color',[.5 .5 .5])
plot(mean(temp,1),'k','linewidth',2)
set(gca,'xdir','reverse')
xlim([.5 3.5]);xticks(1:3);xticklabels({'40','20','0'});
ylim(yl);

t = table(temp(:,3),temp(:,2),temp(:,1),...
    'VariableNames',{'meas1','meas2','meas3'});
Meas = table(1:3,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas3~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title({'absolute spk pos',titleText})

% % - Ele position
% subplot(1,5,5); hold on
% clear temp
% plot([1 4],[0 0],'k:')
% temp(:,1) = mean(dataPos_perMice(:,6),2);
% temp(:,2) = mean(dataPos_perMice(:,1),2);
% temp(:,3) = mean(dataPos_perMice(:,7),2);
% plot(temp','color',[.5 .5 .5])
% plot(mean(temp,1),'k','linewidth',2)
% set(gca,'xdir','reverse')
% xlim([.5 3.5]);xticks(1:4);xticklabels({'-20','0','20'});
% ylim(yl); yticks(-.1:.05:.1)
% title('spk ele pos')

% % - finish plotting, (save)
set(gcf, 'units', 'normalized','position',[.05 .4 .9 .3])
sgtitle('per mice')
saveas(gcf,[saveFolder filesep 'AVperMice_allQuantif.svg'])
saveas(gcf,[saveFolder filesep 'AVperMice_allQuantif.jpg'])
savefig([saveFolder filesep 'AVperMice_allQuantif.fig'])

%% plot all BRs and LOUDs -- not presented in Mazo et al., 2024
yl = [-.05 .1];
dataPos_perMice = NaN(n_animals,7,3,3);
figure;
for mouse = 1:n_animals
    for BR = 1:3
        for LOUD = 1:3
            dataPos_perMice(mouse,:,BR,LOUD) = squeeze(nanmean(AVmod{mouse}(:,BR,:,LOUD),1));
            subplot(3,3,LOUD+3*(BR-1)); hold on
            plot(dataPos_perMice(mouse,pos_vec,BR,LOUD),'color',[.5 .5 .5])
        end
    end
end
for BR = 1:3
    for LOUD = 1:3
        subplot(3,3,LOUD+3*(BR-1)); hold on
        % % - relative distance
        toPlot = dataPos_perMice(:,pos_vec,BR,LOUD);
        plot(mean(toPlot),'k','linewidth',2)
        ylim(yl)
        
        t = table(toPlot(:,1),toPlot(:,2),toPlot(:,3),toPlot(:,4),toPlot(:,5),...
            'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
        Meas = table(1:5,'VariableNames',{'Measurements'});
        rm = fitrm(t,'meas1-meas5~1');
        ranovatbl = ranova(rm);
        p = table2array(ranovatbl(1,5));
        titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
        
        % % - absolute distance
        %                 toPlot_abs(:,1) =  dataPos_perMice(:,1,BR,LOUD);
        %                 toPlot_abs(:,2) =  mean(dataPos_perMice(:,[2,5],BR,LOUD),2);
        %                 toPlot_abs(:,3) =  mean(dataPos_perMice(:,[3,4],BR,LOUD),2);
        %                  plot(mean(toPlot_abs),'k','linewidth',2)
        %                              ylim(yl)
        %              t = table(toPlot_abs(:,1),toPlot_abs(:,2),toPlot_abs(:,3),...
        %                 'VariableNames',{'meas1','meas2','meas3'});
        %             Meas = table(1:3,'VariableNames',{'Measurements'});
        %             rm = fitrm(t,'meas1-meas3~1');
        %             ranovatbl = ranova(rm);
        %             p = table2array(ranovatbl(1,5));
        %             titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
        
        title({['BR=' num2str(BR) ', LOUD=' num2str(LOUD)],titleText})
    end
end

%% AV @ 55dB, averaged across BR -- Supp Fig 10i
dataPos_perMice2 = NaN(n_animals,5);
figure;
subplot(1,2,1); hold on
for mouse = 1:n_animals
    dataPos_perMice2(mouse,:) = squeeze(nanmean(mean(AVmod{mouse}(:,:,pos_vec,3),2),1));
    plot(dataPos_perMice2(mouse,:),'color',[.5 .5 .5])
end

% % - relative distance
toPlot = dataPos_perMice2;
plot([1 5],[0 0],'k:')
plot(mean(toPlot),'k','linewidth',2)
ylim([-.05 .05])
xlabel('V1RF-speaker distance');xticks(1:5); xticklabels({'-40','-20','0','20','40'}); xlim([.5 5.5]);

t = table(toPlot(:,1),toPlot(:,2),toPlot(:,3),toPlot(:,4),toPlot(:,5),...
    'VariableNames',{'meas1','meas2','meas3','meas4','meas5'});
Meas = table(1:5,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas5~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));
titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title(titleText)

subplot(1,2,2) ; hold on
toPlot_abs(:,1) = dataPos_perMice2(:,3);
toPlot_abs(:,2) = mean(dataPos_perMice2(:,[2 4]),2);
toPlot_abs(:,3) = mean(dataPos_perMice2(:,[1 5]),2);
plot(toPlot_abs','color',[.5 .5 .5])
plot(mean(toPlot_abs),'k','linewidth',2)
ylim([-.05 .05])
xlabel('rel. V1RF-speaker distance');xticks(1:3); xticklabels({'0','20','40'}); xlim([.5 3.5]);

t = table(toPlot_abs(:,1),toPlot_abs(:,2),toPlot_abs(:,3),...
    'VariableNames',{'meas1','meas2','meas3'});
Meas = table(1:5,'VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas3~1');
ranovatbl = ranova(rm);
p = table2array(ranovatbl(1,5));

titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
title(titleText)

set(gcf,'units','normalized','position',[.1 .1 .5 .4])
saveas(gcf,[saveFolder filesep 'AV_55db_AvgAcrossBR.svg'])
saveas(gcf,[saveFolder filesep 'AV_55db_AvgAcrossBR.jpg'])
savefig([saveFolder filesep 'AV_55db_AvgAcrossBR.fig'])

%% 3W RM ANOVA
azOnly = true;
AVmod_meanPerMouse = NaN(n_animals,3,7,3);
for mouse = 1:n_animals
    AVmod_meanPerMouse(mouse,:,:,:) = squeeze(nanmean(AVmod{mouse},1));
end
if azOnly
    nSpkPos = 5;
    AVmod_meanPerMouse = AVmod_meanPerMouse(:,:,[2,3,1,4,5],:);
else
    nSpkPos = 7 ;
end

y = reshape(AVmod_meanPerMouse,n_animals,3*nSpkPos*3);
% temp  = squeeze(mean(allCell(:,tAna,2:4,2:8,2:4),2));
% y = reshape(temp,nROIs,3*7*3);
varNames = cell(3*nSpkPos*3,1);
for i = 1:3*nSpkPos*3
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
tbiases = array2table(y, 'VariableNames',varNames);

Brightness = cell(3*nSpkPos*3,1);
Loudness = cell(3*nSpkPos*3,1);
SPKpos = cell(3*nSpkPos*3,1);

c1 = cell(1,1); c1{1} = '1'; Loudness(1:nSpkPos*3,1) = c1;
c1 = cell(1,1); c1{1} = '2'; Loudness(nSpkPos*3+1:(nSpkPos*3)*2,1) = c1;% Loudness(16:(5*3)*2,1) = c1;
c1 = cell(1,1); c1{1} = '3'; Loudness((nSpkPos*3)*2+1:(nSpkPos*3)*3,1) = c1;

for i =1:nSpkPos
    for ii = 1:3
        c1 = cell(1,1); c1{1} = num2str(i); c1 = repmat(c1,3,1);
        SPKpos((nSpkPos*3)*(ii-1)+(i-1)*3+1:3*i+(nSpkPos*3)*(ii-1),1) = c1;
    end
end

for i = 1:3*nSpkPos*3
    toFill = mod(i,3);
    if toFill == 0; toFill = 3; end
    stringToFill = num2str(toFill);
    c1 = cell(1,1); c1{1} = stringToFill; Brightness(i,1) = c1;
end

factorNames = {'Brightness','SPKpos', 'Loudness'};
within = table(Brightness, SPKpos, Loudness, 'VariableNames', factorNames);

if azOnly
    rm = fitrm(tbiases,'V1-V45~1','WithinDesign',within);
else
    rm = fitrm(tbiases,'V1-V63~1','WithinDesign',within);
end
[ranovatblb] = ranova(rm, 'WithinModel','Brightness*SPKpos*Loudness');
save([saveFolder filesep 'ThreeWRManova.mat'],'ranovatblb')

% Mrm1 = multcompare(rm,'SPKpos','By','Brightness','ComparisonType','dunn-sidak');
% Mrm2 = multcompare(rm,'SPKpos','By','Loudness','ComparisonType','dunn-sidak');

if doAllcells
%% average magnitude of AVmod, across BR and POS or LOUD and POS - all cells
grandAvg  = squeeze(nanmean(dataALL(:,:,pos_vec,:)));
cmax = max([mean(grandAvg(:,1:5,:),3),squeeze(mean(grandAvg,1))'],[],'all');
figure;
subplot(1,3,1); hold on
imagesc(mean(grandAvg,3),[0 cmax]); colormap gray
xticks(1:5);xlabel('az pos');xticklabels(pos_labels);xlim([.5 5.5])
yticks(1:3);ylabel('BR')
axis square
title('AV')

subplot(1,3,2); hold on
imagesc(squeeze(mean(grandAvg,1))',[0 cmax]);
xticks(1:5);xlabel('az pos');xticklabels(pos_labels);xlim([.5 5.5])
yticks(1:3);ylabel('LOUD');yticklabels(loudness_vec)
set(gca,'YDir','reverse')
axis square
title('AV')

subplot(1,3,3); hold on
temp = cat(1,visTrace{:});
AudGrandAvg  = squeeze(mean(mean(temp(:,tAna,1,2:6,2:4),2)));
imagesc(AudGrandAvg(pos_vec,:)',[0 .15]);
xticks(1:5);xlabel('az pos');xticklabels(pos_labels);xlim([.5 5.5])
yticks(1:3);ylabel('LOUD');yticklabels(loudness_vec)
set(gca,'YDir','reverse')
axis square
title('A only')


%% %---
% LOUD = 3;
figure
for BR=1:3
    for LOUD = 1:3
        data = grandAvg(BR,:,LOUD);
        subplot(3,3,LOUD+3*(BR-1)); hold on
        plot(data)
        xticks(1:5);xlabel('az pos');xticklabels(pos_labels);xlim([.5 5.5])
        title(['BR=' num2str(BR) ', LOUD=' num2str(loudness_vec(LOUD))])
        ylim([-.01 0.05])
    end
end
% sgtitle(['LOUD=' num2str(loudness_vec(LOUD))])

%% 3W RM ANOVA all cells
azOnly = true;
nSpkPos = 5;
AVmod_meanPerCell = dataALL(:,:,pos_vec,:);
nCells = size(AVmod_meanPerCell,1);

y = reshape(AVmod_meanPerCell,nCells,3*nSpkPos*3);

% temp  = squeeze(mean(allCell(:,tAna,2:4,2:8,2:4),2));
% y = reshape(temp,nROIs,3*7*3);
varNames = cell(3*nSpkPos*3,1);
for i = 1:3*nSpkPos*3
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
tbiases = array2table(y, 'VariableNames',varNames);

Brightness = cell(3*nSpkPos*3,1);
Loudness = cell(3*nSpkPos*3,1);
SPKpos = cell(3*nSpkPos*3,1);

c1 = cell(1,1); c1{1} = '1'; Loudness(1:nSpkPos*3,1) = c1;
c1 = cell(1,1); c1{1} = '2'; Loudness(nSpkPos*3+1:(nSpkPos*3)*2,1) = c1;% Loudness(16:(5*3)*2,1) = c1;
c1 = cell(1,1); c1{1} = '3'; Loudness((nSpkPos*3)*2+1:(nSpkPos*3)*3,1) = c1;

for i =1:nSpkPos
    for ii = 1:3
        c1 = cell(1,1); c1{1} = num2str(i); c1 = repmat(c1,3,1);
        SPKpos((nSpkPos*3)*(ii-1)+(i-1)*3+1:3*i+(nSpkPos*3)*(ii-1),1) = c1;
    end
end

for i = 1:3*nSpkPos*3
    toFill = mod(i,3);
    if toFill == 0; toFill = 3; end
    stringToFill = num2str(toFill);
    c1 = cell(1,1); c1{1} = stringToFill; Brightness(i,1) = c1;
end

factorNames = {'Brightness','SPKpos', 'Loudness'};
within = table(Brightness, SPKpos, Loudness, 'VariableNames', factorNames);

if azOnly
    rm = fitrm(tbiases,'V1-V45~1','WithinDesign',within);
else
    rm = fitrm(tbiases,'V1-V63~1','WithinDesign',within);
end
[ranovatblb] = ranova(rm, 'WithinModel','Brightness*SPKpos*Loudness');

Mrm1 = multcompare(rm,'SPKpos','By','Brightness','ComparisonType','dunn-sidak');
Mrm2 = multcompare(rm,'SPKpos','By','Loudness','ComparisonType','dunn-sidak');
end
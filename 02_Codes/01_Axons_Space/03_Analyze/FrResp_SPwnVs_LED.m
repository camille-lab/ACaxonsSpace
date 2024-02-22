%Plots Fig2a
clear all

MainFolder = 'D:\AVspace_final';

dataset = 'LMaxons';           % 'ACaxons', 'LMaxons'
soundStim = '2-20';            %'2-20', '2-80' max fqcy (in kHz). Only if ACaxons
bckgd = 'bl6';                 % 'bl6'. Only if 2-80 kHz sound

metric = 'dFF';
selection_method = 'wilcoxon';  % 'ttest'; 'wilcoxon'; 'bootstrap'
alpha_threshold = 0.01;
amp_threshold = 0;
dFF_threshold = 0.15;

%%  saveFolder config
dt = datestr(now,'yyyymmdd_HHMM');
stimType = [metric soundStim 'kHz' '_' bckgd];
temp = num2str(alpha_threshold); temp2 = strfind(temp,'.'); alpha = temp(temp2+1:end);
temp = num2str(dFF_threshold); temp2 = strfind(temp,'.'); dFF = temp(temp2+1:end);
SaveFolder = [MainFolder filesep 'Plots' filesep dataset filesep stimType filesep ...
    dt '_' selection_method '_' alpha '_' num2str(amp_threshold) filesep ];

if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end

%% load
DataFolder = [MainFolder filesep '01_Data'  filesep '03_ResponsiveData' filesep];
dataname = [dataset '_LED_' metric '_' selection_method '_dot' alpha '_' num2str(amp_threshold) '_dot' dFF];
fprintf(1,['\n loading ' dataname '...'])
LED = load([DataFolder dataname '.mat']);
dataname = [dataset '_SPKwn_' metric '_' selection_method '_dot' alpha '_' num2str(amp_threshold)  '_dot' dFF];
fprintf(1,['\n loading ' dataname '...'])
SPK = load([DataFolder dataname '.mat']);

fprintf(1,'done \n')

disp(LED.RespROIs.info);

%%
n_animals = size(LED.RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,LED.RespROIs.info.data_details);
n_pos = sum(ispos,2);
colors = distinguishable_colors(n_animals);

% % - fraction responsive
figure; subplot(1,5,[1 2]); hold on

% shuffle mean across mice
FrResp_sh = mean(LED.RespROIs.fractionResp.nResp_sh./LED.RespROIs.fractionResp.nROIs,2,'omitnan');
sh_mean = mean(FrResp_sh);
sh_sd = std(FrResp_sh)/sqrt(length(FrResp_sh));

fill([.5 2.5 2.5 .5],[sh_mean-sh_sd sh_mean-sh_sd sh_mean+sh_sd sh_mean+sh_sd],...
    [.5 .5 .5],'edgecolor','none','facealpha',.2)

for i = 1:n_animals
    nROIs = LED.RespROIs.fractionResp.nROIs(i,ispos(i,:));
    frResp_SPK = SPK.RespROIs.fractionResp.nResp(i,ispos(i,:))./nROIs;
    frResp_LED = LED.RespROIs.fractionResp.nResp(i,ispos(i,:))./nROIs;
    
    % % - mouse average
    jitter = rand(1)*0.1-.05;
    plot([1 2]+jitter,[mean(frResp_SPK) mean(frResp_LED)],'-s','color',[.5 .5 .5])
    plot([1 1]+jitter,[mean(frResp_SPK)+(std(frResp_SPK)./sqrt(length(frResp_LED))) mean(frResp_SPK)-(std(frResp_SPK)./sqrt(length(frResp_LED)))],'color',[.5 .5 .5]);
    plot([2 2]+jitter,[mean(frResp_LED)+(std(frResp_LED)./sqrt(length(frResp_LED))) mean(frResp_LED)-(std(frResp_LED)./sqrt(length(frResp_LED)))],'color',[.5 .5 .5]);
end

% % - mean per mice
temp = mean(LED.RespROIs.fractionResp.nResp./LED.RespROIs.fractionResp.nROIs,2,'omitnan');
LED_mean = mean(temp);
LED_sd = std(temp)./sqrt(length(temp));
temp = mean(SPK.RespROIs.fractionResp.nResp./SPK.RespROIs.fractionResp.nROIs,2,'omitnan');
SPK_mean = mean(temp);
SPK_sd = std(temp)./sqrt(length(temp));
plot([2.1 2.3],[LED_mean LED_mean],'k','linewidth',3);plot([2.2 2.2],[LED_mean-LED_sd LED_mean+LED_sd],'k','linewidth',2)
plot([0.7 0.9],[SPK_mean SPK_mean],'k','linewidth',3);plot([0.8 0.8],[SPK_mean-SPK_sd SPK_mean+SPK_sd],'k','linewidth',2)

xticks([1 2]);xticklabels({'SPK','LED'});xtickangle(45);xlim([0 3])
frResp_SPK = SPK.RespROIs.fractionResp.nResp./LED.RespROIs.fractionResp.nROIs;
frResp_SPK = mean(frResp_SPK,2,'omitnan');
frResp_LED = LED.RespROIs.fractionResp.nResp./LED.RespROIs.fractionResp.nROIs;
frResp_LED = mean(frResp_LED,2,'omitnan');
[~,p]=ttest(frResp_SPK,frResp_LED);
title({['p. ttest, p=' num2str(p,3)],'(animal level)'})
ylim([0 .6])

% % - tests
subplot(1,5,3);
[~,p2]=ttest(frResp_SPK,FrResp_sh);
[~,p3]=ttest(frResp_LED,FrResp_sh);
text(0,1,{['SPK vs sh, p = ' num2str(p2,4)] ['LED vs sh, p = ' num2str(p3,4)]},'horizontalalignment','left','verticalalignment','top')
xlim([0 1]);ylim([0 1]);axis off

%% exclude the LED responsive in SPK responsive for V2L axons and vice versa for AC axons
% % -
frResp_SPK_notLED = NaN(n_animals,max(n_pos));
frResp_LED_notSPK = NaN(n_animals,max(n_pos));
FrOverlap = NaN(n_animals,max(n_pos));
nResp_SPK_notLED_2 = NaN(n_animals,max(n_pos));
nResp_LED_notSPK_2 = NaN(n_animals,max(n_pos));
for i = 1:n_animals
    for ii = 1:n_pos(i)
        clear temp;
        nROIs = LED.RespROIs.fractionResp.nROIs(i,ii);
        temp(:,1) = SPK.RespROIs.fractionResp.responsive_SPK{i,ii};
        temp(:,2) = LED.RespROIs.fractionResp.responsive_LED{i,ii};
        
        nResp_SPK_notLED = false(length(temp),1);
        nResp_SPK_notLED(temp(:,1)&~temp(:,2)) = true;
        frResp_SPK_notLED(i,ii) = sum(nResp_SPK_notLED)./nROIs;
        nResp_SPK_notLED_2(i,ii) = sum(nResp_SPK_notLED);
        
        nResp_LED_notSPK = false(length(temp),1);
        nResp_LED_notSPK(temp(:,2)&~temp(:,1)) = true;
        frResp_LED_notSPK(i,ii) = sum(nResp_LED_notSPK)./nROIs;
        nResp_LED_notSPK_2(i,ii) = sum(nResp_LED_notSPK);
        
        nOverlap = false(length(temp),1);
        nOverlap(temp(:,1)&temp(:,2)) = true;
        FrOverlap(i,ii) = sum(nOverlap)./nROIs;
    end
end
SPKnotLED_perMice = mean(frResp_SPK_notLED,2,'omitnan');
% mean(SPKnotLED_perMice)
% std(SPKnotLED_perMice)
LEDnotSPK_perMice = mean(frResp_LED_notSPK,2,'omitnan');
% mean(LEDnotSPK_perMice)
% std(LEDnotSPK_perMice)
LEDandSPK_perMice = mean(FrOverlap,2,'omitnan');
% mean(LEDandSPK_perMice)
% std(LEDandSPK_perMice)
subplot(1,5,4);
pie([mean(LEDnotSPK_perMice) mean(LEDandSPK_perMice) mean(SPKnotLED_perMice),...
    1-(mean(LEDnotSPK_perMice)+mean(SPKnotLED_perMice)+mean(LEDandSPK_perMice))])

sum(nResp_LED_notSPK_2,'all','omitnan')

%% resample test - not shown in Mazo et al., 2024
test1 = cat(1,SPK.RespROIs.fractionResp.responsive_SPK{:});
test2 = cat(1,LED.RespROIs.fractionResp.responsive_LED{:});
test = false(length(test1),1);
test(test1&test2) = true;
FrOverlap = sum(test)/length(test);

% - resample test
for i = 1:1000
    test_sh1 = test1(randperm(length(test1)));
    test_sh2 = test2(randperm(length(test1)));
    test_rand = false(34270,1);
    test_rand(test_sh1&test_sh2) = true;
    result(i) = sum(test_rand)/length(test_rand);
end

subplot(1,5,5); hold on
histogram(result);
yl = ylim;
plot([FrOverlap FrOverlap],yl,'r')

%% save
SaveFolder_FrResp = [SaveFolder filesep 'FrResp' filesep];
if ~exist(SaveFolder_FrResp,'dir');mkdir(SaveFolder_FrResp);end
saveas(gcf,[SaveFolder_FrResp filesep 'FrResp.tif'])
saveas(gcf,[SaveFolder_FrResp filesep 'FrResp.svg'])
savefig([SaveFolder_FrResp filesep 'FrResp'])


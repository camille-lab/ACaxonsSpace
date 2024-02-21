%% LMaxons
Dir = 'C:\Users\camil\Data\Histology\AtlasReg\LMaxons';
animalID ={'CMad85','CMad86','RD10278'};

%% load data 
for i = 1:length(animalID)
  path = dir([Dir filesep animalID{i} '_quantification' filesep 'Quantification' filesep 'Reports' filesep,...
        animalID{i} '_RefAtlasRegions' filesep '*csv']);
    clear value value_area
    for ii = 2:length(path) % first one is god knows what
        opts = detectImportOptions([path(ii).folder filesep path(ii).name]);
        opts.SelectedVariableNames ={'RegionName','Load','RegionArea'};
        data = readtable([path(ii).folder filesep path(ii).name],opts);
        value(:,ii-1) = table2array(data(:,2));
        value_area(:,ii-1) = table2array(data(:,3));
    end
    values_all{i} = sum(value(:,2:end),2);  
    region_area{i} = sum(value_area(:,2:end),2);  
end
regions = table2cell(data(:,1));

%% to check the main areas
% figure
% for i = 1:length(animalID)
%     subplot(1,6,i); hold on
%     [B,I]=sort(values_all{i},'descend');
%     areas = table2cell(data(:,1));
%     areas_sorted = areas(I);
%     bar(B)
%     xticks([1:30]);
%     xticklabels(areas_sorted(1:30)')
%     xtickangle(45)
%     xlim([0 21])
% end

%%
clear values_grouped
LM = false(1328,1);
LI=LM;POR=LM;TEA=LM;VISal=LM;V1=LM;VISpl=LM;
MGB=LM;LGN=LM;LP=LM;
for i = 1:1328
    if contains(regions{i},'Lateral visual') %&& ~contains(areas{i},'layer 1')
        LM(i) = true;
    elseif contains(regions{i},'Laterointermediate') %&& ~contains(areas{i},'layer 1')
        LI(i) = true;
    elseif contains(regions{i},'Postrhinal') %&& ~contains(areas{i},'layer 1')
        POR(i) = true;
    elseif contains(regions{i},'Primary visual')
        V1(i) = true;
    elseif contains(regions{i},'Temporal association') %&& ~contains(areas{i},'layer 1')
        TEA(i) = true;
    elseif contains(regions{i},'Anterolateral') %&& ~contains(areas{i},'layer 1')
        VISal(i) = true;
   elseif contains(regions{i},'Posterolateral') %&& ~contains(areas{i},'layer 1')
        VISpl(i) = true;

    elseif contains(regions{i},'Medial geniculate complex')
        MGB(i) = true;
    elseif contains(regions{i},'Dorsal part of the lateral geniculate complex')
        LGN(i) = true;
    elseif contains(regions{i},'Lateral posterior')
        LP(i)=true;
    end
end
%
for i = 1:length(animalID)
    values_grouped(5,i) = sum(values_all{i}(LM));
    values_grouped(4,i) = sum(values_all{i}(LI));
    values_grouped(2,i) = sum(values_all{i}(POR));
    values_grouped(7,i) = sum(values_all{i}(V1));
    values_grouped(1,i) = sum(values_all{i}(TEA));
    values_grouped(3,i) = sum(values_all{i}(VISal));
    values_grouped(6,i) = sum(values_all{i}(VISpl));
       
    RegionsArea(1,i) = sum(region_area{i}(LM));
    RegionsArea(2,i) = sum(region_area{i}(LI));
    RegionsArea(3,i) = sum(region_area{i}(V1));
    RegionsArea(4,i) = sum(region_area{i}(VISal));
    RegionsArea(5,i) = sum(region_area{i}(TEA));
    RegionsArea(6,i) = sum(region_area{i}(POR));
    
    RegionsArea(10,i) = sum(region_area{i}(MGB));
    RegionsArea(11,i) = sum(region_area{i}(LGN));
    RegionsArea(12,i) = sum(region_area{i}(LP));
end
%


%%
figure; hold  on
p=plot(1:7,values_grouped(1:7,:)./max(values_grouped(1:7,:)),'-');
plot(1:7,mean(values_grouped(1:7,:)./max(values_grouped(1:7,:)),2),'k','linewidth',3)
normMean_acrossAreas = mean(values_grouped(1:7,:)./max(values_grouped(1:7,:)),2);
sd_acrossAreas = std(values_grouped(1:7,:)./max(values_grouped(1:7,:)),[],2)./sqrt(length(animalID));
for i = 1:7
   plot([i i],[normMean_acrossAreas(i)+sd_acrossAreas(i) normMean_acrossAreas(i)-sd_acrossAreas(i)],...
       'k','linewidth',3) 
end
xticks(1:7);
xticklabels({'TEA','VISpor','VISal','LI','LM','VISpl','V1'})
xtickangle(45); xlim([0 10])
ylabel('norm. load');
legend(p,animalID)

set(gcf,'units','normalized','position',[.2 .4 .7 .4])
savefig([Dir filesep 'PerAreas'])
saveas(gcf,[Dir filesep 'PerAreas.svg'])
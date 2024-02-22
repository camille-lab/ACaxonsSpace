%% ACaxons
Dir = 'D:\AVspace_final\01_Data\20_Histology\AtlasReg\ACaxons';
animalID ={'CMad50','CMad56','CMad58','CMad62','CMad65','CMad67','CMad68'};


%% load data 
for i = 1:length(animalID)
%     disp
    path = dir([Dir filesep animalID{i} '_newquantification' filesep 'Quantification' filesep 'Reports' filesep,...
        animalID{i} '_newquantification_RefAtlasRegions' filesep '*csv']);
    if isempty(path)
    path = dir([Dir filesep animalID{i} '_newquantification' filesep 'Quantification' filesep 'Reports' filesep,...
        animalID{i} '_RefAtlasRegions' filesep '*csv']);
    end
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

%% to quickly check the main areas
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
A1 = false(1328,1);
Ap=A1;Av=A1;TEA=A1;VISal=A1;Ad=A1;VISl=A1;PTLp=A1;ECT=A1;MGB=A1;LGN=A1;LP=A1;
for i = 1:1328
    if contains(regions{i},'Primary auditory') %&& ~contains(areas{i},'layer 1')
        A1(i) = true;
    elseif contains(regions{i},'Posterior auditory') %&& ~contains(areas{i},'layer 1')
        Ap(i) = true;
    elseif contains(regions{i},'Ventral auditory') %&& ~contains(areas{i},'layer 1')
        Av(i) = true;
    elseif contains(regions{i},'Dorsal auditory')
        Ad(i) = true;
    elseif contains(regions{i},'Temporal association') %&& ~contains(areas{i},'layer 1')
        TEA(i) = true;
    elseif contains(regions{i},'Anterolateral visual') %&& ~contains(areas{i},'layer 1')
        VISal(i) = true;
    elseif contains(regions{i},'Lateral visual') %&& ~contains(areas{i},'layer 1')
        VISl(i) = true;
    elseif contains(regions{i},'Posterior parietal') %&& ~contains(areas{i},'layer 1')
        PTLp(i) = true;
    elseif contains(regions{i},'Ectorhinal') %&& ~contains(areas{i},'layer 1')
        ECT(i) = true;
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
    values_grouped(4,i) = sum(values_all{i}(A1));
    values_grouped(2,i) = sum(values_all{i}(Ap));
    values_grouped(3,i) = sum(values_all{i}(Av));
    values_grouped(5,i) = sum(values_all{i}(Ad));
    values_grouped(1,i) = sum(values_all{i}(TEA));
    values_grouped(6,i) = sum(values_all{i}(VISal));
    values_grouped(7,i) = sum(values_all{i}(VISl));
    values_grouped(8,i) = sum(values_all{i}(PTLp));
    values_grouped(9,i) = sum(values_all{i}(ECT));
    values_grouped(10,i) = sum(values_all{i}(MGB));
    values_grouped(11,i) = sum(values_all{i}(LGN));
    values_grouped(12,i) = sum(values_all{i}(LP));
    
    RegionsArea(1,i) = sum(region_area{i}(TEA));
    RegionsArea(2,i) = sum(region_area{i}(Ap));
    RegionsArea(3,i) = sum(region_area{i}(Av));
    RegionsArea(4,i) = sum(region_area{i}(A1));
    RegionsArea(5,i) = sum(region_area{i}(Ad));
    RegionsArea(6,i) = sum(region_area{i}(VISal));
    RegionsArea(7,i) = sum(region_area{i}(VISl));
    RegionsArea(8,i) = sum(region_area{i}(PTLp));
    RegionsArea(9,i) = sum(region_area{i}(ECT));
    RegionsArea(10,i) = sum(region_area{i}(MGB));
    RegionsArea(11,i) = sum(region_area{i}(LGN));
    RegionsArea(12,i) = sum(region_area{i}(LP));
end
%


%%
figure;hold  on
p=plot(1:9,values_grouped(1:9,:)./max(values_grouped(1:9,:)),'-');
plot(1:9,mean(values_grouped(1:9,:)./max(values_grouped(1:9,:)),2),'k','linewidth',3)
normMean_acrossAreas = mean(values_grouped(1:9,:)./max(values_grouped(1:9,:)),2);
sd_acrossAreas = std(values_grouped(1:9,:)./max(values_grouped(1:9,:)),[],2)./sqrt(length(animalID));
for i = 1:9
   plot([i i],[normMean_acrossAreas(i)+sd_acrossAreas(i) normMean_acrossAreas(i)-sd_acrossAreas(i)],...
       'k','linewidth',3) 
end
xticks(1:9);
xticklabels({'TEA','AUDp','AUDv','A1','AUDd','VISal','VISl','PTLp','ECT'})
xtickangle(45); xlim([0 10])
ylabel('norm. load');
legend(p,animalID)

set(gcf,'units','normalized','position',[.2 .4 .7 .4])
savefig([Dir filesep 'PerAreas'])
saveas(gcf,[Dir filesep 'PerAreas.svg'])
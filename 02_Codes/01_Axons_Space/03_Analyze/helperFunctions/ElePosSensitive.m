function IsEleMod = ElePosSensitive(RespROIs,SaveFolder)
% function designed to balance out azimuth and elevation sampling for 2W anova
%% some definitions
animalID = RespROIs.info.animalID;
base_window = RespROIs.info.base_window;
n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
timeVect = 0:1/6.0962:7;
colors = distinguishable_colors(n_animals);

IsAzModulated = cell(n_animals,max(n_pos));    
IsEleModulated  = cell(n_animals,max(n_pos)); 
IsInteraction = cell(n_animals,max(n_pos));

AzResampled = [1,3,5];
for i = 1:8
    AzResampled = [AzResampled;AzResampled(i,:)+1];   
end

%% do 2W ANOVA on all ROIs, for all resampled azimuth
for animal = 1:n_animals
    for pos = 1:n_pos(animal)
        if ~isempty(RespROIs.data{animal,pos})
            data = RespROIs.data{animal,pos};
            nROIs = size(data,1);
            if nROIs >= 10
                nRep = size(data,4);
                tBase = false(size(timeVect));
                tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
                tBase = tBase(1:size(data,2));
                try
                    resp_window = RespROIs.info.resp_window{animal,pos};
                catch
                    resp_window = RespROIs.info.resp_window;
                end
                tAna = false(size(timeVect));
                tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
                tAna = tAna(1:size(data,2));
                
                % -- do F0
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
                
                meanResp = squeeze(nanmean(dFF(:,tAna,:,:),2));
                temp = dFF(:,randperm(size(dFF,2)),:,:);
                
                isAzModulated = false(nROIs,size(AzResampled,1));    
                isEleModulated = false(nROIs,size(AzResampled,1));  
                isInteraction = false(nROIs,size(AzResampled,1));
                
                nRep = size(meanResp,3);
                for ROI = 1:nROIs
                    temp = squeeze(meanResp(ROI,:,:));
                    temp = reshape(temp,13,3,nRep);
                    for i = 1:size(AzResampled,1)
                        nBoutons = temp(AzResampled(i,:),:,:);
                        nBoutons = permute(nBoutons,[3 2 1]);
                        anova_data = reshape(nBoutons,nRep*3,3);
                        p = anova2(anova_data,nRep,'off');
%                        
                        if p(1)<0.05
                            isAzModulated(ROI,i) = true;
                        end
                        if p(2)<0.05
                            isEleModulated(ROI,i) = true;
                        end
                        if p(3)<0.05
                            isInteraction(ROI,i) = true;
                        end
                        
                    end
                end % end loop through ROIs
                
                multipleBoutonsPerROIs = true;
                if multipleBoutonsPerROIs
                    % - transform ROIs into boutons
                    nBoutons_max = max(RespROIs.nBoutonsPerROI{animal,pos});
                    for BoutonPerROI = 2:nBoutons_max
                        multipleBoutons = find(RespROIs.nBoutonsPerROI{animal,pos}==BoutonPerROI);
                        isInteraction = [isInteraction ; repmat(isInteraction(multipleBoutons,:),BoutonPerROI-1,1)];
                        isAzModulated  = [isAzModulated ; repmat(isAzModulated(multipleBoutons,:),BoutonPerROI-1,1)];
                        isEleModulated  = [isEleModulated ; repmat(isEleModulated(multipleBoutons,:),BoutonPerROI-1,1)];

                    end
                end
                
                IsInteraction{animal,pos}  = isInteraction; 
                IsAzModulated{animal,pos}  = isAzModulated; 
                IsEleModulated{animal,pos} = isEleModulated;
                
            end
        end
    end % end loop through sessions
end % end loop through mice

IsMod.IsAzModulated  = IsAzModulated;
IsMod.IsEleModulated = IsEleModulated;
IsMod.IsInteraction  = IsInteraction;

%% Plot summary fig
nBoutons = cellfun(@length,IsAzModulated);
for i = 1:n_animals
    for ii = 1:max(n_pos)
        if ~ispos(i,ii)
           IsAzModulated{i,ii} = NaN(1,9); 
           IsEleModulated{i,ii} = NaN(1,9);
           IsInteraction{i,ii} = NaN(1,9);
        end
    end
end
% - per mice

temp1 = cellfun(@sum,IsAzModulated,'UniformOutput',false);
sum_acrossAzBin = cellfun(@mean,temp1);
FractionMod_az = sum(sum_acrossAzBin,2,'omitnan')./sum(nBoutons,2);
temp1 = cellfun(@sum,IsEleModulated,'UniformOutput',false);
sum_acrossEleBin = cellfun(@mean,temp1);
FractionMod_ele = sum(sum_acrossEleBin,2,'omitnan')./sum(nBoutons,2);
temp1 = cellfun(@sum,IsInteraction,'UniformOutput',false);
sum_acrossIntBin = cellfun(@mean,temp1);
FractionMod_int = sum(sum_acrossIntBin,2,'omitnan')./sum(nBoutons,2);


figure; hold on
scatter(ones(n_animals,1),FractionMod_az,...
    'markerfacecolor',[.5 .5 .5],'markeredgecolor','none','markerfacealpha',.3)
scatter(2*ones(n_animals,1),FractionMod_ele,...
    'markerfacecolor',[.5 .5 .5],'markeredgecolor','none','markerfacealpha',.3)
scatter(3*ones(n_animals,1),FractionMod_int,...
    'markerfacecolor',[.5 .5 .5],'markeredgecolor','none','markerfacealpha',.3)

frMod = [FractionMod_az FractionMod_ele FractionMod_int];
plot([1:3],[FractionMod_az FractionMod_ele FractionMod_int]','color',[.5 .5 .5])
plot([0.9:2.9;1.1:3.1],[mean(frMod);mean(frMod)],'k','linewidth',5)
plot([1:3;1:3],[mean(frMod)+std(frMod)./sqrt(n_animals);mean(frMod)-std(frMod)./sqrt(n_animals)],...
    'k','linewidth',2)

ylim([0 1]); xlim([0 4])
xticks([1:3]);xticklabels({'Azimuth','Elevation','Interaction'})
ylabel('Fraction Modulated')

[~,p]=ttest(frMod(:,1),frMod(:,2));
title([RespROIs.info.stim_toAnalyze newline 'p.t-test, p=' num2str(p,3)],'interpreter','none')

if ~isempty(SaveFolder)
    savefig([SaveFolder 'FrResp' RespROIs.info.stim_toAnalyze])
    saveas(gcf,[SaveFolder 'FrResp_' RespROIs.info.stim_toAnalyze '.tif']);
    saveas(gcf,[SaveFolder 'FrResp_' RespROIs.info.stim_toAnalyze '.svg']);
    pause(0.1); close
end
end % end of the function
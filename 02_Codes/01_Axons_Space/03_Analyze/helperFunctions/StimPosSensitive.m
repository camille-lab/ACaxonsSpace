function IsMod = StimPosSensitive(RespROIs,SaveFolder)

animalID = RespROIs.info.animalID;
base_window = RespROIs.info.base_window;
n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
timeVect = 0:1/6.0962:7;
colors = distinguishable_colors(n_animals);

IsSpaceModulated = cell(n_animals,max(n_pos)); 
IsAzModulated = cell(n_animals,max(n_pos)); 
IsEleModulated  = cell(n_animals,max(n_pos));
IsInteraction = cell(n_animals,max(n_pos));

IsAzModulated_sh = cell(n_animals,max(n_pos));    
IsEleModulated_sh  = cell(n_animals,max(n_pos)); 
IsInteraction_sh = cell(n_animals,max(n_pos));

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
            meanResp_sh = squeeze(nanmean(temp(:,tAna,:,:),2));
            
            isSpaceModulated = false(nROIs,1); 
            isAzModulated = false(nROIs,1);
            isEleModulated = false(nROIs,1);
            isInteraction = false(nROIs,1);
            isAzModulated_sh = false(nROIs,1); 
            isEleModulated_sh = false(nROIs,1);
            isInteraction_sh = false(nROIs,1);
            
            nRep = size(meanResp,3);
            for ROI = 1:nROIs
                temp = squeeze(meanResp(ROI,:,:));
                temp = reshape(temp,13,3,nRep);
                temp = permute(temp,[3 2 1]);
                anova_data = reshape(temp,nRep*3,13);
                p = anova2(anova_data,nRep,'off');
                temp = squeeze(meanResp_sh(ROI,:,:));
                temp = reshape(temp,13,3,nRep);
                temp = permute(temp,[3 2 1]);
                anova_data_sh = reshape(temp,nRep*3,13);
                p_sh = anova2(anova_data_sh,nRep,'off');
                if p(1)<0.05
                    isAzModulated(ROI) = true;
                end
                if p(2)<0.05
                    isEleModulated(ROI) = true;
                end
                if p(3)<0.05
                    isInteraction(ROI) = true;
                end
                
                if p_sh(1)<0.05
                    isAzModulated_sh(ROI) = true;
                end
                if p_sh(2)<0.05
                    isEleModulated_sh(ROI) = true;
                end
                if p_sh(3)<0.05
                    isInteraction_sh(ROI) = true;
                end
                
            end % end loop through ROIs
            
            multipleBoutonsPerROIs = true;
            if multipleBoutonsPerROIs
                % - transform ROIs into boutons
                nBoutons_max = max(RespROIs.nBoutonsPerROI{animal,pos});
                for BoutonPerROI = 2:nBoutons_max
                    multipleBoutons = find(RespROIs.nBoutonsPerROI{animal,pos}==BoutonPerROI);
                    isInteraction = [isInteraction ; repmat(isInteraction(multipleBoutons),BoutonPerROI-1,1)];
                    isAzModulated  = [isAzModulated ; repmat(isAzModulated(multipleBoutons),BoutonPerROI-1,1)];
                    isEleModulated  = [isEleModulated ; repmat(isEleModulated(multipleBoutons),BoutonPerROI-1,1)];
                    
                    isInteraction_sh = [isInteraction_sh ; repmat(isInteraction_sh(multipleBoutons),BoutonPerROI-1,1)];
                    isAzModulated_sh  = [isAzModulated_sh ; repmat(isAzModulated_sh(multipleBoutons),BoutonPerROI-1,1)];
                    isEleModulated_sh  = [isEleModulated_sh ; repmat(isEleModulated_sh(multipleBoutons),BoutonPerROI-1,1)];
                end
            end
            
            IsInteraction{animal,pos}  = isInteraction;
            IsAzModulated{animal,pos}  = isAzModulated;   
            IsEleModulated{animal,pos} = isEleModulated;  
            
            IsInteraction_sh{animal,pos}  = isInteraction_sh;
            IsAzModulated_sh{animal,pos}  = isAzModulated_sh;
            IsEleModulated_sh{animal,pos} = isEleModulated_sh;
            end
        end
    end % end loop through sessions
end % end loop through mice

IsMod.IsAzModulated  = IsAzModulated;
IsMod.IsEleModulated = IsEleModulated;
IsMod.IsInteraction  = IsInteraction;

%% Plot summary fig
% - per mice
temp1 = cellfun(@sum,IsAzModulated);
temp2 = cellfun(@length,IsAzModulated);
FractionMod_az = sum(temp1,2)./sum(temp2,2);
temp1 = cellfun(@sum,IsEleModulated);
temp2 = cellfun(@length,IsEleModulated);
FractionMod_ele = sum(temp1,2)./sum(temp2,2);
temp1 = cellfun(@sum,IsInteraction);
temp2 = cellfun(@length,IsInteraction);
FractionMod_int = sum(temp1,2)./sum(temp2,2);

% - shuffle
temp1 = cellfun(@sum,IsAzModulated_sh);
temp2 = cellfun(@length,IsAzModulated_sh);
FractionMod_az_sh = sum(temp1,2)./sum(temp2,2);
temp1 = cellfun(@sum,IsEleModulated_sh);
temp2 = cellfun(@length,IsEleModulated_sh);
FractionMod_ele_sh = sum(temp1,2)./sum(temp2,2);
temp1 = cellfun(@sum,IsInteraction_sh);
temp2 = cellfun(@length,IsInteraction_sh);
FractionMod_int_sh = sum(temp1,2)./sum(temp2,2);
sdAz_sh = std(FractionMod_az_sh)/sqrt(length(FractionMod_az_sh));
sdEle_sh = std(FractionMod_ele_sh)/sqrt(length(FractionMod_ele_sh));
sdInt_sh = std(FractionMod_int_sh)/sqrt(length(FractionMod_int_sh));

figure; hold on
% - shuffles
fill([0.8 1.2 1.2 0.8],[mean(FractionMod_az_sh)-sdAz_sh mean(FractionMod_az_sh)-sdAz_sh flip([mean(FractionMod_az_sh)-sdAz_sh mean(FractionMod_az_sh)-sdAz_sh])],[.5 .5 .5])
fill([1.8 2.2 2.2 1.8],[mean(FractionMod_ele_sh)-sdEle_sh mean(FractionMod_ele_sh)-sdEle_sh flip([mean(FractionMod_ele_sh)-sdEle_sh mean(FractionMod_ele_sh)-sdEle_sh])],[.5 .5 .5])
fill([2.8 3.2 3.2 2.8],[mean(FractionMod_int_sh)-sdInt_sh mean(FractionMod_int_sh)-sdInt_sh flip([mean(FractionMod_int_sh)-sdInt_sh mean(FractionMod_int_sh)-sdInt_sh])],[.5 .5 .5])

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

t = table(frMod(:,1),frMod(:,2),frMod(:,3),...
    'VariableNames',{'g1','g2','g3'});
Meas = table([1:3]','VariableNames',{'groups'});
rm = fitrm(t,'g1-g3~1','WithinDesign',Meas);
[rmTable] = ranova(rm);
p = table2array(rmTable(1,5));
c = multcompare(rm,'groups');

title([RespROIs.info.stim_toAnalyze newline '1W ANOVA,p=' num2str(p,3) ',n=' num2str(n_animals) ' mice'],'interpreter','none')

if ~isempty(SaveFolder)
savefig([SaveFolder 'FrResp' RespROIs.info.stim_toAnalyze])
saveas(gcf,[SaveFolder 'FrResp_' RespROIs.info.stim_toAnalyze '.tif']);
saveas(gcf,[SaveFolder 'FrResp_' RespROIs.info.stim_toAnalyze '.svg']);
pause(0.1); close
end

end % end of the function
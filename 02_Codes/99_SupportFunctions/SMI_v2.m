% can use "SMI_plotALL.m" after this function to plot them together
function SMI = SMI_v2(RespROIs,SaveFolder)

animalID = RespROIs.info.animalID;

base_window = RespROIs.info.base_window;

n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);

smi_persession=cell(n_animals,max(n_pos));
smi_permice = cell(1,n_animals);

timeVect = 0:1/6.0962:7;

colors = distinguishable_colors(n_animals);
figure;
for animal = 1:n_animals
    if n_animals>2
        subplot(3,n_animals,animal); hold on
    else
        subplot(3,3,animal); hold on
    end
    for pos = 1:n_pos(animal)
        try
        data = RespROIs.data{animal,pos};
        catch
            data = RespROIs.OFFresp{animal,pos};
        end
        nROIs = size(data,1);
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
        
        meanResp = squeeze(nanmean(nanmean(dFF(:,tAna,:,:),2),4));
        meanResp_perROI = nanmean(meanResp,2);
        SqrResp = meanResp.^2;
        
        respDiff = (meanResp-meanResp_perROI).^2;
        smi = sum(respDiff,2)./sum(SqrResp,2);
        
        % - % transform ROIs into boutons
        nBoutons_max = max(RespROIs.nBoutonsPerROI{animal,pos});
        for BoutonPerROI = 2:nBoutons_max
            multipleBoutons = find(RespROIs.nBoutonsPerROI{animal,pos}==BoutonPerROI);
            smi = [smi;repmat(smi(multipleBoutons),BoutonPerROI-1,1)];
        end
        
        smi_persession{animal,pos} = smi;
        p1 = cdfplot(smi);
        p1.Color = [.5 .5 .5];
        
    end % end loop through positions
    smi_persession_all = cat(1,smi_persession{animal,:});
    p2=cdfplot(smi_persession_all);
    p2.Color = colors(animal,:);
    p2.LineWidth = 2;
    if animal == 1
        xlabel('SMI')
        ylabel('Cumulative Probability')
    else
        xlabel('SMI')
        ylabel('');yticks([])
    end
    xlim([0 1]);ylim([0 1]);
    title(animalID{animal},'interpreter','none')
    smi_permice{animal} = smi_persession_all;
end % end loop through animals

% -- finish plotting
smi_all = cat(1,smi_persession{:});
if n_animals == 2
    subplot(3,3,4); hold on
elseif n_animals>1
    subplot(3,n_animals,[n_animals+1 n_animals+floor(n_animals/2)]); hold on
else
    subplot(3,1,2); hold on
end
for animal = 1:n_animals
    p2 = cdfplot(smi_permice{animal});
    p2.Color = colors(animal,:);
    p2.LineWidth = 0.5;
end
p2=cdfplot(smi_all);
p2.Color = [0 0 0];
p2.LineWidth = 3;
xlabel('SMI');ylabel('Cumulative Probability')
title('responsive boutons - per animal')

smi_all = cat(1,smi_persession{:});
if n_animals == 2
    subplot(3,3,5); hold on
elseif n_animals>1
    subplot(3,n_animals,[n_animals+floor(n_animals/2)+1 n_animals*2-1]); hold on
else
    subplot(3,1,3); hold on
end
histogram(smi_all,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
xlabel('SMI'); ylabel('n boutons')
title('responsive boutons - all')

% -- legend
if n_animals>2
    subplot(3,n_animals,2*n_animals); hold on
else
    subplot(3,3,6); hold on
end
for animal = 1:n_animals
    scatter(0.2,n_animals-(animal-1),'markerFacecolor',colors(animal,:),'markeredgecolor','none')
    text(0.3,n_animals-(animal-1),[animalID{animal} ', n=' num2str(length(smi_permice{animal}))],'horizontalalignment','left','verticalalignment','middle')
end
scatter(0.2,0,'markerFacecolor',[0 0 0],'markeredgecolor',[0 0 0])
text(0.3,0,['all, n=' num2str(length(smi_all))],'horizontalalignment','left','verticalalignment','middle')
xlim([0 1]); ylim([-1 n_animals+1])
axis off

% -- info
if n_animals>2
    subplot(3,n_animals,[2*n_animals+1 3*n_animals]); hold on
else
    subplot(3,3,[7,8]); hold on
end
xl=xlim; yl=ylim;
dt = datestr(now,'yyyymmdd');
text(xl(1),yl(2),{dt,RespROIs.info.stim_toAnalyze,RespROIs.info.metric},...
    'horizontalalignment','left','verticalalignment','top')
text((xl(2)-xl(1))/2,yl(2),{[RespROIs.info.selection_method,', alpha=', num2str(RespROIs.info.alpha_threshold)]},...
    'horizontalalignment','left','verticalalignment','top')
axis off

set(gcf,'units','normalized','position',[.4 .4 .5 .5])

SMI = smi_persession;

savefig([SaveFolder 'SMI'])
saveas(gcf,[SaveFolder 'SMI.tif']);
end
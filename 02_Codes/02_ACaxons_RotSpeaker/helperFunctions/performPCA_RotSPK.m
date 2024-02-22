function performPCA_RotSPK(RespROIs,az_vector,SaveFolder)
n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
timeVect = 0:1/6.0962:10;

nAzPos = length(az_vector);
if length(az_vector) == 12
    azPos_toUse = true(16,1);
   azPos_toUse(1:4) = false;
   [coeff,score,latent,~,explained] = pca(all_data_norm); 
elseif length(az_vector) == 10
    azPos_toUse = false(16,1);
    azPos_toUse([1:5,7:2:end]) = true;
elseif length(az_vector) == 16
    azPos_toUse = true(16,1);
end
%%
counter = 0;
for animal = 1:n_animals
    for pos = 1:n_pos(animal)
        if ~isempty(RespROIs.data{animal,pos})
            counter = counter+1;
            
            data  = RespROIs.data{animal,pos};
            nROIs = size(data,1);
            if nROIs>10
            try
                resp_window = RespROIs.info.resp_window{animal,pos};
            catch
                resp_window = RespROIs.info.resp_window;
            end
            base_window = RespROIs.info.base_window;
            tAna = false(size(timeVect));
            tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
            tAna = tAna(1:size(data,2));
            tBase = false(size(timeVect));
            tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
            tBase = tBase(1:size(data,2));
            
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
            
            dFF_data{counter} = squeeze(median(mean(dFF(:,tAna,:,:),2),4));
            V1az(counter) = RespROIs.V1az(animal,pos);
            end
        end
    end
end

%% - PCA on concatenated data
all_data = cat(1,dFF_data{:});
max_allData = max(all_data,[],2);
all_data_norm = all_data./max_allData;
% all_data_norm(find(max_allData<0),:)=[];
[coeff,score,latent,~,explained] = pca(all_data_norm(:,azPos_toUse));

cumvar = cumsum(explained);
figure;
for i = 1:nAzPos
    subplot(nAzPos,1,i);
    p=plot(coeff(:,i),'k');
    title(['Cum.var=' num2str(cumvar(i),4)])
    if i == 36
        xlabel('Spk az. pos')
        ylabel('PC coeff.')
    else
        xticks('')
    end
end
lg=legend(p,{'+20','0','-20'},'location','northeastoutside');
title(lg,'elevation')
set(gcf,'units','normalized','position',[.4 0 .2 .9])
sgtitle('PCA')
savefig([SaveFolder 'PCA_Traces'])
saveas(gcf,[SaveFolder 'PCA_Traces.tif']);
pause(0.1); close

figure
for i = 1:nAzPos
    subplot(nAzPos+1,1,i);
    imagesc(coeff(:,i)');
    colormap gray
    axis off
end
subplot(nAzPos+1,1,nAzPos+1);
plot(cumvar)

set(gcf,'units','normalized','position',[.4 0.05 .2 .9])
sgtitle('PCA')
savefig([SaveFolder 'PCA_HeatMap'])
saveas(gcf,[SaveFolder 'PCA_HeatMap.tif']);
saveas(gcf,[SaveFolder 'PCA_HeatMap.svg']);
pause(0.1); close

%% PCA per session
figure; hold on
for i = 1:counter
    data_pca = dFF_data{i};
    max_allData = max(data_pca,[],2);
    all_data_norm = data_pca./max_allData;
    [coeff,score,latent,~,explained] = pca(all_data_norm);
    if size(coeff,2)>10
        for ii = 1:2
            subplot(floor(counter/2)+1,6,ii+(i-1)*3);
            imagesc(coeff(:,ii)');
            colormap gray
        end
        subplot(floor(counter/2)+1,6,3+(i-1)*3);
        plot(cumvar)
        pc1(:,i) =  coeff(:,1);
        pc2(:,i) =  coeff(:,2);
    else
        
    end
end
    figure; hold on
    [~,idx_pc1] = max(pc1);
    max_pc1 = idx_pc1;
    [~,idx_pc2] = max(pc2);
    max_pc2 = idx_pc2;
%     V1az = RespROIs.V1az';
    V1az = V1az(:); V1az(isnan(V1az))=[];
    plot(V1az,max_pc1*10-20,'k*')
    plot(V1az,max_pc2*10-20,'b*')
    plot([-20 100],[-20 100],'k:')
end
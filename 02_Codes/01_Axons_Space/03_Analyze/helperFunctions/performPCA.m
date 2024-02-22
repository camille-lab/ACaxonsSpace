function cumvar_all = performPCA(RespROIs,SaveFolder)
n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
timeVect = 0:1/6.0962:7;

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
[coeff,score,latent,~,explained] = pca(all_data_norm);

cumvar = cumsum(explained);
figure;
for i = 1:39
    subplot(6,7,i);
    p=plot(reshape(coeff(:,i),13,3));
    p(1).Color = [0.8500 0.3250 0.0980 1];
    p(2).Color = [0.8500 0.3250 0.0980 0.7];
    p(3).Color = [0.8500 0.3250 0.0980 .5];
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
set(gcf,'units','normalized','position',[.05 0 .9 .9])
sgtitle('PCA')
savefig([SaveFolder 'PCA_Traces'])
saveas(gcf,[SaveFolder 'PCA_Traces.tif']);
pause(0.1); close

figure
for i = 1:39
    subplot(6,7,i);
    imagesc(reshape(coeff(:,i),13,3)');
    colormap gray
end
subplot(6,7,[41 42]);
plot(cumvar)

set(gcf,'units','normalized','position',[.05 0 .9 .9])
sgtitle('PCA')
savefig([SaveFolder 'PCA_HeatMap'])
saveas(gcf,[SaveFolder 'PCA_HeatMap.tif']);
saveas(gcf,[SaveFolder 'PCA_HeatMap.svg']);
pause(0.1); close

%% - PC1 peak as a function of V1RF
figure; hold on
[~,idx_pc1] = max(pc1);
max_pc1 = mod(idx_pc1,13);
max_pc1(max_pc1==0)=13;
[~,idx_pc2] = max(pc2);
max_pc2 = mod(idx_pc2,13);
max_pc2(max_pc2==0)=13;
%     V1az = RespROIs.V1az';
V1az = V1az(:); V1az(isnan(V1az))=[];
plot(V1az,max_pc1*10-20,'k*')
plot(V1az,max_pc2*10-20,'b*')
plot([-20 100],[-20 100],'k:')

title('PC1 (k) and PC2 (b) best azimuth as a function of V1RF')
savefig([SaveFolder 'PC1peak'])
saveas(gcf,[SaveFolder 'PC1peak.tif']);
saveas(gcf,[SaveFolder 'PC1peak.svg']);

end
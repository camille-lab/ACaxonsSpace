%% Peak Az in boutons vs V1RF -- rotating speaker
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
pearsonR_th = 0.3; 
az_vector = -20:10:100;
ele_vector = 20:-20:-20;
timeVect = 0:1/6.0962:10;

% % - do the calculation
idx_perMouse_az = cell(n_animals,max(n_pos),n_subsampling);
idx_perMouse_ele = cell(n_animals,max(n_pos),n_subsampling);
for mouse = 1:n_animals
    for pos = 1:n_pos(mouse)
        data = cat(1,RespROIs.data{mouse,pos});
        nRep = size(data,4);
        
        % % - analysis time windows
            try
                resp_window = RespROIs.info.resp_window{mouse,pos};
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
        
        % - do a dFF
        F0 = nanmean(data(:,tBase,:,:),2);
        dFF = (data-F0)./F0;
        
        % - transform ROIs into boutons
        nBoutons_max = max(RespROIs.nBoutonsPerROI{mouse,pos});
        for BoutonPerROI = 2:nBoutons_max
            multipleBoutons = RespROIs.nBoutonsPerROI{mouse,pos}==BoutonPerROI;
            dFF = cat(1,dFF,repmat(dFF(multipleBoutons,:,:,:),BoutonPerROI-1,1,1,1));
        end
        nROIs = size(dFF,1);
        
        balanced_Resp = squeeze(nanmean(dFF(:,tAna,:,:),2));
        
        for Rand = 1:n_subsampling
            Xvalidated = false(1,nROIs);
            temp = randperm(nRep,nRep/2);
            halfTrials = false(nRep,1); halfTrials(temp) = true;
            A = mean(balanced_Resp(:,:,halfTrials),3);
            B = mean(balanced_Resp(:,:,~halfTrials),3);
            pearson_R = NaN(1,nROIs);
            for i = 1:nROIs
                A2 = squeeze(A(i,:,:)); A2=A2(:);
                B2 = squeeze(B(i,:,:)); B2=B2(:);
                R_temp = corrcoef(A2,B2);
                pearson_R(i) = R_temp(1,2);
                if pearson_R(i)>pearsonR_th
                    Xvalidated(i) = true;
                end
            end
%             balanced_meanResp = mean(balanced_Resp(:,:,:),3);
            balanced_meanResp = mean(balanced_Resp(Xvalidated,:,:),3);

            temp = reshape(balanced_meanResp,size(balanced_meanResp,1),13,3);
            az_avged = mean(temp,3);
            ele_avged = squeeze(mean(temp,2));
            
            [~,idxMax_az]=max(az_avged,[],2);
            [~,idxMax_ele]=max(ele_avged,[],2);
            
            idx_perMouse_az{mouse,pos,Rand} = idxMax_az;
            idx_perMouse_ele{mouse,pos,Rand} = idxMax_ele;
        end
    end
end

%% distribution - all boutons
N_az = NaN(13,n_subsampling);
N_ele = NaN(3,n_subsampling);
for Rand = 1:n_subsampling
    distrib_az = cat(1,idx_perMouse_az{:,:,Rand});
    N_az(:,Rand) = histcounts(distrib_az,[.5:1:13.5],'Normalization','probability');
    
    distrib_ele = cat(1,idx_perMouse_ele{:,:,Rand});
    N_ele(:,Rand) = histcounts(distrib_ele,[.5:1:3.5],'Normalization','probability');
    
    nROIs(Rand) = length(distrib_ele);
end
% - azimtuh
mean_N_az = median(N_az,2);
CI_N_az = prctile(N_az,[5 95],2);

% - elevation
mean_N_ele = median(N_ele,2);
CI_N_ele = prctile(N_ele,[5 95],2);

%% Plot
figure;
subplot(1,4,[1 3]);hold on
plot([1 13],[1/13 1/13],'k:')
fill([[1:13] [13:-1:1]],[CI_N_az(:,1);flip(CI_N_az(:,2))]','r','edgecolor','none','facealpha',0.5)
plot(mean_N_az,'k')
xlabel('Peak azimtuh');xlim([1,13]);xticks([1:13]);xticklabels(az_vector)
ylabel('Fraction');ylim([0 .25])

subplot(1,4,4);hold on
plot([1 3],[1/3 1/3],'k:')
fill([[1:3] [3:-1:1]],[CI_N_ele(:,1);flip(CI_N_ele(:,2))]','r','edgecolor','none','facealpha',0.5)
plot(mean_N_ele,'k')
xlabel('Peak elevation');xlim([1,3]);xticks([1:3]);xticklabels(ele_vector)
ylim([0 1])


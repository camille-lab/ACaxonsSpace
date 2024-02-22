function peakAz = TuningCurves_v2(RespROIs,XVal,nShuffles,pearsonR_th,saveFolder,datatype)
n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
timeVect = 0:1/6.0962:7;
n_resample = nShuffles;

%% - simply calculate the meanResp across positions and concatenate all the data
data_all = cell(n_animals,max(n_pos));
dFF_all = cell(n_animals,max(n_pos));
nROIs = cell(n_animals,max(n_pos));
for animal = 1:n_animals
    for pos = 1:n_pos(animal)
        data = RespROIs.data{animal,pos};
        nROIs{animal,pos} = size(data,1);
        if nROIs{animal,pos}>1
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
            
            % - do a dFF
            F0 = nanmean(data(:,tBase,:,:),2);
            ROI_wNegBase = 0;
            for ROI = 1:nROIs{animal,pos}
                Fo=squeeze(F0(ROI,:,:,:));
                if any(Fo<0,'all')
                    ROI_wNegBase= ROI_wNegBase+1;
                    minF = min(data(ROI,:,:,:),[],'all');
                    data(ROI,:,:,:) = data(ROI,:,:,:)-minF;
                    F0(ROI,1,:,:) = nanmean(data(ROI,tBase,:,:),2);
                end
            end
            dFF = (data-F0)./F0;
            
            % - transform ROIs into boutons
            nBoutons_max = max(RespROIs.nBoutonsPerROI{animal,pos});
            for BoutonPerROI = 2:nBoutons_max
                multipleBoutons = RespROIs.nBoutonsPerROI{animal,pos}==BoutonPerROI;
                dFF = cat(1,dFF,repmat(dFF(multipleBoutons,:,:,:),BoutonPerROI-1,1,1,1));
            end
            
            data_all{animal,pos} = squeeze(mean(mean(dFF(:,tAna,:,:),2),4));
            
            dFF_all{animal,pos} = dFF;
            
            [~,idx_temp] = max(data_all{animal,pos},[],2);
            idx_all{animal,pos} = idx_temp;
        end
    end
end

% % - concatenate
idx_conc = cat(1,idx_all{:});
idex_az = false(length(idx_conc),39);
for i = 1:39
    idex_az(idx_conc==i,i)=true;
end

%% cross validate (plot the 2nd half, if r2>pearsonR_th)
if XVal
    clear idx_all idex
    xVal_half = cell(n_animals,max(n_pos),n_resample);
    xVal_half_ele = cell(n_animals,max(n_pos),n_resample);
    
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            try
                resp_window = RespROIs.info.resp_window{animal,pos};
            catch
                resp_window = RespROIs.info.resp_window;
            end
            tAna = false(size(timeVect));
            tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
            tAna = tAna(1:size(dFF_all{animal,pos},2));
            
            % % - template
            nRep = size(dFF_all{animal,pos},4);
            for Rand = 1:n_resample
                [template_data,idx] = datasample(dFF_all{animal,pos},nRep/2,4,'Replace',false);
                template_mean = squeeze(mean(mean(template_data(:,tAna,:,:),2),4));
                
                % % - data to plot
                test_idx = true(nRep,1);
                test_idx(idx) = false;
                test_data = dFF_all{animal,pos}(:,:,:,test_idx);
                test_mean = squeeze(mean(mean(test_data(:,tAna,:,:),2),4));
                
                % % - correlation
                nROIs = size(template_mean,1);
                R = NaN(1,nROIs);
                for i = 1:nROIs
                    Rtemp = corrcoef(template_mean(i,:),test_mean(i,:));
                    R(i) = Rtemp(1,2);
                end
                correlatedBoutons = false(1,nROIs);
                correlatedBoutons(R>pearsonR_th) = true;
                
                % % - to plot the other half
                % % -- avg across ele and normalize to Peak(az)
                template_mean_reshaped = reshape(template_mean,size(test_mean,1),13,3);
                template_mean_reshaped_az = mean(template_mean_reshaped,3); % average across ele
                [max_template_mean_az,idx_template_az] = max(template_mean_reshaped_az,[],2);
                
                template_mean_reshaped_ele = squeeze(mean(template_mean_reshaped,2)); % average across az
                [max_template_mean_ele,idx_template_ele] = max(template_mean_reshaped_ele,[],2);
                
                max_template_az = NaN(nROIs,1);
                max_template_ele = NaN(nROIs,1);
                for i = 1:nROIs
                    if  max_template_mean_az(i)>0.15
                        max_template_az(i) = max_template_mean_az(i);
                    else
                        max_template_az(i) = NaN;
                    end
                    if  max_template_mean_ele(i)>0.15
                        max_template_ele(i) = max_template_mean_ele(i);
                    else
                        max_template_ele(i) = NaN;
                    end
                end
                temp = test_mean./max_template_az;
                temp2 = test_mean./max_template_ele;
                
                
                if ~isempty(correlatedBoutons)
                    xVal_half{animal,pos,Rand} = temp(correlatedBoutons,:);
                    idx_az{animal,pos,Rand} = idx_template_az(correlatedBoutons);
                    xVal_half_ele{animal,pos,Rand} = temp2(correlatedBoutons,:);
                    idx_ele{animal,pos,Rand} = idx_template_ele(correlatedBoutons);
                end               
            end % end loop through random subsampling
        end % end loop through positions
    end % end loop through animals
    
    %% concatenate data across sessions and mice
    data_conc_xVal_az = NaN(2000,39,100);
    data_conc_xVal_ele = NaN(2000,39,100);
    idex_az = false(2000,13,100);
    idex_ele = false(2000,3,100);
    for Rand = 1:n_resample
        idx_conc = cat(1,idx_az{:,:,Rand});
        for i = 1:13
            idex_az(idx_conc==i,i,Rand)=true;
        end
        data_conc_xVal_az(1:length(idx_conc),:,Rand) = cat(1,xVal_half{:,:,Rand});
        
        idx_conc_ele = cat(1,idx_ele{:,:,Rand});
        for i = 1:3
            idex_ele(idx_conc_ele==i,i,Rand)=true;
        end
        data_conc_xVal_ele(1:length(idx_conc_ele),:,Rand) = cat(1,xVal_half_ele{:,:,Rand});
    end
       
    %% 1) Plot the elevation-averaged tuning curve
    figure; hold on
    for i = 1:13
        dataToPlot_mean = NaN(13,n_resample);
        for Rand = 1:n_resample
            % % - 1) average across elevation
            dataToPlot = [data_conc_xVal_az(idex_az(:,i,Rand),:,Rand)];
            dataToPlot_reshaped = reshape(dataToPlot,size(dataToPlot,1),13,3);
            dataToPlot_mean(:,Rand) = nanmean(nanmean(dataToPlot_reshaped,3),1);
            
            % % - 2) or plot only the data from the max ele
            %              dataToPlot = [data_conc_xVal(idex_az(:,i,Rand),1:13,Rand);...
            %              data_conc_xVal(idex_az(:,i+13,Rand),14:26,Rand);...
            %              data_conc_xVal(idex_az(:,i+26,Rand),27:39,Rand)];
            %              dataToPlot_mean(:,Rand) = nanmean(dataToPlot,1);
            
            nROIs(i,Rand) = size(dataToPlot,1);
        end
               
        DataToPlot = median(dataToPlot_mean,2,'omitnan');
        CI = prctile(dataToPlot_mean,[5 95],2);
        nROIs_avg(i) = nanmean(nROIs(i,:),2);
        nROIs_sd(i) = std(nROIs(i,:),[],2);
        
        % %  - plot the other half (xval) - ele averaged
        subtightplot(1,13,i); hold on
        plot(DataToPlot,'k');
        fill([1:13 13:-1:1],[CI(:,1);flip(CI(:,2))]',...
            [1 0 0],'edgecolor','none','facealpha',.2)
        text(0,1,['n=' num2str(nROIs_avg(i),3) '+-' num2str(nROIs_sd(i),3)],...
            'horizontalalignment','left','verticalalignment','top','fontsize',6)
        
        if i == 1
            xlim([1 13]); xticks([1 13])
            yticks([0 .5 1]); ylim([0 1.1]);
        else
            xlim([1 13]);xticks([]);yticks([]); ylim([0 1.1]);
        end
        
    end
    
    set(gcf,'units','normalized','position',[.05 .5 .9 .4])
    sgtitle(['Xval, R>' num2str(pearsonR_th) ' | median +- 95%CI, ' num2str(n_resample) ' resample'])
    
    % % - save
    TuningFullSaveName=[saveFolder datatype];
    saveas(gcf,[TuningFullSaveName '_2ndHalf_EleAvg.tif'])
    saveas(gcf,[TuningFullSaveName '_2ndHalf_EleAvg.svg'])
    savefig([TuningFullSaveName '_2ndHalf_EleAvg'])
    
    %% 2) Plot the azimuth-averaged tuning curve
    figure; hold on
    dataToPlot_mean = NaN(3,n_resample);
    plotVector = 3:-1:1;
    
    for i = 1:3
        subplot(1,3,plotVector(i)); hold on
        for Rand = 1:n_resample
            select = any(idex_ele(:,i,Rand),2);
            dataToPlot = data_conc_xVal_ele(select,:,Rand);
            dataToPlot_reshaped = reshape(dataToPlot,size(dataToPlot,1),13,3);
            dataToPlot_mean(:,Rand) = squeeze(nanmean(nanmean(dataToPlot_reshaped,2),1));
            nROIs(i,Rand) = size(dataToPlot,1);
        end
        
        DataToPlot = median(dataToPlot_mean,2,'omitnan');
        CI = prctile(dataToPlot_mean,[5 95],2);
        nROIs_avg(i) = nanmean(nROIs(i,:),2);
        nROIs_sd(i) = std(nROIs(i,:),[],2);
        
        plot(flip(DataToPlot),'k');
        fill([3:-1:1 1:3],[CI(:,1);flip(CI(:,2))]',...
            [1 0 0],'edgecolor','none','facealpha',.2)
        text(0,1,['n=' num2str(nROIs_avg(i),3) '+-' num2str(nROIs_sd(i),3)],...
            'horizontalalignment','left','verticalalignment','top','fontsize',6)
        xlim([1 3]); xticks([1:3])
        yticks([0 .5 1]); ylim([0 1.1]);
    end
    
    set(gcf,'units','normalized','position',[.1 .5 .9 .4])
    sgtitle(['Xval, R>' num2str(pearsonR_th) ' | median +- 95%CI, ' num2str(n_resample) ' resample'])
    
    % % - save
    TuningFullSaveName=[saveFolder datatype];
    saveas(gcf,[TuningFullSaveName '_2ndHalf_AzAvg.tif'])
    saveas(gcf,[TuningFullSaveName '_2ndHalf_AzAvg.svg'])
    savefig([TuningFullSaveName '_2ndHalf_AzAvg'])    
end % end if Xval
end
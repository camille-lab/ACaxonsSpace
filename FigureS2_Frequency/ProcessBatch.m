%%
% Uses F values
% analysis done per mouse on responsive boutons

%% paramterize
SaveFolder = ['C:\Users\camil\AudSpace\Plots\FrequencyTuning' filesep];
dataDir = ['C:\Users\camil\AudSpace\Data\FrequencyTuning' filesep];
animalID = {'CMad58','CMad62','CMad65','CMad67','CMad68'};

response = [3 4];
baseline = [0 1];

% % - to determine who's responsive
apha_resp = 0.01;
resp_magn = 0.1;

% % - to determine who's x-validated
pearsonR_th = 0.3;

%% some definitions
n_mice = length(animalID);

time_vect = 0:1/6:4;
time_vect(25) = [];

tAna = false(size(time_vect));
tAna(time_vect>response(1) & time_vect<=response(2)) = true;
tBase = false(size(time_vect));
tBase(time_vect>=baseline(1) & time_vect<baseline(2)) = true;

%% main loop. concatenate data per mice
data_all = cell(n_mice,3); respROI_all = cell(n_mice,3); correlatedBoutons_all = cell(n_mice,3);
isTuned_all = cell(n_mice,3); 
for mouse = 1:n_mice
    mouseDir = [dataDir animalID{mouse}];
    positions = dir([mouseDir filesep 'dataSorted_' animalID{mouse} '*.mat']);
    for pos = 1:length(positions)
        load([positions(pos).folder filesep positions(pos).name]);
        data = data_sorted.tracesDist;
        nROIs = size(data,1);
        
        % % - dFF
        F0 = mean(data(:,tBase,:,:),2);
        F0 = repmat(F0,1,size(data,2),1,1);
        dFF = (data-F0)./F0;
        dFF = dFF(:,1:24,:,:);
        
        % % - best response
        resp = squeeze(mean(mean(dFF(:,tAna,:,:),2),4));
        [maxResp,maxResp_idx] = max(resp,[],2);
        
        respROI = false(nROIs,1); isTuned = false(nROIs,1);
        for i = 1:nROIs
            bestRespDist = squeeze(mean(dFF(i,tAna,maxResp_idx(i),:),2));
            bestBaseDist = squeeze(mean(dFF(i,tBase,maxResp_idx(i),:),2));
            
            h = ttest(bestRespDist,bestBaseDist,apha_resp);
            
            if h && maxResp(i)>resp_magn
                respROI(i,1) = true;
                
                allResp = squeeze(mean(dFF(i,tAna,:,:),2));
                
                t = table(allResp(1,:)',allResp(2,:)',allResp(3,:)',allResp(4,:)',allResp(5,:)',allResp(6,:)',allResp(7,:)',allResp(8,:)',allResp(9,:)',...
                    allResp(10,:)',allResp(11,:)',allResp(12,:)',allResp(13,:)',allResp(14,:)',allResp(15,:)',allResp(16,:)',allResp(17,:)',allResp(18,:)',...
                    'VariableNames',{'meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9',...
                        'meas10','meas11','meas12','meas13','meas14','meas15','meas16','meas17','meas18'});
                Meas = table(1:18,'VariableNames',{'Measurements'});
                rm = fitrm(t,'meas1-meas18~1');
                ranovatbl = ranova(rm);
                p = table2array(ranovatbl(1,5));

                if p<0.05
                    isTuned(i,1) = true;
                end
            end
            
        end % end loop thru ROIs

     
        % - transform ROIs into boutons
        nBoutons_max = max(data_sorted.nBoutonsPerROI);
        for BoutonPerROI = 2:nBoutons_max
            multipleBoutons = find(data_sorted.nBoutonsPerROI==BoutonPerROI);
            respROI = [respROI;repmat(respROI(multipleBoutons),BoutonPerROI-1,1)];
            isTuned = [isTuned;repmat(isTuned(multipleBoutons),BoutonPerROI-1,1)];
            dFF = [dFF;repmat(dFF(multipleBoutons,:,:,:),BoutonPerROI-1,1,1,1)];
        end
        
        % % - save variables
        data_all{mouse,pos} = dFF;
        respROI_all{mouse,pos} = respROI;
        isTuned_all{mouse,pos} = isTuned;

    end % end loop thru positions
end % end loop thru mice

%% distribution per mice  -- not presented in Mazo et al., 2024
figure
for mouse = 1:n_mice
    data = cat(1,data_all{mouse,:});
    respROI = cat(1,respROI_all{mouse,:});
    isTuned = cat(1,isTuned_all{mouse,:});
    
        data_sel = data(respROI&isTuned,:,:,:);

    AvgResp = squeeze(mean(mean(data_sel(:,tAna,:,:),2),4));
    [~,maxResp_idx] = max(AvgResp,[],2); 
    subplot(1,n_mice+1,mouse)
    histogram(maxResp_idx)
    maxRespIdx_all{mouse} = maxResp_idx;
    title(animalID{mouse})
end
subplot(1,n_mice+1,n_mice+1)
histogram(cat(1,maxRespIdx_all{:}))
title('agregated data')

%% distribution - agregate data only
figure
histogram(cat(1,maxRespIdx_all{:}));
title('agregated data')
xlabel('Freq. (kHz)');xticks([1:5:18]); xticklabels({'2','7','12','17'})

savefig([SaveFolder 'TuningDistrib'])
saveas(gcf,[SaveFolder 'TuningDistrib.tif']);

%% all cells, sorted based on BF -- not presented in Mazo et al., 2024
figure
data = cat(1,data_all{:});
respROI = cat(1,respROI_all{:});
isTuned = cat(1,isTuned_all{:});
meandata = squeeze(mean(mean(data(respROI&isTuned,tAna,:,:),2),4));
[maxResp,maxResp_idx] = max(meandata,[],2); 
for i = 1:size(meandata,1)
meandata_norm(i,:) = meandata(i,:)./maxResp(i);
end
[~,I] = sort(maxResp_idx);
imagesc(meandata_norm(I,:),[0 1])
xlabel('Freq. (kHz)');xticks([1:5:18]); xticklabels({'2','7','12','17'})
colormap gray

%% Tuning curves - averaged across cells w/ same BF
figure
for i = 1:18
subplot(18,18,[1:18]+(i-1)*18); hold on
dataaa = meandata(maxResp_idx==i,:)./meandata(maxResp_idx==i,i);
sem = std(dataaa)./sqrt(size(dataaa,1));
fill([1:18 18:-1:1],[mean(dataaa)+sem flip(mean(dataaa)-sem)],'k','edgecolor','none','facealpha',.2)
plot(mean(dataaa),'k');
temp(i,:) = mean(dataaa);
ylim([0 1]); yticks([0 1]); yticklabels({'',''});
text(18,1,['n=' num2str(size(dataaa,1))],'horizontalalignment','left','verticalalignment','top')
if i ~= 18
    xticklabels('')
else
    xlabel('Freq. (kHz)');xticks([1:5:18]); xticklabels({'2','7','12','17'})
end
end
%
set(gcf,'units','normalized','position',[.4 .05 .3 .9])
savefig([SaveFolder 'TuningCurves'])
saveas(gcf,[SaveFolder 'TuningCurves.tif']);

%% Pie chart - Fraction of tuned cells/total resp. cells
figure
x=sum(cat(1,isTuned_all{:}))/sum(cat(1,respROI_all{:}));
pie([x 1-x])

savefig([SaveFolder 'PieChart'])
saveas(gcf,[SaveFolder 'PieChart.tif']);
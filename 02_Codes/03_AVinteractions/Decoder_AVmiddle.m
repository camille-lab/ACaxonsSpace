%% Used to plot spk pos decoding error in A trials, July 2023

%% parametrize
clearvars -except data_all
resp_window = [1.5 3];

MainFolder = 'C:\Users\camil\AudSpace\Data\AV_Cells\L2_3middle';
animals  = {'CMad76','CMad79','CMad80','CMad81',...         % May-June 21
    'CAMP2','CAMP3','CMad133','CMad135'}; % October 23

loudness_vec = [75 60 45];
brightness_vec = [1 3 5];
pos_vect = -40:20:40;

nTrainingTrials = 6;
nRep = 7;
az_vector = [-40:20:40];
nAzPos = length(az_vector);
nElePos = 3;
nTrialTypes = 1*3+7*3+7*3; % use only one LED, 7 As, 7 AVs
nTrials =nTrialTypes*nRep;

n_perm = 100; % for the shuffle

saveDir = 'C:\Users\camil\AudSpace\Plots\AVcells';

%% Load data
if ~exist('data_all','var')
disp('loading data...');
load([MainFolder filesep 'data_all_v5.mat']);
else
    disp('reusing data from workspace')
end
% SaveFolder = ['C:\Users\camil\Data\Face\SpeakerSpace\' dataset];

%% Select the training and test trials - jacknife leave one out
nDecodedTrials = nRep-nTrainingTrials;
nRuns = nRep-nDecodedTrials+1; % number of loop to decode all the trials

% ----- randomize trials used for trainig
decodedTrials = false(nRuns,nRep);
for run = 1:nRuns
    decodedTrials(run,run) = true;
end
decodedTrials = logical(decodedTrials);
trainingTrials = ~decodedTrials;

%% decoder per se
n_animals = length(animals);
yo = cell(n_animals,size(data_all,2));
yo_AVspkPos = cell(n_animals,size(data_all,2));
yo_AspkPos = cell(n_animals,size(data_all,2));
yo_LOUD = cell(n_animals,size(data_all,2));
yo_LOUD_AV = cell(n_animals,size(data_all,2));
yo_BR_AV = cell(n_animals,size(data_all,2));

yo_AVspkPos_sh = cell(n_animals,size(data_all,2),3,3,n_perm);
yo_AspkPos_sh = cell(n_animals,size(data_all,2),3,3,n_perm);

error_AVspkPos = cell(n_animals,size(data_all,2),3,3);
error_AVspkPos_sh = cell(n_animals,size(data_all,2),3,3,n_perm);

error_AspkPos = cell(n_animals,size(data_all,2),3);
error_AspkPos_sh = cell(n_animals,size(data_all,2),3,n_perm);

for mice = 1:n_animals
    for pos = 1:size(data_all,2)
       if ~isempty(data_all{mice,pos}) 
        data = data_all{mice,pos};
        nROIs = size(data,1);
        
        nFrames = size(data,2);
        x_vect = 0:1/6:3; x_vect=x_vect(1:nFrames);
        tAna = false(size(x_vect));
        tAna(x_vect>=resp_window(1) & x_vect<=resp_window(2))=true;
        tAna = tAna(1:size(data,2));
        
        bestLED = squeeze(mean(mean(mean(data(:,tAna,:,2:4,1,1,:),2),4),7));
        [~,idx_bestLED] = max(bestLED,[],2);
        
        respDist = NaN(nROIs,4,8,4,nRep);
        for i = 1:nROIs
            respDist(i,:,:,:,:) = squeeze(nanmean(data(i,tAna,idx_bestLED(i),:,:,:,:),2));
            % respDist has now 5Ds: ROIs, brightness, SPKID, loudness, rep
        end
        
        %% Built the decoder and the response matrix
        PosteriorMatrix_Brightness = cell(nRuns,length(brightness_vec));
        PosteriorMatrix_Brightness_AV = cell(nRuns,length(brightness_vec));
        PosteriorMatrix_AVspkPos = cell(nRuns,length(az_vector),length(brightness_vec),length(loudness_vec));
        PosteriorMatrix_AspkPos = cell(nRuns,length(az_vector),length(loudness_vec));
        PosteriorMatrix_Loudness = cell(nRuns,length(loudness_vec));
        PosteriorMatrix_Loudness_AV = cell(nRuns,length(loudness_vec));

        PosteriorMatrix_AVspkPos_sh = cell(nRuns,length(az_vector),length(brightness_vec),length(loudness_vec),n_perm);
        PosteriorMatrix_AspkPos_sh = cell(nRuns,length(az_vector),length(loudness_vec),n_perm);

        for iii = 1:nRuns
            
            sdMatrix = std(respDist(:,:,:,:,trainingTrials(iii,:)),[],5,'omitnan');
            RespMatrix = nanmean(respDist(:,:,:,:,trainingTrials(iii,:)),5);
            
            meanResp = respDist(:,:,:,:,decodedTrials(iii,:));
                                    
            % % - decode brightness in V and AV trilas (middle speaker)
            Posterior_perStimType = NaN(1,length(brightness_vec));
            Posterior_perStimType_AV = NaN(length(brightness_vec),length(loudness_vec));
            for ii = 1:length(brightness_vec)
                for j = 1:length(brightness_vec)
                    Posterior_perROI = NaN(1,nROIs);
                    Posterior_perROI_AV = NaN(1,nROIs);                    
                    for i = 1:nROIs
                        Posterior_perROI(i) = -log(sdMatrix(i,j+1,1,1)*sqrt(2*pi)) - (meanResp(i,ii+1,1,1)-RespMatrix(i,j+1,1,1))^2./(2*sdMatrix(i,j+1,1,1)^2);
                    for LOUD = 1:3
                        Posterior_perROI_AV(i) = ...
                             -log(sdMatrix(i,j+1,2,LOUD+1)*sqrt(2*pi)) - (meanResp(i,ii+1,2,LOUD+1)-RespMatrix(i,j+1,2,LOUD+1))^2./(2*sdMatrix(i,j+1,2,LOUD+1)^2);
                    end
                    end
                    Posterior_perStimType(j) = nansum(Posterior_perROI);
                    Posterior_perStimType_AV(j,:) = nansum(Posterior_perROI_AV);
                end
                
                PosteriorMatrix_Brightness{iii,ii} = Posterior_perStimType;
                PosteriorMatrix_Brightness_AV{iii,ii} = Posterior_perStimType_AV;
            end
            
            % % - decode position in AV trials all BR and LOUD
            for BR = 1:3
                for LOUD = 1:3
                    Posterior_perStimType = NaN(1,length(az_vector));
                    for ii = 1:length(az_vector) % all trial types
                        for j = 1:length(az_vector)
                            Posterior_perROI = NaN(1,nROIs);
                            for i = 1:nROIs
                                Posterior_perROI(i) = -log(sdMatrix(i,BR+1,j+1,LOUD+1)*sqrt(2*pi)) - (meanResp(i,BR+1,ii+1,LOUD+1)-RespMatrix(i,BR+1,j+1,LOUD+1))^2./(2*sdMatrix(i,BR+1,j+1,LOUD+1)^2);
                            end
                            Posterior_perStimType(j) = nansum(Posterior_perROI);
                        end
                        PosteriorMatrix_AVspkPos{iii,ii,BR,LOUD} = Posterior_perStimType;
                    end
                end
            end
 
            % % - shuffle --------------
            for perm = 1:n_perm
            temp = randperm(5); k = [1 temp+1 7 8];
            meanResp_sh = meanResp(:,:,k,:);
            for BR = 1:3
                for LOUD = 1:3
                    Posterior_perStimType_sh = NaN(1,length(az_vector));
                    for ii = 1:length(az_vector)
                        for j = 1:length(az_vector)
                            Posterior_perROI_sh = NaN(1,nROIs);
                            for i = 1:nROIs
                                Posterior_perROI_sh(i) = -log(sdMatrix(i,BR+1,j+1,LOUD+1)*sqrt(2*pi)) - (meanResp_sh(i,BR+1,ii+1,LOUD+1)-RespMatrix(i,BR+1,j+1,LOUD+1))^2./(2*sdMatrix(i,BR+1,j+1,LOUD+1)^2);
                            end
                            Posterior_perStimType_sh(j) = nansum(Posterior_perROI_sh);
                        end
                        PosteriorMatrix_AVspkPos_sh{iii,ii,BR,LOUD,perm} = Posterior_perStimType_sh;
                    end
                end
            end
            end
             % % --------------------
             
               % % - decode position in A trials
                for LOUD = 1:3
                    Posterior_perStimType = NaN(1,length(az_vector));
                    for ii = 1:length(az_vector)
                        for j = 1:length(az_vector)
                            Posterior_perROI = NaN(1,nROIs);
                            for i = 1:nROIs
                                Posterior_perROI(i) = -log(sdMatrix(i,1,j+1,LOUD+1)*sqrt(2*pi)) - (meanResp(i,1,ii+1,LOUD+1)-RespMatrix(i,1,j+1,LOUD+1))^2./(2*sdMatrix(i,1,j+1,LOUD+1)^2);
                            end
                            Posterior_perStimType(j) = nansum(Posterior_perROI);
                        end
                        PosteriorMatrix_AspkPos{iii,ii,LOUD} = Posterior_perStimType;
                    end
                end
            
                % % - shuffle --------------
            for perm = 1:n_perm
            temp = randperm(5); k = [1 temp+1 7 8];
            meanResp_sh = meanResp(:,:,k,:);
                for LOUD = 1:3
                    Posterior_perStimType_sh = NaN(1,length(az_vector));
                    for ii = 1:length(az_vector)
                        for j = 1:length(az_vector)
                            Posterior_perROI_sh = NaN(1,nROIs);
                            for i = 1:nROIs
                                Posterior_perROI_sh(i) = -log(sdMatrix(i,1,j+1,LOUD+1)*sqrt(2*pi)) - (meanResp_sh(i,1,ii+1,LOUD+1)-RespMatrix(i,1,j+1,LOUD+1))^2./(2*sdMatrix(i,1,j+1,LOUD+1)^2);
                            end
                            Posterior_perStimType_sh(j) = nansum(Posterior_perROI_sh);
                        end
                        PosteriorMatrix_AspkPos_sh{iii,ii,LOUD,perm} = Posterior_perStimType_sh;
                    end
                end
            end
            
            % % - decode loudness in A and AV trials USING MIDDLE SPEAKER
            Posterior_perStimType = NaN(1,length(loudness_vec));
            Posterior_perStimType_AV = NaN(length(loudness_vec),length(brightness_vec));
            for ii = 1:length(loudness_vec)
                for j = 1:length(loudness_vec)
                    Posterior_perROI = NaN(1,nROIs);
                    Posterior_perROI_AV = NaN(nROIs,length(brightness_vec));
                    for i = 1:nROIs
                        Posterior_perROI(i) = -log(sdMatrix(i,1,2,j+1)*sqrt(2*pi)) - (meanResp(i,1,2,ii+1)-RespMatrix(i,1,2,j+1))^2./(2*sdMatrix(i,1,2,j+1)^2);
                        for BR = 1:length(brightness_vec)
                            Posterior_perROI_AV(i,BR) =...
                                -log(sdMatrix(i,BR+1,2,j+1)*sqrt(2*pi)) - (meanResp(i,BR+1,2,ii+1)-RespMatrix(i,BR+1,2,j+1))^2./(2*sdMatrix(i,BR+1,2,j+1)^2);
                        end
                    end
                    Posterior_perStimType(j) = nansum(Posterior_perROI);
                    Posterior_perStimType_AV(j,:) = squeeze(nansum(Posterior_perROI_AV));
                end
                PosteriorMatrix_Loudness{iii,ii} = Posterior_perStimType;
                PosteriorMatrix_Loudness_AV{iii,ii} = Posterior_perStimType_AV;
            end
            
        end % end loop through runs
        
        %% brightness in V trials
        decoded = NaN(length(brightness_vec),length(brightness_vec));
        decodedAV = NaN(length(brightness_vec),length(brightness_vec),length(loudness_vec));
        for ii = 1:length(brightness_vec)
             % % - AV trials
            temp = cat(1,PosteriorMatrix_Brightness{:,ii});
            [~,decodedID] = max(temp,[],2);
            for i = 1:length(brightness_vec)
                decoded(i,ii) = length(find(decodedID==i));
            end
             % % - AV trials
            temp = cat(3,PosteriorMatrix_Brightness_AV{:,ii});
            for LOUD = 1:length(brightness_vec)
            [~,decodedID] = max(squeeze(temp(:,LOUD,:)),[],1);
            for i = 1:length(loudness_vec)
                decodedAV(i,ii,LOUD) = length(find(decodedID==i));
            end
            end
        end
        yo{mice,pos} = decoded;
        yo_BR_AV{mice,pos} = decodedAV;
        
        %% decode position in AV trials (abs and error)
        ordered_azVec = [3,1,2,4,5];
        for BR = 1:3
            for LOUD = 1:3
                decoded = NaN(length(az_vector),length(az_vector));
                delta = NaN(1,length(az_vector));
                for ii = 1:length(az_vector)
                    temp = cat(1,PosteriorMatrix_AVspkPos{:,ii,BR,LOUD});
                    [~,decodedID] = max(temp,[],2);
                    diff = NaN(1,nRuns);
                    for jj = 1:nRuns
                    diff(jj) = abs(ordered_azVec(decodedID(jj)) - ordered_azVec(ii));
                    end
                    delta(1,ii) = mean(diff);
                    for i = 1:length(az_vector)
                        decoded(i,ii) = length(find(decodedID==i));
                    end
                end
                yo_AVspkPos{mice,pos,BR,LOUD} = decoded;
                error_AVspkPos{mice,pos,BR,LOUD} = delta;
            end
        end
        
        % % - shuffle --------------
        for perm = 1:n_perm
        for BR = 1:3
            for LOUD = 1:3
                decoded_sh = NaN(length(az_vector),length(az_vector));
                delta_sh = NaN(1,length(az_vector));
                for ii = 1:length(az_vector)
                    temp = cat(1,PosteriorMatrix_AVspkPos_sh{:,ii,BR,LOUD,perm});
                    [~,decodedID] = max(temp,[],2);
                    diff = NaN(1,nRuns);
                    for jj = 1:nRuns
                    diff(jj) = abs(ordered_azVec(decodedID(jj)) - ordered_azVec(ii));
                    end
                    delta_sh(1,ii) = mean(diff);
                    for i = 1:length(az_vector)
                        decoded_sh(i,ii) = length(find(decodedID==i));
                    end
                end
                yo_AVspkPos_sh{mice,pos,BR,LOUD,perm} = decoded_sh;
                error_AVspkPos_sh{mice,pos,BR,LOUD,perm} = delta_sh;
            end
        end
        end
        
        %% decode position in A trials (abs and error)
        for LOUD = 1:3
            decoded = NaN(length(az_vector),length(az_vector));
            delta = NaN(1,length(az_vector));
            for ii = 1:length(az_vector)
                temp = cat(1,PosteriorMatrix_AspkPos{:,ii,LOUD});
                [~,decodedID] = max(temp,[],2);
                diff = NaN(1,nRuns);
                for jj = 1:nRuns
                    diff(jj) = abs(ordered_azVec(decodedID(jj)) - ordered_azVec(ii));
                end
                delta(1,ii) = mean(diff);                    
                for i = 1:length(az_vector)
                    decoded(i,ii) = length(find(decodedID==i));
                end
            end
            yo_AspkPos{mice,pos,LOUD} = decoded;
            error_AspkPos{mice,pos,LOUD} = delta;
        end
        
        % % - shuffle --------------
        for perm = 1:n_perm
            for LOUD = 1:3
                decoded_sh = NaN(length(az_vector),length(az_vector));
                delta_sh = NaN(1,length(az_vector));
                for ii = 1:length(az_vector)
                    temp = cat(1,PosteriorMatrix_AspkPos_sh{:,ii,LOUD,perm});
                    [~,decodedID] = max(temp,[],2);
                    diff = NaN(1,nRuns);
                    for jj = 1:nRuns
                        diff(jj) = abs(ordered_azVec(decodedID(jj)) - ordered_azVec(ii));
                    end
                    delta_sh(1,ii) = mean(diff);
                    for i = 1:length(az_vector)
                        decoded_sh(i,ii) = length(find(decodedID==i));
                    end
                end
                yo_AspkPos_sh{mice,pos,LOUD,perm} = decoded_sh;
                error_AspkPos_sh{mice,pos,LOUD,perm} = delta_sh;
            end
        end
        
        %% decode loudness in A vs AV trials USING MIDDLE SPEAKER
        decoded = NaN(length(loudness_vec),length(loudness_vec));
        decodedAV = NaN(length(loudness_vec),length(loudness_vec),length(brightness_vec));
        for ii = 1:length(loudness_vec)
            % % - A trials
            temp = cat(1,PosteriorMatrix_Loudness{:,ii});
            [~,decodedID] = max(temp,[],2);
            for i = 1:length(loudness_vec)
                decoded(i,ii) = length(find(decodedID==i));
            end
            % % - AV trials
            temp = cat(3,PosteriorMatrix_Loudness_AV{:,ii});
            for BR = 1:length(brightness_vec)
            [~,decodedID] = max(squeeze(temp(:,BR,:)),[],1);
            for i = 1:length(loudness_vec)
                decodedAV(i,ii,BR) = length(find(decodedID==i));
            end
            end
        end
                    yo_LOUD{mice,pos} = decoded;
                    yo_LOUD_AV{mice,pos} = decodedAV;
       end
    end % end loop through recording position
end % end loop through mice


%% decode brightness. average across sessions, then mice
% % -- AV
yop = NaN(n_animals,3,3);
for mice = 1:n_animals
    temp = NaN(3,3,2);
    for pos = 1:2
         if ~isempty(yo{mice,pos})
        temp(:,:,pos) = yo{mice,pos};
         end
    end
    yop(mice,:,:) = mean(temp,3);
end
yo_mice = squeeze(mean(yop,1));

% % -- AV
yop_AV = NaN(n_animals,length(brightness_vec),length(brightness_vec));
for mice = 1:n_animals
    temp = NaN(length(brightness_vec),length(brightness_vec),length(loudness_vec),2);
    for pos = 1:2
        if ~isempty(yo_BR_AV{mice,pos})
        temp(:,:,:,pos) = yo_BR_AV{mice,pos};
        end
    end
    yop_AV(mice,:,:) = mean(mean(temp,4),3); % average across positions and brightnesses
end

% % - calculate performance, i.e., number of times where the decoded BR
% was the actual BR
perf_permice = NaN(n_animals,length(loudness_vec));
perf_permice_AV = NaN(n_animals,length(loudness_vec));
% subplot(1,2,1);imagesc(yo_mice);
figure; subplot(1,2,1); hold on
for i = 1:n_animals
    perf = NaN(1,length(loudness_vec));
    perfAV = NaN(1,length(loudness_vec));
    for ii = 1:length(loudness_vec) % brightness
        perf(ii) = squeeze(yop(i,ii,ii))/7;
        perfAV(ii) = squeeze(yop_AV(i,ii,ii))/7;
    end
    perf_permice(i,:)=perf;
    perf_permice_AV(i,:)=perfAV;
%     scatter([1:length(loudness_vec)],perf,'k','filled')
%     scatter([1:length(loudness_vec)],perfAV,'r','filled')

    p1 = plot(perf,'color',[.5 .5 .5]); p1.Color = [.5 .5 .5 0.5];
    p2 = plot(perfAV,'color',[1 0 0]); p2.Color = [1 0 0 0.5];
end
xlim([.5 length(loudness_vec)+.5]); ylim([0 1])
plot([1 3],[1/3 1/3],'k:')
plot([.9 1.9 2.9;1.1 2.1 3.1],[mean(perf_permice);mean(perf_permice)],'k','linewidth',2)
plot([.9 1.9 2.9;1.1 2.1 3.1],...
    [nanmean(perf_permice);nanmean(perf_permice)],'k', 'linewidth',2)
plot([.9 1.9 2.9;1.1 2.1 3.1],...
    [nanmean(perf_permice_AV);nanmean(perf_permice_AV)],'r', 'linewidth',2)

subplot(1,2,2); hold on
plot([.9 1.9;1.1 2.1],...
    [mean(nanmean(perf_permice)) mean(nanmean(perf_permice_AV));mean(nanmean(perf_permice)) mean(nanmean(perf_permice_AV))],'k')
ylim([0 1])
plot([1 2],[.33 .33],'k:')

%% Plot decoded speaker position, AV trials
yop = NaN(n_animals,5,5,3,3);
error = NaN(n_animals,5,3,3);
for BR = 1:3
    for LOUD = 1:3
        for mice = 1:n_animals
            temp = NaN(length(az_vector),length(az_vector),2);
            temp2 = NaN(5,2);
            for pos = 1:2
                if ~isempty(yo_AVspkPos{mice,pos,BR,LOUD})
                    temp(:,:,pos) = yo_AVspkPos{mice,pos,BR,LOUD};
                    temp2(:,pos) = error_AVspkPos{mice,pos,BR,LOUD};
                end
            end
            yop(mice,:,:,BR,LOUD) = nanmean(temp,3);
            error(mice,:,BR,LOUD) = nanmean(temp2,2);
        end
    end
end
yop_mice = squeeze(mean(mean(mean(yop,4),5),1));
figure;
imagesc(yop_mice([2,3,1,4,5],[2,3,1,4,5]));

 % % - shuffle --------------
yop_sh = NaN(n_animals,5,5,3,3,n_perm);
error_sh = NaN(n_animals,5,3,3,n_perm);
for perm = 1:n_perm
for BR = 1:3
    for LOUD = 1:3
        for mice = 1:n_animals
            temp_sh = NaN(length(az_vector),length(az_vector),2);
            temp2 = NaN(5,2);
            for pos = 1:2
                if ~isempty(yo_AVspkPos{mice,pos,BR,LOUD})
                temp_sh(:,:,pos) = yo_AVspkPos_sh{mice,pos,BR,LOUD,perm};
                temp2(:,pos) = error_AVspkPos_sh{mice,pos,BR,LOUD,perm};
                end
            end
            yop_sh(mice,:,:,BR,LOUD,perm) = nanmean(temp_sh,3);
            error_sh(mice,:,BR,LOUD,perm) = nanmean(temp2,2);
        end
    end
end
end

% % - plot all BR and LOUD 
% perf_permice = NaN(4,length(az_vector));
% figure;
% for BR = 1:3
%     for LOUD = 1:3
%         subplot(3,3,LOUD+3*(BR-1)); hold on
%         for i = 1:4
%             for ii = 1:length(az_vector) % brightness
%                 perf(ii) = squeeze(yop(i,ii,ii,BR,LOUD))/7;
%             end
%             perf_permice(i,:)=perf;
%             scatter([1:5],perf)
%             
%             plot(perf_permice(i,:),'color',[.5 .5 .5])
%         end
%         xlim([.5 5.5]); ylim([0 .5])
%         plot([1 5],[1/5 1/5],'k:')
%         plot([.9 1.9 2.9 3.9 4.9;1.1 2.1 3.1 4.1 5.1],[mean(perf_permice);mean(perf_permice)],'k','linewidth',2)
%         [p,~,stats] = anova1(perf_permice,[],'off');
%         title(['BR=' num2str(BR) ',LOUD=' num2str(LOUD) ', p=' num2str(p,4)])
%     end
% end

%% error in decoding spk pos (degrees)
error_allMice = mean(mean(error,3),4);
error_allMice_sh = squeeze(mean(mean(mean(error_sh,3),4),1));
CI_sh = prctile(error_allMice_sh,[5 95],2);

figure; hold on
fill([1:5 5:-1:1],[CI_sh([2,3,1,4,5],1);flip(CI_sh([2,3,1,4,5],2))]*20,'k','edgecolor','none','facealpha',.2)
plot(mean(error_allMice(:,[2,3,1,4,5]))*20,'r','linewidth',2)
for i = 1:n_animals
plot(error_allMice(i,[2,3,1,4,5])*20,'color',[.5 .5 .5])
end
ylim([15 55])
xlim([.5 5.5])
chance(1) = (10+20+30+40)/5;
chance(2) = (10+10+20+30)/5;
chance(3) = (20+10+10+20)/5;
chance(4) = chance(2);
chance(5) = chance(1);
plot(chance*2,'k:')
ylabel('distance btwn decoded and actual speaker position (degrees)')
xlabel('spk-led distance')
xticks(1:5); xticklabels({'-40','-20','0','+20','+40'})
title('AV trials')
ylim([20 45])

saveas(gcf,[saveDir filesep 'DeocderError_AV.tif'])
saveas(gcf,[saveDir filesep 'DeocderError_AV.svg'])
savefig(gcf,[saveDir filesep 'DeocderError_AV'])

%% subtract chance
figure; hold on
toPlot = mean(error_allMice(:,[2,3,1,4,5]))*20 - chance*2;
plot(toPlot,'r','linewidth',2)
toPlot = [mean(error_allMice_sh([2,3,1,4,5],:),2)*20]' - chance*2;
plot(toPlot,'k','linewidth',2)
toPlot = error_allMice_sh([2,3,1,4,5],:).*20 - [chance*2]';
CI_sh = prctile(toPlot,[5 95],2);
fill([1:5 5:-1:1],[CI_sh(:,1);flip(CI_sh(:,2))],'k','edgecolor','none','facealpha',.2)

%%
% % - average across BR and LOUD
perf_permice = NaN(4,length(az_vector),3,3);
for BR = 1:3
    for LOUD = 1:3
        for i = 1:n_animals
            perf = NaN(1,length(az_vector));
            for ii = 1:length(az_vector) % brightness
                perf(ii) = squeeze(yop(i,ii,ii,BR,LOUD))/nRuns;
            end
            perf_permice(i,:,BR,LOUD)=perf;
        end
    end
end

perf_permice_sh = NaN(4,length(az_vector),3,3,n_perm);
for perm = 1:n_perm
    for BR = 1:3
        for LOUD = 1:3
            for i = 1:n_animals
                perf_sh = NaN(1,length(az_vector));
                for ii = 1:length(az_vector) % brightness
                    perf_sh(ii) = squeeze(yop_sh(i,ii,ii,BR,LOUD,perm))/nRuns;
                end
                perf_permice_sh(i,:,BR,LOUD,perm)=perf_sh;
            end
        end
    end
end

% % - plot
figure; hold on
temp = squeeze(mean(mean(mean(perf_permice_sh,3),4),1));
CI = prctile(temp,[5 95],2);

fill([1:5 5:-1:1],[CI([2,3,1,4,5],1);flip(CI([2,3,1,4,5],2))]','k','edgecolor','none','facealpha',.2)
plot(mean(mean(perf_permice(:,[2,3,1,4,5],:,:),3),4)','-','color',[.5 .5 .5])
plot(mean(mean(mean(perf_permice(:,[2,3,1,4,5],:,:),3),4),1),'r-','linewidth',3)
plot([1 5],[1/5 1/5],'k:')
xlim([.5 5.5]); ylim([0 .5])
yticks([0.25 .5]); ylabel('decoding accuracy')
xticks([1:5]); xticklabels({'-40','-20','0','+20','+40'}); xlabel('speaker-led distance')
anovadata = mean(mean(perf_permice(:,[2,3,1,4,5],:,:),3),4);
[p,~,stats] = anova1(anovadata,[],'off');
c = multcompare(stats,'display','off');
xticks([1:5]);xticklabels({'-40','-20','0','20','40'});xlabel('spk-led distance')

anovadata_relative(:,1) = anovadata(:,3);
anovadata_relative(:,2) = mean(anovadata(:,[2,4]),2);
anovadata_relative(:,3) = mean(anovadata(:,[1,5]),2);
[p,~,stats] = anova1(anovadata_relative,[],'off');
c = multcompare(stats,'display','off');

sh(1,:) = temp(3,:);
sh(2,:) = mean(temp([2,4],:));
sh(3,:) = mean(temp([1,5],:));
CI_relative = prctile(sh,[9 95],2);

figure; hold on
fill([1:3 3:-1:1],[CI_relative(:,1);flip(CI_relative(:,2))]','k','edgecolor','none','facealpha',.2)
for i =1:4
plot(anovadata_relative(i,:)','color',[.5 .5 .5])
end
plot(mean(anovadata_relative),'r','linewidth',2)
xticks([1:3]);xticklabels({'0','20','40'});xlabel('relative spk-led distance')
yticks([0 .25 .5]);ylim([0 .5])
xticks([1:3]);xticklabels({'0','20','40'});xlabel('rel. speaker-led distance');xlim([.5 3.5])
plot([1 3],[1/5 1/5],'k:')

%% Plot decoded speaker position, A only trials

% % - calculate decoder result per mice
yop = NaN(n_animals,length(az_vector),length(az_vector),length(loudness_vec));
errorAonly = NaN(n_animals,5,3);

for LOUD = 1:3
    for mice = 1:n_animals
        temp = NaN(length(az_vector),length(az_vector),2);
        temp2 = NaN(5,2);
        for pos = 1:2
            if ~isempty(yo_AspkPos{mice,pos,LOUD})
            temp(:,:,pos) = yo_AspkPos{mice,pos,LOUD};
            temp2(:,pos) = error_AspkPos{mice,pos,LOUD};
            end
        end
        yop(mice,:,:,LOUD) = nanmean(temp,3);
        errorAonly(mice,:,LOUD) = nanmean(temp2,2);
    end
end

% % - plot all BR and LOUD
% perf_permice = NaN(4,length(az_vector));
% figure;
%     for LOUD = 1:3
%         subplot(3,3,LOUD+3*(BR-1)); hold on
%         for i = 1:4
%             for ii = 1:length(az_vector) % brightness
%                 perf(ii) = squeeze(yop(i,ii,ii,BR,LOUD))/7;
%             end
%             perf_permice(i,:)=perf;
%             scatter([1:5],perf)
%             
%             plot(perf_permice(i,:),'color',[.5 .5 .5])
%         end
%         xlim([.5 5.5]); ylim([0 1])
%         plot([1 5],[1/5 1/5],'k:')
%         plot([.9 1.9 2.9 3.9 4.9;1.1 2.1 3.1 4.1 5.1],[mean(perf_permice);mean(perf_permice)],'k','linewidth',2)
%     end

% % - average across BR and LOUD
% colorToUse = [1 0 0];
% alpha = [.5 .3 0];
% markersize = [10 5 3];

% % - shuffle --------------
yop_sh = NaN(n_animals,5,5,3,n_perm);
error_sh = NaN(n_animals,5,3,n_perm);
for perm = 1:n_perm
    for LOUD = 1:3
        for mice = 1:n_animals
            temp_sh = NaN(length(az_vector),length(az_vector),2);
            temp2 = NaN(5,2);
            for pos = 1:2
                if ~isempty(yo_AspkPos{mice,pos,LOUD})
                temp_sh(:,:,pos) = yo_AspkPos_sh{mice,pos,LOUD,perm};
                temp2(:,pos) = error_AspkPos_sh{mice,pos,LOUD,perm};
                end
            end
            yop_sh(mice,:,:,LOUD,perm) = nanmean(temp_sh,3);
            error_sh(mice,:,LOUD,perm) = nanmean(temp2,2);
        end
    end
end

% % - Plot the accuracy
% perf_permice = NaN(4,length(az_vector));
% figure; hold on
%     for LOUD = 1:3
%         for i = 1:4
%             for ii = 1:length(az_vector) % brightness
%                 perf(ii) = squeeze(yop(i,ii,ii,LOUD))/7;
%             end
%             perf_permice(i,:,LOUD)=perf;
%             
% %             p1 = plot(perf_permice(i,:),'o-','MarkerSize',markersize(LOUD));
% %             p1.Color = [1 0 0 alpha(BR)];
%         end
%         
%         %         plot([.9 1.9 2.9 3.9 4.9;1.1 2.1 3.1 4.1 5.1],[mean(perf_permice);mean(perf_permice)],'k','linewidth',2)
%     end
% plot(mean(perf_permice,3)','o-','color',[.5 .5 .5])
% plot(mean(mean(perf_permice,3),1),'k-','linewidth',3)
% plot([1 5],[1/5 1/5],'k:')
% xlim([.5 5.5]); ylim([0 1])
% anovadata = mean(perf_permice,3);
% [p,~,stats] = anova1(anovadata,[],'off');
% c = multcompare(stats,'Display','off');
% 
% anovadata_relative(:,1) = anovadata(:,3);
% anovadata_relative(:,2) = mean(anovadata(:,[2,4]),2);
% anovadata_relative(:,3) = mean(anovadata(:,[1,5]),2);
% [p,~,stats] = anova1(anovadata_relative,[],'off');
% c = multcompare(stats,'Display','off');


% % - plot decoder error (degrees)
error_allMice = mean(errorAonly,3);
error_allMice_sh = squeeze(mean(mean(error_sh,3),1));
CI_sh = prctile(error_allMice_sh,[5 95],2);

figure; hold on
fill([1:5 5:-1:1],[CI_sh([2,3,1,4,5],1);flip(CI_sh([2,3,1,4,5],2))]*20,'k','edgecolor','none','facealpha',.2)
plot(mean(error_allMice(:,[2,3,1,4,5]))*20,'r','linewidth',2)
% ----- plot indiv mouse -------
%     plot(error_allMice(i,[2,3,1,4,5])*20,'color',[.5 .5 .5])
% end
ylim([15 55])
xlim([.5 5.5])
chance(1) = (10+20+30+40)/5;
chance(2) = (10+10+20+30)/5;
chance(3) = (20+10+10+20)/5;
chance(4) = chance(2);
chance(5) = chance(1);
plot(chance*2,'k:')
ylabel('distance btwn decoded and actual speaker position (degrees)');
xlabel('spk-led distance');xticks([1:5]); xticklabels({'-40','-20','0','+20','+40'})
ylim([20 45])

title('A only trials')
saveas(gcf,[saveDir filesep 'DeocderError_Aonly.tif'])
saveas(gcf,[saveDir filesep 'DeocderError_Aonly.svg'])
savefig(gcf,[saveDir filesep 'DeocderError_Aonly'])

%% LOUDNESS in A and AV trials
% % - this is just to average data across pos of the same mouse
% % -- A
yop = NaN(n_animals,length(loudness_vec),length(loudness_vec));
for mice = 1:n_animals
    temp = NaN(length(loudness_vec),length(loudness_vec),2);
    for pos = 1:2
        temp(:,:,pos) = yo_LOUD{mice,pos};
    end
    yop(mice,:,:) = mean(temp,3);
end
yo_mice = squeeze(mean(yop,1));

% % - AV
yop_AV = NaN(n_animals,length(loudness_vec),length(loudness_vec));
for mice = 1:n_animals
    temp = NaN(length(loudness_vec),length(loudness_vec),length(brightness_vec),2);
    for pos = 1:2
        temp(:,:,:,pos) = yo_LOUD_AV{mice,pos};
    end
    yop_AV(mice,:,:) = mean(mean(temp,4),3); % average across positions and brightnesses
end

% % - calculate performance, i.e., number of times where the decoded LOUD
% was the actual LOUD
perf_permice = NaN(n_animals,length(loudness_vec));
perf_permice_AV = NaN(n_animals,length(loudness_vec));
% subplot(1,2,1);imagesc(yo_mice);
figure; subplot(1,2,1); hold on
for i = 1:n_animals
    perf = NaN(1,length(loudness_vec));
    perfAV = NaN(1,length(loudness_vec));
    for ii = 1:length(loudness_vec) % brightness
        perf(ii) = squeeze(yop(i,ii,ii))/7;
        perfAV(ii) = squeeze(yop_AV(i,ii,ii))/7;
    end
    perf_permice(i,:)=perf;
    perf_permice_AV(i,:)=perfAV;
    scatter([1:length(loudness_vec)],perf,'k','filled')
    scatter([1:length(loudness_vec)],perfAV,'r','filled')

    plot(perf,'color',[.5 .5 .5])
    plot(perfAV,'color',[1 0 0])
end
xlim([.5 length(loudness_vec)+.5]); ylim([0 1])
plot([1 3],[1/3 1/3],'k:')
plot([.9 1.9 2.9;1.1 2.1 3.1],[mean(perf_permice);mean(perf_permice)],'k','linewidth',2)

subplot(1,2,2); hold on
plot([.9 1.9;1.1 2.1],...
    [mean(mean(perf_permice)) mean(mean(perf_permice_AV));mean(mean(perf_permice)) mean(mean(perf_permice_AV))],'k')
ylim([0 1])
plot([1 2],[.33 .33],'k:')
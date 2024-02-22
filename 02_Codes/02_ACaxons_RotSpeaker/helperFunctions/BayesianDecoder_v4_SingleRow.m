
function Posterior = BayesianDecoder_v4_SingleRow(data,training_fraction,V1az,V1ele,axonRois,info,name,az_vector,MainFolder)
% Plot the posterior distribution for all stimuli locations
% training_fraction: fraction of trials used for training
%% Parametrize
allPositions = 1; % train the decoder using both az and ele (=1), az only (=2) or "chuncks" (=3)
savefigFlag = 1;
base_window = info.base_window;
resp_window = info.resp_window;

%% Define some variables
% az_vector = [-90:20:-10,0:10:100];
nAzPos = length(az_vector);
azPos_toUse = false(39,1);
azPos_toUse(15:26) = true;

nTrialTypes = nAzPos;
nRep = size(data,4);
nTrials =nTrialTypes*nRep;
nROIs = size(data,1);

timeVect = 0:1/6.0962:10;
tBase = false(size(timeVect));
tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
tBase = tBase(1:size(data,2));
tAna = false(size(timeVect));
try
    tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
catch
    tAna(timeVect>=resp_window{1,1}(1) & timeVect<=resp_window{1,1}(2))=true;
end
tAna = tAna(1:size(data,2));

%% Use axonized data
if ~isempty(axonRois)
    for i = 1:length(axonRois) % loop through positions
        axonRois{1,i}(axonRois{1,i}==0)=[];
        if length(axonRois{1,i})>1
            meandFF = mean(mean(mean(data(axonRois{1,i},:,:,:),2),3),4);
            [~,idx_maxdFF] = max(meandFF);
            RepROI(i) = axonRois{1,i}(idx_maxdFF); % use the one that has overall the highest fluo
        else
            RepROI(i) = axonRois{1,i};
        end
    end
    selected_RepRois = false(1,nROIs);
    selected_RepRois(RepROI)=true;
else
    selected_RepRois = true(1,nROIs);
end
nROIs = sum(selected_RepRois);

%%
if sum(selected_RepRois)>=10
    
    %% Select the training and test trials
    nTrainingTrials = training_fraction*nRep;
    nDecodedTrials = round((1-training_fraction)*nTrials);
    if nRep < 20
        %     keyboard
        nTrainingTrials = 13;
        nDecodedTrials = round((1-0.8125)*nTrials);
    end
    nDecodedTrials_perTrialtype = nDecodedTrials/nTrialTypes;
    [~,nRuns] = rat(training_fraction); % number of loop to decode all the trials
    
    V1ele = -V1ele; % because image and plot objects have inversed conventions
    
    % ----- randomize trials used for trainig
    decodedTrials = zeros(nRuns,nRep);
    decodedTrials_order = randperm(nRep);
    c = 0;
    for run = 1:nRuns
        decodedTrials(run,decodedTrials_order(1+c:nDecodedTrials_perTrialtype+c)) = true;
        c = c+nDecodedTrials_perTrialtype;
    end
    decodedTrials = logical(decodedTrials);
    trainingTrials = ~decodedTrials;
    
    %% Built the decoder and the response matrix
    for iii = 1:nRuns % because 80% of the trials are used for the training 4/5th
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
        respDist = squeeze(mean(dFF(:,tAna,azPos_toUse,:),2));
        
        sdMatrix = squeeze(std(respDist(selected_RepRois,:,trainingTrials(iii,:)),[],3,'omitnan'));
        RespMatrix = squeeze(nanmean(respDist(selected_RepRois,:,trainingTrials(iii,:)),3));
        meanResp = respDist(selected_RepRois,:,decodedTrials(iii,:));
        meanResp = reshape(meanResp,nROIs,nDecodedTrials);
        
        %%%%%%%%%%%%%%%%%%%%%%% SHUFFLE Axon label %%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:nTrialTypes
            for j = 1:nTrainingTrials
                shuffle_trainingTrials(:,i,j) = respDist(randperm(nROIs),i,j);
            end
        end
        %     shuffle_trainingTrials = reshape(shuffle_trainingTrials,nROIs,39,16);
        RespMatrix_shuffled = squeeze(nanmean(shuffle_trainingTrials,3));
        sdMatrix_shuffled = squeeze(std(shuffle_trainingTrials,[],3));
        for i = 1:nDecodedTrials
            shuffle_meanResp(:,i) = meanResp(randperm(size(meanResp,1)),i);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 BigMean = mean(mean(meanResp,1),2);
        for ii = 1:nDecodedTrials
            for j = 1:nTrialTypes
                for i = 1:nROIs
                    Posterior_perROI(i) = -log(sdMatrix(i,j)*sqrt(2*pi)) - (meanResp(i,ii)-RespMatrix(i,j))^2./(2*sdMatrix(i,j)^2);
                    Posterior_perROI_shuffled(i) = -log(sdMatrix_shuffled(i,j)*sqrt(2*pi)) - (shuffle_meanResp(i,ii)-RespMatrix_shuffled(i,j))^2./(2*sdMatrix_shuffled(i,j)^2);
                end
                Posterior_perStimType(j) = nansum(Posterior_perROI);
                Posterior_perStimType_shuffled(j) = nansum(Posterior_perROI_shuffled);
            end
            PosteriorMatrix{iii}(:,ii) = Posterior_perStimType;
            PosteriorMatrix_shuffled{iii}(:,ii) = Posterior_perStimType_shuffled;
        end
    end
    
    
    %% Stack into trials of same type
    if allPositions == 1
        for iii = 1:nRuns
            PosteriorMatrix_perStimType{iii} = reshape(PosteriorMatrix{iii},...
                nTrialTypes,nTrialTypes,nDecodedTrials_perTrialtype);
            PosteriorMatrix_perStimType_shuffled{iii} = reshape(PosteriorMatrix_shuffled{iii},...
                nTrialTypes,nTrialTypes,nDecodedTrials_perTrialtype);
        end
    elseif allPositions == 2
        for iii = 1:nRuns
            PosteriorMatrix_perStimType{iii} = reshape(PosteriorMatrix{iii},...
                nAzPos,nAzPos,nDecodedTrials_perTrialtype*3);
            PosteriorMatrix_perStimType_shuffled{iii} = reshape(PosteriorMatrix_shuffled{iii},...
                nAzPos,nAzPos,nDecodedTrials_perTrialtype*3);
        end
    else
        for iii = 1:nRuns
            PosteriorMatrix_perStimType{iii} = reshape(PosteriorMatrix{iii},...
                nPos,nPos,nDecodedTrials_perTrialtype*13);
            PosteriorMatrix_perStimType_shuffled{iii} = reshape(PosteriorMatrix_shuffled{iii},...
                nPos,nPos,nDecodedTrials_perTrialtype*13);
        end
    end
    c = cat(3,PosteriorMatrix_perStimType{:}); % concatenate as many matrices as their are in the cell array "PosteriorMatrix_perStimType"
    avg_PerStimType = squeeze(median(c,3,'omitnan')); % likelihood are too different to take the mean
    d = cat(3,PosteriorMatrix_perStimType_shuffled{:}); % concatenate as many matrices as their are in the cell array "PosteriorMatrix_perStimType"
    avg_PerStimType_shuffle = squeeze(nanmean(d,3));
    
    %% decoding the az position
    for iii = 1:nRuns
        for j = 1:nTrialTypes
            for i = 1:nDecodedTrials_perTrialtype
                [~, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
                [~, I2] = max(PosteriorMatrix_perStimType_shuffled{iii}(:,j,i));
                if  I==j ||  I==j-1 || I==j+1
                    DecodedAz(iii,j,i) = 1;
                else
                    DecodedAz(iii,j,i) = 0;
                end
                if I2==j || I2==j+1 || I2==j-1
                    DecodedAz_shuffle(iii,j,i) = 1;
                else
                    DecodedAz_shuffle(iii,j,i) = 0;
                end
            end
        end
    end
    
    b = squeeze(sum(DecodedAz,3));
    b2 = squeeze(sum(DecodedAz_shuffle,3));
    d = mean(b,1)./nDecodedTrials_perTrialtype;
    collapsed(:,1) = d;
    collapsed(:,2) = az_vector;
    collapsed(:,2) = abs(collapsed(:,2)-V1az);
    [delta_az, II]= sort(collapsed(:,2));
    collapsed2(:,2) = collapsed(II,1);
    collapsed2(:,1) = delta_az;
    
    d2 = mean(b2,1)./nDecodedTrials_perTrialtype;
    collapsed_shuffled(:,1) = d2;
    collapsed2_shuffle(:,2) = collapsed_shuffled(II,1);
    collapsed2_shuffle(:,1) = delta_az;
    
    if ~isnan(V1az)
        binnum = ceil(collapsed2(:,1)/10);
        if binnum(1) == 0; binnum=binnum+1;end
        binsum = accumarray(binnum(:), collapsed2(:,2),[],@mean);
        binsum_shuffle = accumarray(binnum(:), collapsed2_shuffle(:,2),[],@mean);
        
        DecoderFig = figure('Name','DecoderFig');
        h=subplot(2,3,5);
        set(h,'position',[0.5 0.1100 0.3 0.3412]); hold on
        
        p1=plot(((1:size(binsum,1))-1)*10, binsum,'o-','Color',[1 0 0]);
        p2=plot(((1:size(binsum_shuffle,1))-1)*10, binsum_shuffle,'o-','Color',[0 0 0]);
        
        collapsed_sd(:,1) = std(b,[],1)./(sqrt(size(b,1))*nDecodedTrials_perTrialtype);
        collapsed_sd(:,2) = collapsed(:,2);
        collapsed2_sd(:,2) = collapsed_sd(II,1);
        collapsed2_sd(:,1) = delta_az;
        binnum_sd = ceil(collapsed2_sd(:,1)/10);
        if binnum_sd(1) == 0; binnum_sd=binnum_sd+1;end
        binsum_sd = accumarray(binnum_sd(:), collapsed2_sd(:,2),[],@mean);
        
        p3=plot(((1:size(binsum,1))-1)*10, binsum+binsum_sd,'--','Color',[0.75 0.75 0.75]);
        plot(((1:size(binsum,1))-1)*10, binsum-binsum_sd,'--','Color',[0.75 0.75 0.75]);
        p4=plot([0 120],[3/13 3/13],'k:');
        
        ylim([0 1]); yticks([])
        xlim([-5 125]);xlabel('Degrees away from V1 azimuth')
        lgd=legend([p1,p2,p3,p4],{'Data','Shuffle','sem ac. runs','chance'});
        lgd.Position = [0.8 0.2 0.19 0.14];
        lgd.FontSize = 7;
    end
    
    if ~isempty(MainFolder)
        if isnan(V1az)
            DecoderFig = figure('Name','DecoderFig');
        end
        h=subplot(2,3,4); hold on
        set(h,'position',[0.13 0.1100 0.3 0.3412]); hold on
        plot(az_vector,collapsed(:,1),'o-','Color',[1 0 0]);
        plot(az_vector,collapsed_shuffled,'o-','Color',[0 0 0]);
        xlim([-100 110]); xlabel('Speaker Az position');ylabel('Decoding accuracy')
        ylim([0 1]);
        yl = ylim;yticks([0 .5 1])
        plot([V1az V1az],[yl(1) yl(2)],'k:')
        text(V1az,1,'\itV1az',...
            'horizontalalignment','left','verticalalignment','top','FontSize',6)
        plot([az_vector(1) az_vector(end)],[3/nAzPos 3/nAzPos],'k:');
        ylabel('Decoding Accuracy')
        xlabel('Stimulus Position')
    end
   
      
    %% save
    Posterior.collapsed          = collapsed;
    Posterior.collapsed_shuffled = collapsed_shuffled;
    if ~contains(name,'pos') %&& ~strcmp(name,'MickeyMouse')
        Posterior.delta_decod2actual      = delta_decod2actual;
        Posterior.delta_decodShuff2actual = delta_decodShuff2actual;
        Posterior.delta_decod2actual_abs      = delta_decod2actual_abs;
        Posterior.delta_decodShuff2actual_abs = delta_decodShuff2actual_abs;
        
        Posterior.deltaDecoded2Actual = deltaDecoded2Actual;
        
        Posterior.avg_PerStimType = avg_PerStimType;
    elseif ~strcmp(name,'MickeyMouse')
        Posterior.binsum         = binsum;
        Posterior.binsum_shuffle = binsum_shuffle;
    end
else
    disp([name ' skipped'])
    Posterior=[];
end
end
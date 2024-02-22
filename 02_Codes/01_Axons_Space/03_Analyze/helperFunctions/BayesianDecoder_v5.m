function Posterior = BayesianDecoder_v5(data,training_fraction,V1az,V1ele,axonRois,info,name,dataset,MainFolder)
% Plot the posterior distribution for all stimuli locations

%% Parametrize
allPositions = 1; % train the decoder using both az and ele (=1), az only (=2) or "chuncks" (=3). 1 used in Mazo et al., 2024
plotMatrix =  0;  % the matrix of all likelihood for all stim positions
savefigFlag = 1;

dataType = info.metric; % 'dFF' vs 'spikes' or 'Fval'
base_window = info.base_window;
resp_window = info.resp_window;

%% Define some variables
nAzPos = 13;
nElePos = 3;
nTrialTypes = 39;
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
if sum(selected_RepRois)>=10 % if at least 10 axons
    
    %% Select the training and test trials
    nTrainingTrials = training_fraction*nRep;
    nDecodedTrials = round((1-training_fraction)*nTrials);
    if nRep < 20
        nTrainingTrials = 13;
        nDecodedTrials = round((1-0.8125)*nTrials);
    end
    nDecodedTrials_perTrialtype = nDecodedTrials/nTrialTypes;
    [~,nRuns] = rat(training_fraction); % number of loop to decode all the trials
    
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
    if allPositions == 1
        for iii = 1:nRuns % because 80% of the trials are used for the training (4/5th)
            
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
            respDist = squeeze(mean(dFF(:,tAna,:,:),2));
            
            sdMatrix = squeeze(std(respDist(selected_RepRois,:,trainingTrials(iii,:)),[],3,'omitnan'));
            RespMatrix = squeeze(nanmean(respDist(selected_RepRois,:,trainingTrials(iii,:)),3));
            meanResp = respDist(selected_RepRois,:,decodedTrials(iii,:));
            
            meanResp = reshape(meanResp,nROIs,nDecodedTrials);
            
            %% %%%%%%%%%%%%%%%%%%%%% SHUFFLE Axon label %%%%%%%%%%%%%%%%%%%%%%%%
            %             for i = 1:nTrialTypes
            %                 for j = 1:nTrainingTrials
            %                     shuffle_trainingTrials(:,i,j) = respDist(randperm(nROIs),i,j);
            %                 end
            %             end
            %             %     shuffle_trainingTrials = reshape(shuffle_trainingTrials,nROIs,39,16);
            %             RespMatrix_shuffled = squeeze(nanmean(shuffle_trainingTrials,3));
            %             sdMatrix_shuffled = squeeze(std(shuffle_trainingTrials,[],3));
            %             for i = 1:nDecodedTrials
            %                 shuffle_meanResp(:,i) = meanResp(randperm(size(meanResp,1)),i);
            %             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% shuffle trial label -- used in Mazo et al., 2024
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%% SHUFFLE trial label %%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:nROIs
                for j = 1:nTrainingTrials
                    shuffle_trainingTrials(i,:,j) = respDist(i,randperm(nTrialTypes),j);
                end
            end
            
            RespMatrix_shuffled = squeeze(nanmean(shuffle_trainingTrials,3));
            sdMatrix_shuffled = squeeze(std(shuffle_trainingTrials,[],3));
            for i = 1:nROIs
                shuffle_meanResp(i,:) = meanResp(i,randperm(nDecodedTrials));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
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
    elseif allPositions == 2
        for iii = 1:nRuns % because 80% of the trials are used for the training 4/5th
            try
                respDist = data_selected.ONresp.respDist;
            catch
                respDist = data_selected.ONresp.RespDist;
                keyboard
            end
            RespMatrix = squeeze(nanmean(reshape(respDist(:,:,trainingTrials(iii,:)),nROIs,13,3*nTrainingTrials),3));
            sdMatrix = squeeze(std(reshape(respDist(:,:,trainingTrials(iii,:)),nROIs,13,3*nTrainingTrials),[],3,'omitnan'));
            
            meanResp = respDist(:,:,decodedTrials(iii,:));
            meanResp = reshape(meanResp,nROIs,nDecodedTrials);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% SHUFFLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:nTrialTypes
                for j = 1:nTrainingTrials
                    shuffle_trainingTrials(:,i,j) = data_selected.ONresp.respDist(randperm(nROIs),i,j);
                end
            end
            %     shuffle_trainingTrials = reshape(shuffle_trainingTrials,nROIs,39,16);
            RespMatrix_shuffled = squeeze(nanmean(shuffle_trainingTrials,3));
            sdMatrix_shuffled = squeeze(std(shuffle_trainingTrials,[],3));
            for i = 1:nDecodedTrials
                shuffle_meanResp(:,i) = meanResp(randperm(size(meanResp,1)),i);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for ii = 1:nDecodedTrials
                for j = 1:nAzPos
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
    elseif allPositions == 3
        nPos = 3;
        for iii = 1:nRuns % because 80% of the trials are used for the training 4/5th
            try
                respDist = data_selected.ONresp.respDist;
            catch
                respDist = data_selected.ONresp.RespDist;
            end
            
            temp = squeeze(nanmean(reshape(respDist(:,:,trainingTrials(iii,:)),nROIs,13,3*nTrainingTrials),3));
            RespMatrix(:,1) = mean(temp(:,1:4),2);
            RespMatrix(:,2) = mean(temp(:,5:9),2);
            RespMatrix(:,3) = mean(temp(:,10:13),2);
            temp = reshape(respDist(:,:,trainingTrials(iii,:)),nROIs,13,3*nTrainingTrials);
            temp2 = reshape(temp(:,1:4,:),nROIs,4*size(temp,3));
            temp3 = reshape(temp(:,5:9,:),nROIs,5*size(temp,3));
            temp4 = reshape(temp(:,10:13,:),nROIs,4*size(temp,3));
            
            sdMatrix(:,1) = std(temp2,[],2);
            sdMatrix(:,2) = std(temp3,[],2);
            sdMatrix(:,3) = std(temp4,[],2);
            
            meanResp = respDist(:,:,decodedTrials(iii,:));
            meanResp = reshape(meanResp,nROIs,nDecodedTrials);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% SHUFFLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:nTrialTypes
                for j = 1:nTrainingTrials
                    shuffle_trainingTrials(:,i,j) = data_selected.ONresp.respDist(randperm(nROIs),i,j);
                end
            end
            %     shuffle_trainingTrials = reshape(shuffle_trainingTrials,nROIs,39,16);
            RespMatrix_shuffled = squeeze(nanmean(shuffle_trainingTrials,3));
            sdMatrix_shuffled = squeeze(std(shuffle_trainingTrials,[],3));
            for i = 1:nDecodedTrials
                shuffle_meanResp(:,i) = meanResp(randperm(size(meanResp,1)),i);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for ii = 1:nDecodedTrials
                for j = 1:nPos
                    for i = 1:nROIs
                        Posterior_perROI(i) = -log(sdMatrix(i,j)*sqrt(2*pi)) - (meanResp(i,ii)-RespMatrix(i,j))^2./(2*sdMatrix(i,j)^2);
                        Posterior_perROI_shuffled(i) = -log(sdMatrix_shuffled(i,j)*sqrt(2*pi)) - (shuffle_meanResp(i,ii)-RespMatrix_shuffled(i,j))^2./(2*sdMatrix_shuffled(i,j)^2);
                    end
                    Posterior_perStimType(j) = nansum(Posterior_perROI);
                    Posterior_perStimType_shuffled(j) = nansum(Posterior_perROI_shuffled);
                end
                %         meanPosterior_perStimType = mean(Posterior_perStimType);
                %         Posterior_perStimType = Posterior_perStimType- meanPosterior_perStimType;
                PosteriorMatrix{iii}(:,ii) = Posterior_perStimType;
                PosteriorMatrix_shuffled{iii}(:,ii) = Posterior_perStimType_shuffled;
            end
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
    c = cat(3,PosteriorMatrix_perStimType{:});          % concatenate as many matrices as their are in the cell array "PosteriorMatrix_perStimType"
    avg_PerStimType = squeeze(median(c,3,'omitnan'));   % likelihood are too different to take the mean
    d = cat(3,PosteriorMatrix_perStimType_shuffled{:});
    avg_PerStimType_shuffle = squeeze(nanmean(d,3));
    
    %% Spatial aspect of decoder -- this section could be combined with the next one "decoding the az position"
    % -- Spatial distribution of posteriors
    if isempty(MainFolder)
        savefigFlag = 0;
    end
    if savefigFlag
        MatrixFig = figure;
        for i = 1:39
            test = reshape(avg_PerStimType(:,i),13,3);
            [~,I]=max(test,[],'all','linear');
            subtightplot(3,13,i);hold on;
            imagesc(test');colormap gray
            x_ax = mod(i,13); if x_ax ==0; x_ax = 13;end
            if i<14; y_ax = 1; elseif i<27; y_ax=2; else; y_ax=3;end
            plot([x_ax-.5 x_ax+.5 x_ax+.5 x_ax-.5 x_ax-.5],[y_ax-.5 y_ax-.5 y_ax+.5 y_ax+.5 y_ax-.5],'r-')
            decodedPos_x = mod(I,13); if decodedPos_x ==0; decodedPos_x=13;end
            if I<14; decodedPos_y = 1; elseif I<27; decodedPos_y=2; else; decodedPos_y=3;end
            plot(decodedPos_x,decodedPos_y,'o','MarkerSize',1)
            set(gca,'ydir','reverse')
            axis equal;   ylim([0.5 3.5]); xlim([0.5 13.5])
            xticks(''); yticks('')
        end
        set(gcf, 'units', 'normalized','position',[.05 .4 .9 .5])
        sgtitle(['Posterior probabilities, avged across sessions' newline name ' - ' dataset],'interpreter','none')
        MatrixFigFullSaveName=[MainFolder name '_MatrixFig'];
        set(gcf,'units','normalized','position',[0.05 0.4 0.9 0.5])
        saveas(gcf,[MatrixFigFullSaveName '.tif'])
        savefig(MatrixFigFullSaveName)
        pause(0.1)
        close
    end
    
    % -- Spatial distribution of decoding
    [X,Y]= meshgrid([-20:10:100],[-20:20:20]);
    Xlin = reshape(X',39,1); Ylin = reshape(Y',39,1);
    decodedMatrix = zeros(39,39);
    deltaDecoded2Actual_sh = NaN(1,780);
    
    counter = 1;
    if allPositions == 1
        for iii = 1:nRuns
            for j = 1:nTrialTypes
                for i = 1:nDecodedTrials_perTrialtype
                    [~, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
                    [~, I2] = max(PosteriorMatrix_perStimType_shuffled{iii}(:,j,i));
                    decoded_pos = mod(I,13);
                    if decoded_pos == 0; decoded_pos = 13; end
                    decoded_pos_shuffle = mod(I2,13);
                    if decoded_pos_shuffle == 0; decoded_pos_shuffle = 13; end
                    actual_pos = mod(j,13);
                    if actual_pos == 0; actual_pos = 13; end
                    delta_decod2actual_abs(counter) = abs(actual_pos-decoded_pos);
                    delta_decodShuff2actual_abs(counter) = abs(actual_pos-decoded_pos_shuffle);
                    delta_decod2actual(counter) = actual_pos-decoded_pos;
                    delta_decodShuff2actual(counter) = actual_pos-decoded_pos_shuffle;
                    decodedMatrix(I,j) = decodedMatrix(I,j)+1;
                    
                    delta_decod2actualEle_abs(counter) = abs(floor((I-1)/13)-floor((j-1)/13));
                    delta_decod2actualEleShuffle_abs(counter) = abs(floor((I2-1)/13)-floor((j-1)/13));
                    
                    % % - distance 2 actual stim
                    deltaDecoded2Actual(counter) = norm([Xlin(I),Ylin(I)]-[Xlin(j),Ylin(j)]);
                    deltaDecoded2Actual_sh(counter) = norm([Xlin(I2),Ylin(I2)]-[Xlin(j),Ylin(j)]);
                    
                    counter = counter+1;
                end
            end
        end
    end
    
    if savefigFlag
        HistFig = figure; hold on
        histogram(deltaDecoded2Actual, 'DisplayStyle', 'stairs','edgecolor','r','linewidth',5);
        histogram(deltaDecoded2Actual_sh, 'DisplayStyle', 'stairs','edgecolor','k','linewidth',5)
        xlabel('Relative distance to actual position (degrees)');
        ylabel('counts')
        legend({'data','shuffle'})
        title(['Decoded position vs. actual position' newline name ' - ' dataset],'interpreter','none')
        HistFigFullSaveName=[MainFolder name '_HistFig'];
        saveas(HistFig,[HistFigFullSaveName '.tif'])
        savefig(HistFigFullSaveName)
        pause(0.1)
        close
        
    end
    %     end % data concatenated per mouse -- spatial plots
    
    %% decoding the az position
    Decoded2Actual = NaN(nRuns,nTrialTypes,nDecodedTrials_perTrialtype);
    Decoded2Actual_sh = NaN(nRuns,nTrialTypes,nDecodedTrials_perTrialtype);
    for iii = 1:nRuns
        for j = 1:nTrialTypes
            for i = 1:nDecodedTrials_perTrialtype
                % % - Accuracy
                [~, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
                [~, I2] = max(PosteriorMatrix_perStimType_shuffled{iii}(:,j,i));
                decodedPos = mod(I,13);
                if decodedPos == 0; decodedPos = 13; end;
                decodedPos_sh = mod(I2,13);
                if decodedPos_sh == 0; decodedPos_sh = 13; end;
                actualPos = mod(j,13);
                if actualPos == 0; actualPos = 13; end;
                if decodedPos==actualPos || decodedPos==actualPos-1 || decodedPos==actualPos+1
                    DecodedAz(iii,j,i) = 1;
                else
                    DecodedAz(iii,j,i) = 0;
                end
                if decodedPos_sh==actualPos || decodedPos_sh==actualPos-1 || decodedPos_sh==actualPos+1
                    DecodedAz_shuffle(iii,j,i) = 1;
                else
                    DecodedAz_shuffle(iii,j,i) = 0;
                end
                % % - Distance 2 actual stim pos (degrees)
                Decoded2Actual(iii,j,i) = norm([Xlin(I),Ylin(I)]-[Xlin(j),Ylin(j)]);
                Decoded2Actual_sh(iii,j,i) = norm([Xlin(I2),Ylin(I2)]-[Xlin(j),Ylin(j)]);
                if isnan(Decoded2Actual_sh(iii,j,i))
                    keyboard
                end
                
            end
        end
    end
    
    
    % % - Accuracy
    b = squeeze(sum(DecodedAz,3));
    b2 = squeeze(sum(DecodedAz_shuffle,3));
    d = mean(b,1)./nDecodedTrials_perTrialtype;
    e = reshape(d,13,3);
    collapsed = mean(e,2);
    d2 = mean(b2,1)./nDecodedTrials_perTrialtype;
    e2 = reshape(d2,13,3);
    collapsed_shuffled = mean(e2,2);
    collapsed(:,2) = [-20:10:100];
    collapsed(:,2) = abs(collapsed(:,2)-V1az);
    [delta_az, II]= sort(collapsed(:,2));
    collapsed2(:,2) = collapsed(II,1);
    collapsed2(:,1) = delta_az;
    collapsed2_shuffle(:,2) = collapsed_shuffled(II,1);
    collapsed2_shuffle(:,1) = delta_az;
    
    % % - Distance 2 actual stim pos (degrees)
    b = mean(Decoded2Actual,3);
    b2 = mean(Decoded2Actual_sh,3);
    d = mean(b,1);
    e = reshape(d,13,3);
    dist2actual_az = mean(e,2);
    dist2actual_ele = mean(e,1);
    d2 = mean(b2,1);
    e2 = reshape(d2,13,3);
    dist2actual_az_sh = mean(e2,2);
    dist2actual_ele_sh = mean(e2,1);
    
    xvect = [-20:10:100];
    xvect = abs(xvect-V1az);
    [delta_az, II]= sort(xvect);
    
    dist2actual_az_V1RF(:,2) = dist2actual_az(II,1);
    dist2actual_az_V1RF(:,1) = delta_az;
    dist2actual_az_sh_V1RF(:,2) = dist2actual_az_sh(II,1);
    dist2actual_az_sh_V1RF(:,1) = delta_az;
    
    if ~isnan(V1az)
        % % - Accuracy
        
        [~,binnum] = sort(dist2actual_az_V1RF(:,1)/10);
        if binnum(1) == 0; binnum=binnum+1;end
        binsum = accumarray(binnum(:), dist2actual_az_V1RF(:,2),[],@mean);
        binsum_shuffle = accumarray(binnum(:), dist2actual_az_sh_V1RF(:,2),[],@mean);
        
        
        f = std(b,[],1)./(sqrt(size(b,1))*nDecodedTrials_perTrialtype);
        g = reshape(f,13,3);
        collapsed_sd = mean(g,2);
        collapsed_sd(:,2) = collapsed(:,2);
        collapsed2_sd(:,2) = collapsed_sd(II,1);
        collapsed2_sd(:,1) = delta_az;
        binnum_sd = ceil(collapsed2_sd(:,1)/10);
        if binnum_sd(1) == 0; binnum_sd=binnum_sd+1;end
        binsum_sd = accumarray(binnum_sd(:), collapsed2_sd(:,2),[],@mean);
    end
    
    if ~isempty(MainFolder)
        DecoderFig = figure('Name','DecoderFig');
        subplot(2,1,1); I = imagesc(avg_PerStimType);
        if allPositions == 1
            xlabel('Stimulus ID'); ylabel('Decoded Stimulus ID')
            I.XData = [1 39];I.YData=[1 39];axis tight; axis equal
        elseif allPositions == 2
            xlabel('Actual Pos.'); ylabel('Decoded Pos.')
            I.XData = [-20 100];I.YData=[-20 100];axis tight; axis equal
        elseif allPositions == 3
            xlabel('Actual Pos.'); ylabel('Decoded Pos.')
        end
        colormap gray
        c = colorbar; c.Label.String = 'log(likelihood)';
        c.Position = [0.7 0.585 0.01 0.342]; %[x y width height]
        avg_PerStimType2 = avg_PerStimType - mean(avg_PerStimType,1);
        avg_PerStimType = avg_PerStimType2 ./ std(avg_PerStimType,1);
        
        h=subplot(2,1,2); hold on
        plot([-20:10:100],collapsed(:,1),'o-','Color',[1 0 0]);
        plot([-20:10:100],collapsed_shuffled,'o-','Color',[0 0 0]);
        xlim([-22 102]); xlabel('Speaker Az position');ylabel('Decoding accuracy')
        ylim([0 1]);
        yl = ylim;yticks([0 .5 1])
        plot([V1az V1az],[yl(1) yl(2)],'k:')
        text(V1az,1,'\itV1az',...
            'horizontalalignment','left','verticalalignment','top','FontSize',6)
        plot([-22 102],[3/13 3/13],'k:');
        ylabel('Decoding Accuracy')
        xlabel('Stimulus Position')
        
        sgtitle([name '| n=' num2str(nROIs) ' axons'],'interpreter','none')
        
        if savefigFlag
            saveas(DecoderFig,[MainFolder filesep name '.tif'])
            pause(0.1)
            close
        end
    end
    
    %% plot single trials (Fig 3a)
    if plotMatrix
        X = -20:10:100; Y = 20:-20:-20;
        [x,y] = meshgrid(X,Y);
        x_reshaped = cat(2,x(1,:),x(2,:),x(3,:));
        y_reshaped = cat(2,y(1,:),y(2,:),y(3,:));
        trialTypesToPlot = [14 20 26];
        
        likelihoodDistrib_all = cat(3,PosteriorMatrix_perStimType{:}); % concatenate as many matrices as their are in the cell array "PosteriorMatrix_perStimType"
        
        for k = trialTypesToPlot
            figure
            for trial = 1:size(likelihoodDistrib_all,3)
                subplot(5,4,trial); hold on
                stimType = likelihoodDistrib_all(:,k,trial);
                r = reshape(stimType,13,3);
                
                I = imagesc(r');
                I.XData = [-20 100];
                I.YData = [20 -20];
                ax=gca;
                ax.YDir = 'normal';
                axis equal; axis tight
                colormap gray
                [~,idx] = max(stimType);
                if idx<14; y_decoded = 20; elseif idx <27; y_decoded = 0; idx = idx-13;else  y_decoded = -20;        idx = idx-26;
                end
                plot([X(idx)-5 X(idx)+5 X(idx)+5],...
                    [y_decoded-10 y_decoded-10 y_decoded+10 ],...
                    'color',[1 0 1],'linewidth',1)
                plot([x_reshaped(k)+5 x_reshaped(k)-5 x_reshaped(k)-5],...
                    [y_reshaped(k)+10 y_reshaped(k)+10 y_reshaped(k)-10],...
                    'color',[0 1 1],'linewidth',1);
                
                axis off
            end % loop through stim type
            zetitle = [name '_' num2str(k)];
            sgtitle(zetitle,'interpreter','none')
            savefig(gcf,['C:\Users\camil\AudSpace\Illustration\Decoder\test2',filesep, zetitle '.fig'])
            saveas(gcf,['C:\Users\camil\AudSpace\Illustration\Decoder\test2',filesep, zetitle '.tif'])
            pause(0.1)
            close
        end % loop through trials
    end
    
    %% save
    % % - Accuracy
    Posterior.collapsed          = dist2actual_az;
    Posterior.collapsed_shuffled = dist2actual_az_sh;
    Posterior.Elecollapsed          = dist2actual_ele;
    Posterior.Elecollapsed_shuffled = dist2actual_ele_sh;
    
    Posterior.delta_decod2actual      = delta_decod2actual;
    Posterior.delta_decodShuff2actual = delta_decodShuff2actual;
    Posterior.delta_decod2actual_abs      = delta_decod2actual_abs;
    Posterior.delta_decodShuff2actual_abs = delta_decodShuff2actual_abs;
    
    Posterior.delta_decod2actualEle_abs = delta_decod2actualEle_abs;
    Posterior.delta_decod2actualEleShuffle_abs = delta_decod2actualEleShuffle_abs;
    Posterior.deltaDecoded2Actual = deltaDecoded2Actual;
    Posterior.deltaDecoded2Actual_sh = deltaDecoded2Actual_sh;
    
    Posterior.avg_PerStimType = avg_PerStimType;
    if ~strcmp(name,'MickeyMouse')
        Posterior.binsum         = binsum;
        Posterior.binsum_shuffle = binsum_shuffle;
    end
    % save([saveDir_variable, '\Posterior.mat'],'Posterior')
else
    disp([name ' skipped'])
    Posterior=[];
end
end



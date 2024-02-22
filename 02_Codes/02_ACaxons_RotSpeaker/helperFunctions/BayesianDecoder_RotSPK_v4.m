% July2022: adapted from _v4 for rotating speaker experiment

function Posterior = BayesianDecoder_RotSPK_v4(data,training_fraction,V1az,V1ele,axonRois,info,name,MainFolder)
% Plot the posterior distribution for all stimuli locations
% training_fraction: fraction of trials used for training

%% Parametrize
% allPositions = 1; % train the decoder using both az and ele (=1), az only (=2) or "chuncks" (=3)
plotMatrix =  0;  % the matrix of all likelihood for all stim positions
savefigFlag = 1;
dataType = info.metric; % 'dFF' vs 'spikes' or 'Fval'
base_window = info.base_window;
resp_window = info.resp_window;

%% Define some variables
nAzPos = 12;
nTrialTypes = 12;
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
for i = 1:length(axonRois)
    axonRois{1,i}(axonRois{1,i}==0)=[];
    if length(axonRois{1,i})>1
        %     for ii = 1:length(axonRois{1,i})
        meandFF = mean(mean(mean(data(axonRois{1,i},:,:,:),2),3),4);
        %     end
        [~,idx_maxdFF] = max(meandFF);
        RepROI(i) = axonRois{1,i}(idx_maxdFF);
    else
        RepROI(i) = axonRois{1,i};
    end
end
selected_RepRois = false(1,nROIs);
selected_RepRois(RepROI)=true;
nROIs = sum(selected_RepRois);

if sum(selected_RepRois)>10
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
        switch dataType
            case 'spikes'
                spikesDist = squeeze(mean(data_selected.ONresp.spikesDist(selected_RepRois,24:30,:,:),2));
                sdMatrix = squeeze(std(spikesDist(selected_RepRois,:,trainingTrials(iii,:)),[],3,'omitnan'));
                RespMatrix = squeeze(nanmean(spikesDist(selected_RepRois,:,trainingTrials(iii,:)),3));
                meanResp = spikesDist(selected_RepRois,:,decodedTrials(iii,:));
                
            case 'dFF'
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
        end
        meanResp = reshape(meanResp,nROIs,nDecodedTrials);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% SHUFFLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i = 1:nTrialTypes
            for j = 1:nTrainingTrials
                shuffle_trainingTrials(:,i,j) = respDist(randperm(nROIs),i,j);
            end
        end

        RespMatrix_shuffled = squeeze(nanmean(shuffle_trainingTrials,3));
        sdMatrix_shuffled = squeeze(std(shuffle_trainingTrials,[],3));
        for i = 1:nDecodedTrials
            shuffle_meanResp(:,i) = meanResp(randperm(size(meanResp,1)),i);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    %% Normalization after each trial. I don't think it is the way to go
    % for i = 1:nRuns
    %     a{i}(:,:) = exp(PosteriorMatrix{i}).*1/39;
    % %     for j = 1:nDecodedTrials
    % %         b{i}(:,j) = a{i}(:,j)/max(a{i}(:,j));
    % %     end
    % end
    
    %% For shuffle control, after each trial
    % for i = 1:nBatches
    % mean_PosteriorMatrix_shuffled{i} = mean(PosteriorMatrix_shuffled{i},1);
    % mean_PosteriorMatrix_shuffled{i} = repmat(mean_PosteriorMatrix_shuffled{i},39,1);
    % b{i} = mean_PosteriorMatrix_shuffled{i}./PosteriorMatrix{i};
    % end
    
    %% Stack into trials of same type
    for iii = 1:nRuns
        PosteriorMatrix_perStimType{iii} = reshape(PosteriorMatrix{iii},...
            nTrialTypes,nTrialTypes,nDecodedTrials_perTrialtype);
        PosteriorMatrix_perStimType_shuffled{iii} = reshape(PosteriorMatrix_shuffled{iii},...
            nTrialTypes,nTrialTypes,nDecodedTrials_perTrialtype);
    end

    c = cat(3,PosteriorMatrix_perStimType{:}); % concatenate as many matrices as their are in the cell array "PosteriorMatrix_perStimType"
    avg_PerStimType = squeeze(nanmean(c,3));
    d = cat(3,PosteriorMatrix_perStimType_shuffled{:}); % concatenate as many matrices as their are in the cell array "PosteriorMatrix_perStimType"
    avg_PerStimType_shuffle = squeeze(nanmean(d,3));
    
    %% decoding the az position
    DecoderFig = figure('Name','DecoderFig');
    if allPositions == 1
        for iii = 1:nRuns
            for j = 1:nTrialTypes
                for i = 1:nDecodedTrials_perTrialtype
                    [~, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
                    [~, I2] = max(PosteriorMatrix_perStimType_shuffled{iii}(:,j,i));
                    if  I==j || I==j+13 || I==j+26 || I==j-13 || I==j-26 || I==j+1 || I==j-1 || I==j+14 || I==j+12 || I==j+27 || I==j+25 || I==j-12 || I==j-14 || I==j-25 || I==j-27
                        DecodedAz(iii,j,i) = 1;
                    else
                        DecodedAz(iii,j,i) = 0;
                    end
                    if I2==j || I2==j+13 || I2==j+26 || I2==j-13 || I2==j-26 || I2==j+1 || I2==j-1 || I2==j+14 || I2==j+12 || I2==j+27 || I2==j+25 || I2==j-12 || I2==j-14 || I2==j-25 || I2==j-27
                        %                 if  I2 == j || I2 == j+13 || I2 == j+26 || I2 == j-13 || I2 == j-26
                        DecodedAz_shuffle(iii,j,i) = 1;
                    else
                        DecodedAz_shuffle(iii,j,i) = 0;
                    end
                end
            end
        end
        
        b = squeeze(sum(DecodedAz,3));
        b2 = squeeze(sum(DecodedAz_shuffle,3));
        % d = sum(b,1);
        % d = d/20;
        % d = d*13;
        d = mean(b,1)./nDecodedTrials_perTrialtype;
        e = reshape(d,13,3);
        collapsed = mean(e,2);
        d2 = mean(b2,1)./nDecodedTrials_perTrialtype;
        e2 = reshape(d2,13,3);
        collapsed_shuffled = mean(e2,2);
        % figure;imagesc(e');
        collapsed(:,2) = [-20:10:100];
        collapsed(:,2) = abs(collapsed(:,2)-V1az);
        [delta_az, II]= sort(collapsed(:,2));
        collapsed2(:,2) = collapsed(II,1);
        collapsed2(:,1) = delta_az;
        collapsed2_shuffle(:,2) = collapsed_shuffled(II,1);
        collapsed2_shuffle(:,1) = delta_az;
        % figure;plot(delta_az,collapsed2(:,2))
        
        binnum = ceil(collapsed2(:,1)/10);
        if binnum(1) == 0; binnum=binnum+1;end
        binsum = accumarray(binnum(:), collapsed2(:,2),[],@mean);
        binsum_shuffle = accumarray(binnum(:), collapsed2_shuffle(:,2),[],@mean);
        
        h=subplot(2,3,5);
        set(h,'position',[0.5 0.1100 0.3 0.3412]); hold on
        
        p1=plot(((1:size(binsum,1))-1)*10, binsum,'o-','Color',[1 0 0]);
        p2=plot(((1:size(binsum_shuffle,1))-1)*10, binsum_shuffle,'o-','Color',[0 0 0]);
        
        f = std(b,[],1)./(sqrt(size(b,1))*nDecodedTrials_perTrialtype);
        g = reshape(f,13,3);
        collapsed_sd = mean(g,2);
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
        
        h=subplot(2,3,4); hold on
        set(h,'position',[0.13 0.1100 0.3 0.3412]); hold on
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
        %     exp = shiftdim(DecodedAz,2);
        %     shuffle = shiftdim(DecodedAz_shuffle,2);
        %     exp = reshape(exp,39*4*5,1);
        %     shuffle = reshape(shuffle,39*4*5,1);
        %     anova_matrix = [exp shuffle];
        %     [p,tbl,stats] = anova2(anova_matrix,20);
        %     multcompare(stats)
        %
        %     exp = shiftdim(DecodedAz,2);
        %     shuffle = shiftdim(DecodedAz_shuffle,2);
        %     exp = reshape(exp,4*5,39);
        %     exp = reshape(exp,20,13,3);
        %     exp = mean(exp,3);
        %     shuffle = reshape(shuffle,4*5,39);
        %     shuffle = reshape(shuffle,20,13,3);
        %     shuffle = mean(shuffle,3);
        %     groups=[ones(20,1);2*ones(20,1)];
        %     anova_matrix = [exp;shuffle];
        % t=table(groups,anova_matrix(:,1),anova_matrix(:,2),anova_matrix(:,3),anova_matrix(:,4),anova_matrix(:,5),anova_matrix(:,6),anova_matrix(:,7),anova_matrix(:,8),anova_matrix(:,9),anova_matrix(:,10),anova_matrix(:,11),anova_matrix(:,12),anova_matrix(:,13),...
        %    anova_matrix(:,14),anova_matrix(:,15),anova_matrix(:,16),anova_matrix(:,17),anova_matrix(:,18),anova_matrix(:,19),anova_matrix(:,20),anova_matrix(:,21),anova_matrix(:,22),anova_matrix(:,23),anova_matrix(:,24),anova_matrix(:,25),anova_matrix(:,26),...
        %    anova_matrix(:,27),anova_matrix(:,28),anova_matrix(:,29),anova_matrix(:,30),anova_matrix(:,31),anova_matrix(:,32),anova_matrix(:,33),anova_matrix(:,34),anova_matrix(:,35),anova_matrix(:,36),anova_matrix(:,37),anova_matrix(:,38),anova_matrix(:,39),...
        %    'VariableNames',{'groups','meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10','meas11','meas12','meas13','meas14','meas15','meas16','meas17','meas18','meas19','meas20','meas21','meas22','meas23','meas24','meas25','meas26','meas27','meas28','meas29','meas30','meas31','meas32','meas33','meas34','meas35','meas36','meas37','meas38','meas39'});
        % Az = table([1:39]','VariableNames',{'Measurements'});
        %     t=table(groups,anova_matrix(:,1),anova_matrix(:,2),anova_matrix(:,3),anova_matrix(:,4),anova_matrix(:,5),anova_matrix(:,6),anova_matrix(:,7),anova_matrix(:,8),anova_matrix(:,9),anova_matrix(:,10),anova_matrix(:,11),anova_matrix(:,12),anova_matrix(:,13),...
        %         'VariableNames',{'groups','meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10','meas11','meas12','meas13'});
        %     Az = table([1:13]','VariableNames',{'Measurements'});
        %     rm = fitrm(t,'meas1-meas13~groups+groups','WithinDesign',Az,'WithinModel','separatemeans');
        %     ranovatbl = ranova(rm);
        
    elseif allPositions == 2
        for iii = 1:nRuns
            for j = 1:nAzPos
                for i = 1:nDecodedTrials_perTrialtype*nElePos
                    [~, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
                    [~, I2] = max(PosteriorMatrix_perStimType_shuffled{iii}(:,j,i));
                    if  I == j
                        DecodedAz(iii,j,i) = 1;
                    else
                        DecodedAz(iii,j,i) = 0;
                    end
                    if  I2 == j
                        DecodedAz_shuffle(iii,j,i) = 1;
                    else
                        DecodedAz_shuffle(iii,j,i) = 0;
                    end
                end
            end
        end
        
        
        b = squeeze(sum(DecodedAz,3));
        b2 = squeeze(sum(DecodedAz_shuffle,3));
        
        collapsed = mean(b,1)./(nDecodedTrials_perTrialtype*3);
        collapsed_shuffled = mean(b2,1)./(nDecodedTrials_perTrialtype*3);
        % figure;imagesc(e');
        collapsed(2,:) = [-20:10:100]';
        collapsed = collapsed';
        collapsed(:,2) = abs(collapsed(:,2)-V1az);
        [delta_az, II]= sort(collapsed(:,2));
        collapsed2(:,2) = collapsed(II,1);
        collapsed2(:,1) = delta_az;
        collapsed_shuffled = collapsed_shuffled';
        collapsed2_shuffle(:,2) = collapsed_shuffled(II,1);
        collapsed2_shuffle(:,1) = delta_az;
        % figure;plot(delta_az,collapsed2(:,2))
        
        binnum = ceil(collapsed2(:,1)/10);%+1;
        if binnum(1) == 0;binnum=binnum+1;end
        binsum = accumarray(binnum(:), collapsed2(:,2),[],@mean);
        binsum_shuffle = accumarray(binnum(:), collapsed2_shuffle(:,2),[],@mean);
        
        
        h=subplot(2,3,5);
        set(h,'position',[0.5 0.1100 0.3 0.3412]); hold on
        p1=plot(((1:size(binsum,1))-1)*10, binsum,'o-','Color',[1 0 0]);
        p2=plot(((1:size(binsum_shuffle,1))-1)*10, binsum_shuffle,'o-','Color',[0 0 0]);
        xlabel('Degrees away from V1 azimuth')
        
        
        collapsed_sd = std(b,[],1)./(sqrt(size(b,1))*nDecodedTrials_perTrialtype*3);
        collapsed_sd = collapsed_sd';
        collapsed_sd(:,2) = collapsed(:,2);
        collapsed2_sd(:,2) = collapsed_sd(II,1);
        collapsed2_sd(:,1) = delta_az;
        binnum_sd = ceil(collapsed2_sd(:,1)/10);%+1;
        if binnum_sd(1) == 0;binnum_sd=binnum_sd+1;end
        binsum_sd = accumarray(binnum_sd(:), collapsed2_sd(:,2),[],@mean);
        
        p3=plot(((1:size(binsum,1))-1)*10, binsum+binsum_sd,'--','Color',[0.75 0.75 0.75]);
        plot(((1:size(binsum,1))-1)*10, binsum-binsum_sd,'--','Color',[0.75 0.75 0.75]);
        xlim([-2 122]);
        xl = xlim;
        % text(102,1/13,'\itchance',...
        %     'horizontalalignment','right','verticalalignment','bottom','FontSize',6)
        ylim([0 1]);yticks([])
        p4 = plot(xl,[1/13 1/13],'k:');
        lgd=legend([p1,p2,p3,p4],{'Data','Shuffle','sem ac. runs','chance'});
        lgd.Position = [0.8 0.2 0.19 0.14];
        lgd.FontSize = 7;
        %     exp = shiftdim(DecodedAz,2);
        %     shuffle = shiftdim(DecodedAz_shuffle,2);
        %     exp = reshape(exp,39*4*5,1);
        %     shuffle = reshape(shuffle,39*4*5,1);
        %     anova_matrix = [exp shuffle];
        %     [p,tbl,stats] = anova2(anova_matrix,60);
        %     multcompare(stats);
        
        % exp = shiftdim(DecodedAz,2);
        % shuffle = shiftdim(DecodedAz_shuffle,2);
        % exp = reshape(exp,size(exp,1)*nRuns,nAzPos);
        % shuffle = reshape(shuffle,size(shuffle,1)*nRuns,nAzPos);
        % groups=[ones(60,1);2*ones(60,1)];
        % anova_matrix = [exp;shuffle];
        % % t=table(groups,anova_matrix(:,1),anova_matrix(:,2),anova_matrix(:,3),anova_matrix(:,4),anova_matrix(:,5),anova_matrix(:,6),anova_matrix(:,7),anova_matrix(:,8),anova_matrix(:,9),anova_matrix(:,10),anova_matrix(:,11),anova_matrix(:,12),anova_matrix(:,13),...
        % %    anova_matrix(:,14),anova_matrix(:,15),anova_matrix(:,16),anova_matrix(:,17),anova_matrix(:,18),anova_matrix(:,19),anova_matrix(:,20),anova_matrix(:,21),anova_matrix(:,22),anova_matrix(:,23),anova_matrix(:,24),anova_matrix(:,25),anova_matrix(:,26),...
        % %    anova_matrix(:,27),anova_matrix(:,28),anova_matrix(:,29),anova_matrix(:,30),anova_matrix(:,31),anova_matrix(:,32),anova_matrix(:,33),anova_matrix(:,34),anova_matrix(:,35),anova_matrix(:,36),anova_matrix(:,37),anova_matrix(:,38),anova_matrix(:,39),...
        % %    'VariableNames',{'groups','meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10','meas11','meas12','meas13','meas14','meas15','meas16','meas17','meas18','meas19','meas20','meas21','meas22','meas23','meas24','meas25','meas26','meas27','meas28','meas29','meas30','meas31','meas32','meas33','meas34','meas35','meas36','meas37','meas38','meas39'});
        % % Az = table([1:39]','VariableNames',{'Measurements'});
        % t=table(groups,anova_matrix(:,1),anova_matrix(:,2),anova_matrix(:,3),anova_matrix(:,4),anova_matrix(:,5),anova_matrix(:,6),anova_matrix(:,7),anova_matrix(:,8),anova_matrix(:,9),anova_matrix(:,10),anova_matrix(:,11),anova_matrix(:,12),anova_matrix(:,13),...
        %     'VariableNames',{'groups','meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10','meas11','meas12','meas13'});
        % Az = table([1:13]','VariableNames',{'Measurements'});
        % rm = fitrm(t,'meas1-meas13~groups+groups','WithinDesign',Az,'WithinModel','separatemeans');
        % ranovatbl = ranova(rm);
        
        h=subplot(2,3,4); hold on
        set(h,'position',[0.13 0.1100 0.3 0.3412]); hold on
        plot([-20:10:100],collapsed(:,1),'o-','Color',[1 0 0]);
        plot([-20:10:100],collapsed_shuffled,'o-','Color',[0 0 0]);
        xlim([-22 102]); xlabel('Speaker Az position');ylabel('Decoding accuracy (+-10deg)')
        ylim([0 1]);
        yl = ylim;yticks([0 .5 1])
        plot([V1az V1az],[yl(1) yl(2)],'k:')
        text(V1az,1,'\itV1az',...
            'horizontalalignment','left','verticalalignment','top','FontSize',6)
        plot([-22 102],[1/13 1/13],'k:');
        
        
    elseif allPositions == 3
        for iii = 1:nRuns
            for j = 1:nPos
                for i = 1:nDecodedTrials_perTrialtype*13
                    [~, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
                    [~, I2] = max(PosteriorMatrix_perStimType_shuffled{iii}(:,j,i));
                    if  I == j
                        DecodedAz(iii,j,i) = 1;
                    else
                        DecodedAz(iii,j,i) = 0;
                    end
                    if  I2 == j
                        DecodedAz_shuffle(iii,j,i) = 1;
                    else
                        DecodedAz_shuffle(iii,j,i) = 0;
                    end
                end
            end
        end
        b = squeeze(sum(DecodedAz,3));
        b2 = squeeze(sum(DecodedAz_shuffle,3));
        
        collapsed = mean(b,1)./(nDecodedTrials_perTrialtype*13);
        collapsed_shuffled = mean(b2,1)./(nDecodedTrials_perTrialtype*13);
        collapsed(2,:) = [-20 40 100]';
        collapsed = collapsed';
        collapsed(:,2) = abs(collapsed(:,2)-V1az);
        [delta_az, II]= sort(collapsed(:,2));
        collapsed2(:,2) = collapsed(II,1);
        collapsed2(:,1) = delta_az;
        collapsed_shuffled = collapsed_shuffled';
        collapsed2_shuffle(:,2) = collapsed_shuffled(II,1);
        collapsed2_shuffle(:,1) = delta_az;
        
        binnum = ceil(collapsed2(:,1)/10);
        if binnum(1) == 0;binnum=binnum+1;end
        binsum = accumarray(binnum(:), collapsed2(:,2),[],@mean);
        binsum_shuffle = accumarray(binnum(:), collapsed2_shuffle(:,2),[],@mean);
        
        
        
        h2 = subplot(2,3,5); hold on
        set(h2,'position',[0.5 0.1100 0.3 0.3412]); hold on
        p1=bar([ collapsed2(:,2) collapsed2_shuffle(:,2)]);
        xlabel('Distance from V1RFaz')
        xticks([1:3]); xticklabels({'same','closest','farthest'})
        ylim([0  1])
        legend({'Data','Shuffle'})
        
        h1= subplot(2,3,4); hold on
        set(h1,'position',[0.13 0.1100 0.3 0.3412]); hold on
        p1=bar([ collapsed(:,1) collapsed_shuffled]);
        xlabel('Speaker Pos');
        xticks([1:3]); xticklabels({'[-20:10]','[20:60]','[70:100]'})
        ylim([0  1])
    end
    % % draft (from 2019)
    % for iii = 1:nBatches
    %     for j = 1:39 % for each trial type
    %         for i = 1:4
    %             [temp, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
    %             IndexDecoded(j,i,iii) = I;
    %         end
    %     end
    % end
    % IndexDecoded = reshape(IndexDecoded,39,4*5);
    % accuracy = zeros(1,39);
    % for k = 1:39
    %     diff = IndexDecoded(k,:)- k;
    %     for i =1:20
    %         if diff(i) == 0 || diff(i) == -13 || diff(i) == +13 ||diff(i) == 1 || diff(i) == -14 || diff(i) == +14 || diff(i) == -26 || diff(i) == +26 || diff(i) == -27 || diff(i) == +27
    %             accuracy(1,k) = accuracy(1,k)+1;
    %         end
    %     end
    % end
    %
    %         accuracy(1,1) = accuracy(1,1)*9/6;
    %         accuracy(1,13) = accuracy(1,13)*9/6;
    %         accuracy(1,14) = accuracy(1,14)*9/6;
    %         accuracy(1,26) = accuracy(1,26)*9/6;
    %         accuracy(1,27) = accuracy(1,27)*9/6;
    %         accuracy(1,39) = accuracy(1,39)*9/6;
    %
    % e = reshape(accuracy,13,3);
    % figure;imagesc(e');
    %  collapsed = mean(e,2);
    % % figure;imagesc(e');
    % collapsed(:,2) = [-20:10:100];
    % collapsed(:,2) = abs(collapsed(:,2)-V1az);
    % [delta_az, II]= sort(collapsed(:,2));
    % collapsed2(:,2) = collapsed(II,1);
    % collapsed2(:,1) = delta_az;
    % % figure;plot(delta_az,collapsed2(:,2))
    %
    % binnum = ceil(collapsed2(:,1)/10)+1;
    % binsum = accumarray(binnum(:), collapsed2(:,2),[],@mean);
    % figure;
    % hold on
    % plot(((1:size(binsum,1))-1)*10, binsum,'o-','Color',[0.75 0.75 0.75]);
    % xlabel('Distance From RF center')
    % ylabel('Normalized accuracy')
    % saveas(gcf,[saveDir_Fig,'\DistanceFromRFcenter.tif'])
    % saveas(gcf,[saveDir_Fig,'\DistanceFromRFcenter.fig'])
    
    %% how far from actual position?
    % AccuratePos = zeros(nRuns,39,i); TenAwayPos = zeros(nRuns,39,i); TwentyAwayPos = zeros(nRuns,39,i); Inaccurate = zeros(nRuns,39,i);
    % for iii = 1:nRuns
    %     for j = 1:39 % for each trial type
    %         for i = 1:4
    %             [~, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
    %             if  I == j
    %                 AccuratePos(iii,j,i) = 1;
    %             elseif I == j+1 || I == j-1
    %                 TenAwayPos(iii,j,i) = 1;
    %             elseif I == j+2 || I == j-2
    %                 TwentyAwayPos(iii,j,i) = 1;
    %             else
    %                 Inaccurate (iii,j,i) = 1;
    %             end
    %         end
    %     end
    % end
    %  b1 = squeeze(sum(AccuratePos,3));
    % b2 = squeeze(sum(TenAwayPos,3));
    %  b3 = squeeze(sum(TwentyAwayPos,3));
    %  b = [b1;b2;b3];
    %  d = sum(b,1);
    %  e =  reshape(d,13,3);
    % %  figure;imagesc(e');
    %   collapsed = mean(e,2);
    %   collapsed(:,2) = [-20:10:100];
    % collapsed(:,2) = abs(collapsed(:,2)-V1az);
    % [delta_az, II]=sort(collapsed(:,2));
    % collapsed2(:,2) = collapsed(II,1);
    % collapsed2(:,1) = delta_az;
    % % figure;plot(delta_az,collapsed2(:,2))
    %
    % binnum = ceil(collapsed2(:,1)/10);%+1;
    % binsum = accumarray(binnum(:), collapsed2(:,2),[],@mean);
    % figure;
    % hold on
    % plot(((1:size(binsum,1))-1)*10, binsum,'o-','Color',[0.75 0.75 0.75]);
    %
    % % binnum = ceil(Posterior2.draft.collapsed2(:,1)/10);%+1;
    % % binsum = accumarray(binnum(:), Posterior2.draft.collapsed2(:,2),[],@mean);
    % % plot(((1:size(binsum,1))-1)*10, binsum,'o-','Color',[0.5 0.5 0.5]);
    %
    % X = [Posterior.draft.collapsed2(:,1); Posterior2.draft.collapsed2(:,1)];
    % Y = [Posterior.draft.collapsed2(:,2); Posterior2.draft.collapsed2(:,2)];
    % binnum = ceil(X/10)+1;
    % binsum = accumarray(binnum(:), Y(:),[],@mean);
    % plot(((1:size(binsum,1))-1)*10, binsum,'k-');
    %
    % plot([0 90], [1 1],'--','Color',[0.7 0.7 0.7]);
    % legend({'CMad37pos5','CMad34pos4','mean','chance'},'Location','northeast')
    % xlabel('distance to V1 azimuth positon')
    % xlim([0 70])
    % ylabel('accuracy above chance')
    
    
    
    % %% Towards plotting single session data
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%% Relative quantification %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% normalize to chance?
    % avg_PerStimType = abs(avg_PerStimType_shuffle)./abs(avg_PerStimType);
    %% normalize after, at the average level to not lose info on the quality of each trial?
    % avg_PerStimType2 = exp(avg_PerStimType)*1/39;
    % avg_PerStimType2 = avg_PerStimType-min(avg_PerStimType,[],1);
    % avg_PerStimType2 = avg_PerStimType2./max(avg_PerStimType2,[],1);
    
    %% "normalize" the distribution: distribution with mean of 0 and sd of 1
    subplot(2,2,1); I = imagesc(avg_PerStimType);
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
    c.Position = [0.44 0.585 0.01 0.342]; %[x y width height]
    avg_PerStimType2 = avg_PerStimType - mean(avg_PerStimType,1);
    avg_PerStimType = avg_PerStimType2 ./ std(avg_PerStimType,1);
    
    % avg_PerStimType2 = exp(avg_PerStimType);
    % avg_PerStimType2 = avg_PerStimType2./max(avg_PerStimType2,[],1);
    subplot(2,2,2); I = imagesc(avg_PerStimType);
    if allPositions == 1
        xlabel('Stimulus ID')
        I.XData = [1 39];I.YData=[1 39];
    elseif allPositions == 2
        xlabel('Actual Pos.')
        I.XData = [-20 100];I.YData=[-20 100];
    elseif allPositions == 3
        xlabel('Actual Pos.')
    end
    axis tight; axis equal
    colormap gray
    c = colorbar; c.Label.String = '~log(norm likelihood)';
    c.Position = [0.87 0.585 0.01 0.342]; %[x y width height]
    try
        sgtitle([name '| n=' num2str(nROIs) ' axons'],'interpreter','none')
    catch
        sp = suptitle( [data_selected.name '|V1az:' num2str(V1az)]);
        sp.Interpreter = 'none';
    end
    if savefigFlag
        % %     saveDir = [MainFolder filesep 'Plots' filesep 'BayesianDecoder'];
        %     if ~exist(MainFolder,'dir')
        %         mkdir(MainFolder)
        %     end
        saveas(DecoderFig,[MainFolder filesep name '.tif'])
        pause(0.1)
        close
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% plot
    if plotMatrix
        if allPositions == 1
            % -- Estimated likelihood averaged across trial, for each position
            figure
            X = -20:10:100; Y = 20:-20:-20;
            [x,y] = meshgrid(X,Y);
            x_reshaped = cat(2,x(1,:),x(2,:),x(3,:));
            y_reshaped = cat(2,y(1,:),y(2,:),y(3,:));
            ymax = max(max(avg_PerStimType));
            ymin = min(min(avg_PerStimType));
            for k =1:nTrialTypes
                subplot(3,13,k); hold on
                stimType = avg_PerStimType(:,k);
                r = reshape(stimType,13,3);
                %     I = imagesc(r',[0 0.6]);
                I = imagesc(r',[ymin ymax]);
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
                
                %             [0 0 0],'edgecolor',[0 1 0],'facecolor','none','edgealpha',1);
                %             [y_reshaped(k)-10 y_reshaped(k)-10 y_reshaped(k)+10 y_reshaped(k)+10],...
                %             [0 0 0],'edgecolor',[0 1 0],'facecolor','none','edgealpha',1);
                %         fill([X(idx)-5 X(idx)+5 X(idx)+5 X(idx)-5],...
                %                     [y_decoded-10 y_decoded-10 y_decoded+10 y_decoded+10],...
                %                     [0 0 0],'edgecolor',[1 0 1],'facecolor','none','edgealpha',1,'linewidth',1)
                %         fill([x_reshaped(k)-5 x_reshaped(k)+5 x_reshaped(k)+5 x_reshaped(k)-5],...
                %             [y_reshaped(k)-10 y_reshaped(k)-10 y_reshaped(k)+10 y_reshaped(k)+10],...
                %             [0 0 0],'edgecolor',[0 1 0],'facecolor','none','edgealpha',1);
                %             if ~strcmp(data_selected.name, 'AllData'); plot(V1az,V1ele,'wo')'; end
                if k == 27
                    xlabel('Azimuth')
                    ylabel('Elevation')
                    ax = gca;
                    % ax.YTick = [-20 0 20];
                    % ax.XTickLabel = [X(1):20:X(end)];
                    ax.YTick = [-20 0 20];
                    ax.YTickLabel = -Y;
                elseif k == 39
                    c = colorbar;
                    c.Label.String = 'Probability';
                    c.Position = [0.96 0.13 0.005 0.135]; %[x y width height]
                    axis off
                else
                    axis off
                end
            end
            sgtitle('Position estimation, for each stim location')
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.9, 0.35]); %[x y width height]
            % saveas(gcf,[saveDir_Fig,'\PositionEstimate_all.tif'])
            % saveas(gcf,[saveDir_Fig,'\PositionEstimate_all.fig'])
            
            % -- Estimated likelihood at correct position averaged across trial, for each position
            for k = 1:nTrialTypes
                DecodingAtActualPos(k) = avg_PerStimType(k,k); % it's already normalized
            end
            r = reshape(DecodingAtActualPos,13,3);
            figure; I = imagesc(r'); axis equal; axis tight
            I.XData = [-20 100];
            I.YData = [-20 20];
            ax = gca;
            % ax.YTick = [-20 0 20];
            % ax.XTickLabel = [X(1):20:X(end)];
            ax.YTick = [-20 0 20];
            ax.YTickLabel = Y;
            colormap gray
            hold on; plot(V1az,V1ele,'ko')
            c = colorbar;
            c.Label.String = 'Probability';
            xlabel('Azimuth')
            ylabel('Elevation')
            title('Decoded probability at the correct location, for each location')
            % saveas(gcf,[saveDir_Fig,'\ProbabilitiesAtCorrectLocation_all.tif'])
            % saveas(gcf,[saveDir_Fig,'\ProbabilitiesAtCorrectLocation_all.fig'])
            
        elseif allPositions == 2
            % -- Estimated likelihood averaged across trial, for each position
            X = -20:10:100;
            figure; hold on
            ymax = max(max(avg_PerStimType));
            ymin = min(min(avg_PerStimType));
            for k =1:nAzPos
                subtightplot(1,13,k)
                stimType = avg_PerStimType(:,k);
                I = imagesc(stimType',[ymin ymax]);
                I.XData = [-20 100];
                axis tight
                colormap gray
                hold on; plot(X(k),1,'b+'); plot(V1az,1,'ro')
                if k == 1
                    xlabel('Azimuth')
                    yticks([])
                elseif k == 13
                    c = colorbar;
                    c.Label.String = 'Probability';
                    c.Position = [0.96 0.13 0.005 0.135]; %[x y width height]
                    axis off
                else
                    axis off
                end
            end
            sgtitle('Position estimation, for each stim location')
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.8, 0.4]); %[x y width height]
            % saveas(gcf,[saveDir_Fig,'\PositionEstimate_all.tif'])
            % saveas(gcf,[saveDir_Fig,'\PositionEstimate_all.fig'])
            
            % -- Estimated likelihood at correct position averaged across trial, for each position
            for k = 1:nAzPos
                DecodingAtActualPos(k) = avg_PerStimType(k,k); % it's already normalized
            end
            figure; I = imagesc(DecodingAtActualPos);
            I.XData = [-20 100];
            colormap gray
            hold on; plot(V1az,V1ele,'ko')
            c = colorbar;
            c.Label.String = 'Probability';
            xlabel('Azimuth')
            ylabel('Elevation')
            title('Decoded probability at the correct location, for each location')
            % saveas(gcf,[saveDir_Fig,'\ProbabilitiesAtCorrectLocation_all.tif'])
        end
    end
    
    %% save
    % Posterior.Matrix = PosteriorMatrix;
    % Posterior.Matrix_shuffled = PosteriorMatrix_shuffled;
    % Posterior.Matrix_perStimType =PosteriorMatrix_perStimType;
    % Posterior.Matrix_perStimType_shuffled = PosteriorMatrix_perStimType_shuffled;
    % Posterior.avg_PerStimType = avg_PerStimType;
    % Posterior.avg_PerStimType_shuffle = avg_PerStimType_shuffle;
    Posterior.binsum = binsum;
    Posterior.binsum_shuffle = binsum_shuffle;
    Posterior.collapsed = collapsed;
    Posterior.collapsed_shuffled = collapsed_shuffled;
    % save([saveDir_variable, '\Posterior.mat'],'Posterior')
else
    disp([name ' skipped'])
    Posterior=[];
end
end



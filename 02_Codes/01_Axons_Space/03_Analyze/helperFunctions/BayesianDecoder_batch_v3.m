function Decoder_Output = BayesianDecoder_batch_v3(RespROIs,dataset,MainFolder,doSingleSession,TwentykHz,plotSingleSession,doMouseagregate,doAccVsV1)
%% a few parameters
corrThresh = 0.3;         % to determine which ROIs belong to the same axon
training_fraction = 0.8;

yl = [15 70];            % ylim to use for the decoding distance 2 actual plots

%% a few definitions
n_animals = size(RespROIs.info.animalID,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
saveFolder = MainFolder;
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

%% do Bayesian decoder on single session (and plot and save)
if doSingleSession
    if plotSingleSession == 0
        saveFolder_temp = [];
    else
        saveFolder_temp = saveFolder;
    end
    
    % -- do the decoder analysis
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            axonRois = getAxonID_CM_v2(RespROIs.pearsonR{animal,pos},corrThresh);
            axonRois_stacked{animal,pos} = axonRois;
            Posterior{animal,pos} = BayesianDecoder_v5(RespROIs.data{animal,pos},training_fraction,RespROIs.V1az(animal,pos),0,axonRois,RespROIs.info,char(RespROIs.info.data_details{animal,pos}),dataset,saveFolder_temp);
        end
    end
    
    % -- Average'em together
    counter=0;
    Accuracy_dRF = NaN(n_animals,max(n_pos),13);
    Accuracy_dRF_shuffle = NaN(n_animals,max(n_pos),13);
    Accuracy = NaN(n_animals,max(n_pos),13);
    Accuracy_shuffle = NaN(n_animals,max(n_pos),13);
    Accuracy_Ele = NaN(n_animals,max(n_pos),3);
    Accuracy_Ele_shuffle = NaN(n_animals,max(n_pos),3);
    delta = NaN(n_animals,max(n_pos));
    delta_sh = NaN(n_animals,max(n_pos));
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            if ~isempty(Posterior{animal,pos})
                counter=counter+1; % just to count the number of sessions
                % % - Azimuth accuracy
                Accuracy(animal,pos,:) = Posterior{animal,pos}.collapsed(:,1);
                Accuracy_shuffle(animal,pos,:) = Posterior{animal,pos}.collapsed_shuffled(:,1);
                
                Accuracy_dRF(animal,pos,1:length(Posterior{animal,pos}.binsum)) = Posterior{animal,pos}.binsum;
                Accuracy_dRF_shuffle(animal,pos,1:length(Posterior{animal,pos}.binsum_shuffle)) = Posterior{animal,pos}.binsum_shuffle;
                
                % % - Elevation accuracy
                Accuracy_Ele(animal,pos,:) = Posterior{animal,pos}.Elecollapsed;
                Accuracy_Ele_shuffle(animal,pos,:) = Posterior{animal,pos}.Elecollapsed_shuffled;
                
                % % - Distance to actual stim position -- Aug2023
                delta(animal,pos) = mean(abs(Posterior{animal,pos}.deltaDecoded2Actual));
                delta_sh(animal,pos) = mean(abs(Posterior{animal,pos}.deltaDecoded2Actual_sh));
                
            end
        end
    end
    
    %% Distance to actual stim position - single session level
    % % - absolute speaker position, az
    accuracy_meanAcrossSessions = squeeze(nanmean(Accuracy,2));
    nNaNs = sum(isnan(accuracy_meanAcrossSessions));
    stdAccuracy = std(accuracy_meanAcrossSessions,[],1,'omitnan')./sqrt(n_animals-nNaNs);
    avgAccuracy = squeeze(nanmean(accuracy_meanAcrossSessions,1));
    
    accuracy_meanAcrossSessions_sh = squeeze(nanmean(Accuracy_shuffle,2));
    stdAccuracy_sh = std(accuracy_meanAcrossSessions_sh,[],1,'omitnan')./sqrt(n_animals-nNaNs);
    avgAccuracy_sh = squeeze(nanmean(accuracy_meanAcrossSessions_sh,1));
    
    % % - absolute speaker position, ele
    accuracyEle_meanAcrossSessions = squeeze(nanmean(Accuracy_Ele,2));
    nNaNs = sum(isnan(accuracyEle_meanAcrossSessions));
    stdAccuracyEle = std(accuracyEle_meanAcrossSessions,[],1,'omitnan')./sqrt(n_animals-nNaNs);
    avgAccuracyEle = squeeze(nanmean(accuracyEle_meanAcrossSessions,1));
    
    accuracyEle_meanAcrossSessions_sh = squeeze(nanmean(Accuracy_Ele_shuffle,2));
    stdAccuracyEle_sh = std(accuracyEle_meanAcrossSessions_sh,[],1,'omitnan')./sqrt(n_animals-nNaNs);
    avgAccuracyEle_sh = squeeze(nanmean(accuracyEle_meanAcrossSessions_sh,1));
    
    % % - stimulus position relative to V1RF
    accuracy_meanAcrossSessions_dRF = squeeze(nanmean(Accuracy_dRF,2));
    avgAccuracy_dRF = squeeze(nanmean(accuracy_meanAcrossSessions_dRF));
    nNaN=sum(isnan(accuracy_meanAcrossSessions_dRF));
    
    accuracy_meanAcrossSessions_dRF_sh = squeeze(nanmean(Accuracy_dRF_shuffle,2));
    avgAccuracy_dRF_sh = squeeze(nanmean(accuracy_meanAcrossSessions_dRF_sh));
    
    deltaRF = [0:10:120];
    stdAccuracy_dRF = std(accuracy_meanAcrossSessions_dRF,[],1,'omitnan')./sqrt(n_animals-nNaN);
    stdAccuracy_dRF_sh = std(accuracy_meanAcrossSessions_dRF_sh,[],1,'omitnan')./sqrt(n_animals-nNaN);
    
    % % - calculate chance
    x=-20:10:100; y = -20:20:20;
    [X,Y] = meshgrid(x,y);
    Xlin = reshape(X',39,1); Ylin = reshape(Y',39,1);
    for i = 1:39
        for ii = 1:39
            chance(i,ii) =  norm([Xlin(i),Ylin(i)]-[Xlin(ii),Ylin(ii)]);
        end
    end
    chance_avg = mean(chance,'all');
    temp = mean(chance);
    chance_2D = reshape(temp,13,3);
    
    % --  Plot'em together
    figure
    subplot(3,4,[1 2]); hold on  % azimuth error
    plot([-20:10:100],mean(chance_2D,2),'k:')
    fill([-20:10:100 100:-10:-20],[avgAccuracy_sh-stdAccuracy_sh flip(avgAccuracy_sh+stdAccuracy_sh)],...
        [0.5 .5 .5],'edgecolor','none','LineStyle','-','facealpha',0.5);
    fill([-20:10:100 100:-10:-20],[avgAccuracy-stdAccuracy flip(avgAccuracy+stdAccuracy)],...
        [0.6350 0.0780 0.1840],'edgecolor','none','LineStyle','-','facealpha',0.5);
    plot([-20:10:100],avgAccuracy,'r','Linewidth',2 )
    plot([-20:10:100],avgAccuracy_sh,'k','Linewidth',2 )
    ylim(yl); ylabel('Distance to Actual Stim pos (degrees)')
    xlim([-25 105]); xlabel('Stim Az Position (deg)')
    
    subplot(3,4,3); hold on % elevation error
    plot([-20:20:20],mean(chance_2D,1),'k:')
    fill([-20:20:20 20:-20:-20],[avgAccuracyEle_sh-stdAccuracyEle_sh...
        flip(avgAccuracyEle_sh+stdAccuracyEle_sh)],...
        [0.5 .5 .5],'edgecolor','none','LineStyle','-','facealpha',0.5);
    fill([-20:20:20 20:-20:-20],[avgAccuracyEle-stdAccuracyEle ...
        flip(avgAccuracyEle+stdAccuracyEle)],...
        [0.6350 0.0780 0.1840],'edgecolor','none','LineStyle','-','facealpha',0.5);
    plot(-20:20:20,avgAccuracyEle,'r','Linewidth',2 )
    plot(-20:20:20,avgAccuracyEle_sh,'k','Linewidth',2 )
    ylim(yl); %ylabel('Ele. Decoding Accuracy')
    xlim([-25 25]); xticks([-20:20:20]);xlabel('Stim Ele Position (deg)')
    
    subplot(3,4,4); hold on % overall error
    animal_avg = nanmean(delta,2);
    animal_avg_sh = nanmean(delta_sh,2);
    scatter(ones(n_animals,1),animal_avg,'filled','markeredgecolor','w')
    scatter(2*ones(n_animals,1),animal_avg_sh,'filled','markeredgecolor','w')
    plot([0.9 1.1],[nanmean(animal_avg) nanmean(animal_avg)],'k')
    plot([1.9 2.1],[nanmean(animal_avg_sh) nanmean(animal_avg_sh)],'k')
    n_mice_temp = sum(~isnan(animal_avg));
    n_sessions_temp = sum(~isnan(delta),'all');
    plot(repmat([1;2],1,n_animals),[animal_avg(:) animal_avg_sh(:)]','k')
    ylim(yl); % ylabel('Distance to Actual Stim pos (degrees)')
    xlim([.5 2.5])
    plot([1 2],[chance_avg chance_avg],'k:')
    [~,p_distance2actual] = ttest(nanmean(delta,2),nanmean(delta_sh,2));
    title(['p. ttest, p=' num2str(p_distance2actual,4)])
    
    subplot(3,4,[5 6]); hold on % error with respect to V1RF
    fill([deltaRF flip(deltaRF)],[avgAccuracy_dRF_sh-stdAccuracy_dRF_sh flip(avgAccuracy_dRF_sh+stdAccuracy_dRF_sh)],...
        [0.5 .5 .5],'edgecolor','none','LineStyle','-','facealpha',0.5);
    fill([deltaRF flip(deltaRF)],[avgAccuracy_dRF-stdAccuracy_dRF flip(avgAccuracy_dRF+stdAccuracy_dRF)],...
        [0.6350 0.0780 0.1840],'edgecolor','none','LineStyle','-','facealpha',0.5);
    p1 = plot(deltaRF,avgAccuracy_dRF_sh,'k','Linewidth',2 );
    p2 = plot(deltaRF,avgAccuracy_dRF,'r','Linewidth',2 );
    ylim(yl); ylabel('Distance to Actual Stim pos (degrees)')
    xlim([-5 125]); xlabel('|Stim-V1RF| distance')
    legend([p1,p2],{'shuffle','data'},'location','nw')
    title('dRF')
    
    subplot(3,4,[7 8]); hold on % error with respect to V1RF
    dRF_sub = accuracy_meanAcrossSessions_dRF-accuracy_meanAcrossSessions_dRF_sh;
    dRF_sub_avg = nanmean(dRF_sub);
    dRF_sub_std = std(dRF_sub,'omitnan')./sqrt(n_mice_temp);
    plot(deltaRF,zeros(size(deltaRF)),'k:' );
    fill([deltaRF flip(deltaRF)],[dRF_sub_avg-dRF_sub_std flip(dRF_sub_avg+dRF_sub_std)],...
        [0.6350 0.0780 0.1840],'edgecolor','none','LineStyle','-','facealpha',0.5);
    plot(deltaRF,dRF_sub_avg,'r','Linewidth',2 );
    ylim([-40 10]); ylabel('delta(Distance to Actual Stim pos) (degrees)')
    xlim([-5 125]); xlabel('|Stim-V1RF| distance')
    
    % % - rm 1w ANOVA
    
    meas = dRF_sub;
    t = table(meas(:,1),meas(:,2),meas(:,3),meas(:,4),meas(:,5),meas(:,6),meas(:,7),meas(:,8),meas(:,9),meas(:,10),...
        'VariableNames',{'meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10'});
    Meas = table([1:10]','VariableNames',{'Measurements'});
    rm = fitrm(t,'meas1-meas10~1');
    ranovatbl = ranova(rm);
    p = table2array(ranovatbl(1,5));
    titleText = sprintf('1w RM ANOVA, p=%.3g, F(%g,%g)=%.3g',p,ranovatbl{1,2},ranovatbl{2,2},ranovatbl{1,4});
    title(titleText)
    
    % % - infos
    subplot(3,4,9)
    text(0,0,['n=' num2str(n_mice_temp) 'mice, ' num2str(n_sessions_temp) ' pos.'],'horizontalalignment','left','verticalalignment','top','FontSize',8);
    axis off
    
    % % --  save azimtuh accuracy fig
    sgtitle(['Average across positions' newline '\itROIs axonized, corr th=' num2str(corrThresh)])
    saveas(gcf,[saveFolder dataset '_avgAcrossSessions.tif']);
    saveas(gcf,[saveFolder dataset '_avgAcrossSessions.svg']);
    savefig([saveFolder dataset '_avgAcrossSessions']);
    pause(0.1); close
    
    % % - save data for "source data"
    DecoderSourceData.az = accuracy_meanAcrossSessions;
    DecoderSourceData.az_sh = accuracy_meanAcrossSessions_sh;
    DecoderSourceData.ele = accuracyEle_meanAcrossSessions;
    DecoderSourceData.ele_sh = accuracyEle_meanAcrossSessions_sh;
    DecoderSourceData.perMouse = delta;
    DecoderSourceData.perMouse_sh = delta_sh;
    DecoderSourceData.dRF = accuracy_meanAcrossSessions_dRF;
    DecoderSourceData.dRF_sh = accuracy_meanAcrossSessions_dRF_sh;
    save([saveFolder 'DecoderSourceData'],'DecoderSourceData')
    
    %% Distribution of the decoding error
    temp = NaN(n_animals,max(n_pos),14);
    temp_sh = NaN(n_animals,max(n_pos),14);
    figure; hold on
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            if ~isempty(Posterior{animal,pos})
                distRel = abs(Posterior{animal,pos}.deltaDecoded2Actual);
                distRel_sh = abs(Posterior{animal,pos}.deltaDecoded2Actual_sh);
                N = histcounts(distRel,-5:10:135); % 10 degrees bin
                temp(animal,pos,:) = N;
                N = histcounts(distRel_sh,-5:10:135);
                temp_sh(animal,pos,:) = N;
            end
        end
        toPlot = squeeze(nanmean(temp(animal,:,:),2))-squeeze(nanmean(temp_sh(animal,:,:),2));
        p1=plot(toPlot./780);
        set(p1,'Color',[1 0 0 0.2]);
    end

    toPlot = squeeze(mean(nanmean(temp,2)))-squeeze(mean(nanmean(temp_sh,2)));
    plot(toPlot./780,'r');
    
    plot([0 14],[0 0],'k:')
    xlim([0 15]); xlabel('Relative decoded dist. to actual sound source location')
    xticks([1:14]);xticklabels([0:10:130])
    ylim([-.1 .25])
    ylabel('fraction decoded')
    
    % % 1W RM ANOVA - effect of distance
    varNames_az = cell(1,14);
    for i=1:14
        varNames_az{i} = ['az',num2str(i)];
    end
    data_anova = [squeeze(nanmean(temp,2))-squeeze(nanmean(temp_sh,2))];
    
    t = table(data_anova(:,1),data_anova(:,2),data_anova(:,3),data_anova(:,4),data_anova(:,5),data_anova(:,6),data_anova(:,7),data_anova(:,8),data_anova(:,9),data_anova(:,10),data_anova(:,11),data_anova(:,12),data_anova(:,13),data_anova(:,14),...
        'VariableNames',varNames_az);
    rm = fitrm(t,'az1-az14~1');
    rmTable1 = ranova(rm);
    title({['Distribution of decoding error'],['1W RM ANOVA,int: F(' num2str(table2array(rmTable1(1,2))) ',' ,num2str(table2array(rmTable1(2,2))) ')=' num2str(table2array(rmTable1(1,4)),4) ,...
        ', p =',num2str(table2array(rmTable1(1,5)),4)]})
    
    % % --  save hist of decoding errors
    saveas(gcf,[saveFolder dataset '_DistributionDecodingError.tif']);
    saveas(gcf,[saveFolder dataset '_DistributionDecodingError.svg']);
    savefig([saveFolder dataset '_DistributionDecodingError']);
    pause(0.1); close
    
    %%  decoding accuracy as a function of V1 - first need to use the same number of axons
    if doAccVsV1
        nperm = 10;        % 10 used in Mazo et al., 2024
        nAxons_min = 20;   % 20 used in Mazo et al., 2024
        disp('Calculating decoding accuracy as a function of V1 retinotopy, will take a few minutes...')
        for animal = 1:n_animals
            for pos = 1:n_pos(animal)
                axonRois = getAxonID_CM_v2(RespROIs.pearsonR{animal,pos},corrThresh);
                axonRois_stacked{animal,pos} = axonRois;
            end
        end
        numAxons_perSession = cellfun(@length,axonRois_stacked);
        numAxons_perSession = numAxons_perSession(:);numAxons_perSession(numAxons_perSession<nAxons_min)=[];
        min_numAxons = min(numAxons_perSession);
        
        Posterior = cell(n_animals,max(n_pos),10);
        for animal = 1:n_animals
            for pos = 1:n_pos(animal)
                nROIs = size(RespROIs.data{animal,pos},1);
                axonRois = axonRois_stacked{animal,pos};
                if length(axonRois)>=nAxons_min
                    RepROI = NaN(1,length(axonRois));
                    for i = 1:length(axonRois) % loop through positions
                        axonRois{1,i}(axonRois{1,i}==0)=[];
                        if length(axonRois{1,i})>1
                            meandFF = mean(mean(mean(RespROIs.data{animal,pos}(axonRois{1,i},:,:,:),2),3),4);
                            [~,idx_maxdFF] = max(meandFF);
                            RepROI(i) = axonRois{1,i}(idx_maxdFF); % use the one that has overall the highest fluo
                        else
                            RepROI(i) = axonRois{1,i};
                        end
                    end
                    
                    for i=1:nperm
                        randAxonSelec = randperm(length(RepROI),min_numAxons);
                        selected_RepRois = false(1,nROIs);
                        selected_RepRois(randAxonSelec)=true;
                        
                        data = RespROIs.data{animal,pos}(selected_RepRois,:,:,:);
                        Posterior{animal,pos,i} = BayesianDecoder_v5(data,training_fraction,RespROIs.V1az(animal,pos),0,[],RespROIs.info,char(RespROIs.info.data_details{animal,pos}),dataset,[]);
                    end
                end
            end
        end
        
        % -- Average'em together
        c=0;
        Accuracy = NaN(n_animals,max(n_pos),nperm);
        for animal = 1:n_animals
            for pos = 1:n_pos(animal)
                for i = 1:nperm
                    if ~isempty(Posterior{animal,pos})
                        c=c+1;
                        Accuracy(animal,pos,i) = nanmean(Posterior{animal,pos,i}.collapsed(:,1));
                    end
                end
            end
        end
        
        % % - calculate averages and plot
        avgAccuracy = nanmean(Accuracy,3);
        V1az_all = RespROIs.V1az;
        V1az_all(V1az_all>=100)=99;
        
        figure; hold on
        plot([0 100],[50.42 50.42],'k:') % - chance
        % % - plot the data
        for i = 1:n_animals
            scatter(V1az_all(i,:),avgAccuracy(i,:),'filled','markerfacecolor',[.5 .5 .5],'markeredgecolor','w')
        end
        
        % % - calculate binned averages
        avgAccuracy_binned = NaN(4,n_animals,6);
        for j = 1:4
            V1azbins=[0:20:100]+5*(j-1);
            for i = 1:n_animals
                [~,~,bin]=histcounts(V1az_all(i,:),V1azbins);
                is_member = NaN(length(V1azbins),length(bin));
                for ii = 1:length(V1azbins)
                    tf = ismember(bin',ii);
                    is_member(ii,:) = tf;
                    avgAccuracy_binned(j,i,ii) = nanmean(avgAccuracy(i,tf));
                end
            end
        end
        
        % % - merge the averages from different binning in a 2D matrix
        movAvg = NaN(n_animals,4*6);
        for j= 1:4
            movAvg(:,1+(j-1):4:4*6)=squeeze(avgAccuracy_binned(j,:,:));
        end
        
        % % - plot
        movAvg_mean = nanmean(movAvg); movAvg_mean(isnan(movAvg_mean))=[];
        nNaNs = sum(isnan(movAvg));
        sd = std(movAvg,'omitnan')./sqrt(n_animals-nNaNs); sd(isnan(sd))=[];
        deltaAz = [10:5:125]; az_vector = deltaAz(1:length(movAvg_mean));
        plot(az_vector,movAvg_mean,'k')
        fill([az_vector flip(az_vector)],[movAvg_mean-sd flip(movAvg_mean+sd)],[0 0 0],'edgecolor','none','facealpha',0.2)
        
        % % - 1w anova with nan => anovan
        deltaAz=movAvg(:,1:20); movAvg_reshape = reshape(deltaAz,n_animals*size(deltaAz,2),1);
        groups=[];
        for  i = 1:size(deltaAz,2)
            groups(:,i)=[i*ones(n_animals,1)];
        end
        g2=reshape(groups,n_animals*size(deltaAz,2),1);
        [p,tbl,stats]=anovan(movAvg_reshape,g2,'model','interaction','varnames',{'az'},'display','off');
        if p<0.05
            c = multcompare(stats,'display','off');
            k = find(c(:,6)<0.05);
            if ~isempty(k)
                writematrix(k,[saveFolder dataset '_avgAcrossV1pos.txt']);
                posthoc = 1;
            else
                posthoc = 0;
            end
            textTitle = sprintf('1W ANOVA F(%.3g,%.3g)=%.3g, p=%.3g | posthoc = %.3g',cell2mat(tbl(2,3)),cell2mat(tbl(3,3)),cell2mat(tbl(2,6)),p,posthoc);
        else
            textTitle = sprintf('1W ANOVA F(%.3g,%.3g)=%.3g, p=%.3g \n n=%g axon/session used',cell2mat(tbl(2,3)),cell2mat(tbl(3,3)),cell2mat(tbl(2,6)),p,min_numAxons);
        end
        
        % % - finalize and save
        xlabel('V1 azimuth')
        ylim(yl); ylabel('Mean distance 2 actual'); yticks([0 .5 1])
        title(textTitle)
        
        saveas(gcf,[saveFolder dataset '_avgAcrossV1pos.tif']);
        saveas(gcf,[saveFolder dataset '_avgAcrossV1pos.svg']);
        savefig([saveFolder dataset '_avgAcrossV1pos']);
        pause(0.1); close
    end
end % if doSingleSession

%% Concatenate positions of each mouse
if doMouseagregate
    
    % % - mega virtual mouse, aka Mickey Mouse
    disp('calculating decoding accuracy = f(n_axons)...')
    % - axonize the data
    finalData = cell(n_animals,max(n_pos));
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            if animal == 3 && pos == 5 && TwentykHz == 1
                disp('do not use CMad56 pos6') % because there is only 16 rep/trial type -- it messes up the decoder (fraction of trials used to generate the model)
                axonRois_stacked{animal,pos} = [];
            else
                clear RpzROI
                session_data = RespROIs.data{animal, pos}(:,1:30,:,:) ;
                SameAxonROIs = getAxonID_CM_v2(RespROIs.pearsonR{animal,pos},corrThresh);
                for i = 1:length(SameAxonROIs)
                    if length(SameAxonROIs{1,i})>1
                        meandFF = mean(mean(mean(session_data(SameAxonROIs{1,i},:,:,:),2),3),4);
                        [~,idx_maxdFF] = max(meandFF);
                        RpzROI(i) = SameAxonROIs{1,i}(idx_maxdFF); % use the one that has overall the highest fluo
                    else
                        RpzROI(i) = SameAxonROIs{1,i};
                    end
                end
                finalData{animal,pos} = session_data(RpzROI,:,:,:);
            end
        end % loop thru pos
    end % loop thru mice
    
    % - concatenate all the data
    data_all=cat(1,finalData{:});
    ideal_nAxons = [10 25,50,100,200,300];
    if size(data_all,1)>300
        nAxons = ideal_nAxons;
    else
        nAxons = 10; next_nAxons = ideal_nAxons(2);
        c = 2;
        while next_nAxons<size(data_all,1)
            nAxon_temp = ideal_nAxons(c);
            nAxons = [nAxons ideal_nAxons(c)];
            c = c+1;
            next_nAxons = ideal_nAxons(c);
        end
        nAxons = [nAxons size(data_all,1)];
    end
    
    nDraws = 10; % 100 in Mazo et al., 2024
    counter = 0;
    for i = 1:length(nAxons)
        for ii = 1:nDraws
            counter=counter+1;
            y = datasample(data_all,nAxons(i),'Replace',false);
            Posterior_MickeyMouse{i,ii} = BayesianDecoder_v5(y,training_fraction,NaN,0,[],RespROIs.info,'MickeyMouse',dataset,[]);
        end
    end
    
    % % - deocding error as a function of number of axons
    delta = NaN(length(nAxons),nDraws);
    for i = 1:length(nAxons)
        for ii = 1:nDraws
            delta(i,ii) = mean(abs(Posterior_MickeyMouse{i,ii}.deltaDecoded2Actual));
        end
    end
    
    % % - Average and CI across rep
    mean_deltaAz = mean(delta,2);
    CI_deltaAz = prctile(delta,[5 95],2);
    
    figure; hold on
    fill([nAxons flip(nAxons)],[CI_deltaAz(:,1); flip(CI_deltaAz(:,2))],'b','edgecolor','none','facealpha',.2)
    plot(nAxons,mean_deltaAz,'b')
    ylim([0 100])
    xticks(nAxons)
    
    x=-20:10:100; y = -20:20:20;
    [X,Y] = meshgrid(x,y);
    Xlin = reshape(X',39,1); Ylin = reshape(Y',39,1);
    for i = 1:39
        for ii = 1:39
            chance(i,ii) =  norm([Xlin(i),Ylin(i)]-[Xlin(ii),Ylin(ii)]);
        end
    end
    chance_avg = mean(chance,'all');
    
    plot([nAxons(1) nAxons(end)],[chance_avg chance_avg],'k:')
end

%% save
if doMouseagregate
    Decoder_Output.MickeyMouse = Posterior_MickeyMouse;
    Decoder_Output.delta = delta;
    Decoder_Output.AcrossMice = Accuracy;
    Decoder_Output.AcrossMice_sh = Accuracy_shuffle;
else
    Decoder_Output.delta = delta;
    Decoder_Output.delta_sh = delta_sh;
end
end
function Decoder_Output = BayesianDecoder_batch_RotSPK_v3(RespROIs,MainFolder,az_vector,doSingleSession,doAnimalData,do_MickeyMouse,PlotAnimalData,plotIndivSession)
%% a few parameters to define
corrThresh = 0.3;         % to determine which ROIs belong to the same axon
training_fraction = 0.8;

%% a few definitions
nAzPos = length(az_vector);
n_animals = size(RespROIs.info.animalID,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);
saveFolder = MainFolder;

%% do Bayesian decoder on single session (and plot and save)
if doSingleSession
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            axonRois = getAxonID_CM_v2(RespROIs.pearsonR{animal,pos},corrThresh);
            axonRois_stacked{animal,pos} = axonRois;
            Posterior{animal,pos} = BayesianDecoder_v4_RotSPK(RespROIs.data{animal,pos},training_fraction,RespROIs.V1az(animal,pos),0,axonRois,RespROIs.info,RespROIs.info.data_details{animal,pos},az_vector,saveFolder,plotIndivSession);
        end
    end
    
    %% decoding error (in 2D)
    chance = NaN(length(az_vector),length(az_vector));
    for i = 1:length(az_vector)
        for ii = 1:length(az_vector)
            chance(i,ii) = norm(az_vector(i)-az_vector(ii));
        end
    end
    chance_avg = mean(chance);
    
    % % - calculate the errors
    mean_error = NaN(n_animals,length(az_vector));
    mean_error_sh = NaN(n_animals,length(az_vector));
    sem_error = NaN(n_animals,length(az_vector));
    mean_error_all = NaN(n_animals,max(n_pos),length(az_vector)); % quick fix for dRF figure instead of re-writting
    mean_error_all_sh = NaN(n_animals,max(n_pos),length(az_vector));
    for animal = 1:n_animals
        temp = NaN(length(az_vector),n_pos(animal));
%         temp_sh = NaN(length(az_vector),n_pos(animal));
        for pos = 1:n_pos(animal)
            if ~isempty(Posterior{animal,pos})
                if isequal(az_vector,-90:20:90)
                    data = abs(Posterior{animal,pos}.delta_decod2actual_error).*20;
                elseif isequal(az_vector,-10:10:100)
                    data = abs(Posterior{animal,pos}.delta_decod2actual_error).*10;
                end
                temp(:,pos) = mean(data,2);
            end
        end
        mean_error(animal,:) = nanmean(abs(temp),2);
        sem_error(animal,:) = std(temp,[],2,'omitnan')./sqrt(n_pos(animal));
        mean_error_all(animal,1:pos,:) = temp';
    end
    
    %% final figure - decoding error
    figure;
    subplot(1,3,1); hold on
    plot(az_vector,chance_avg,'k:')
    for i = 1:n_animals
        plot(az_vector,mean_error(i,:),'color',[.5 .5 .5]);%,'k')
    end
    plot(az_vector,mean(mean_error),'k','linewidth',2)
    ylabel('Decoding error (deg)')
    
    subplot(1,3,2); hold on
    for i = 1:n_animals
        plot(az_vector,mean_error(i,:)-chance_avg,'color',[.5 .5 .5]);%,'k')
    end
    plot(az_vector,mean(mean_error)-chance_avg,'k','linewidth',2)
    ylabel('Decoding error relative to chance (deg)')
    
    if isequal(az_vector,-90:20:90)
        subplot(1,3,3); hold on
        x = mean_error-chance_avg;
        rel_error(:,1) = mean(x(:,[1:3]),2);
        rel_error(:,2) = mean(x(:,[4:7]),2);
        rel_error(:,3) = mean(x(:,[8:10]),2);
        for i = 1:n_animals
            plot([1:3],[rel_error(i,:)],'k')
        end
        scatter(ones(n_animals,1),rel_error(:,1),'k','filled','markeredgecolor','none')
        scatter(2*ones(n_animals,1),rel_error(:,2),'k','filled','markeredgecolor','none')
        scatter(3*ones(n_animals,1),rel_error(:,3),'k','filled','markeredgecolor','none')
        plot([.9 1.9 2.9;1.1 2.1 3.1],...
            [mean(rel_error(:,1)) mean(rel_error(:,2)) mean(rel_error(:,3));mean(rel_error(:,1)) mean(rel_error(:,2)) mean(rel_error(:,3))],...
            'k','linewidth',3)
        xlim([.5 3.5]);xticks([1:3]);xticklabels({'Ipsi','Front','Contra'});
        t = table(rel_error(:,1),rel_error(:,2),rel_error(:,3),...
            'VariableNames',{'meas1','meas2','meas3'});
        rm = fitrm(t,'meas1-meas3~1');
        ranovatbl = ranova(rm);
        p = table2array(ranovatbl(1,5));
        c = multcompare(rm,'Time');
        title({['RM 1W ANOVA,p=' num2str(p,4)],['1vs2,p=' num2str(table2array(c(1,5),4)) '|2vs3,p=' num2str(table2array(c(4,5),4))]})
        
    elseif isequal(az_vector,-10:10:100)
        subplot(1,3,3); hold on
        meanPermice = mean(mean_error,2);
        meanPermice_sh = mean(mean_error_sh,2);
        plot(repmat([1 2]',1,2),[meanPermice meanPermice_sh]','k')
        scatter(ones(n_animals,1),meanPermice,'r','filled','markeredgecolor','none')
        scatter(2*ones(n_animals,1),meanPermice_sh,'k','filled','markeredgecolor','none')
        plot([1 2],[mean(chance_avg) mean(chance_avg)],'k:')
        ylim([20 50])
        xlim([.5 2.5]); xticks([1 2])
    end
    
    set(gcf,'units','normalized','position',[.2 .1 .7 .4])
    saveas(gcf,[saveFolder filesep 'DecodingError.tif']);
    saveas(gcf,[saveFolder filesep 'DecodingError.svg']);
    savefig([saveFolder filesep 'DecodingError']);
    
end % if doSingleSession

%% Animal agregate
if doAnimalData
    clear data
    for animal = 1:n_animals
        for pos = 1:n_pos(animal)
            data{animal,pos} = RespROIs.data{animal, pos}(:,1:30,:,:) ;
            
            axonRois = getAxonID_CM_v2(RespROIs.pearsonR{animal,pos},corrThresh);
            axonRois_stacked{animal,pos} = axonRois;
            nROIs(animal,pos) = size(RespROIs.pearsonR{animal,pos},1);
        end
    end
    for animal = 1:n_animals
        nROIs_cumsum = cumsum(nROIs(animal,:));
        for pos = 2:max(n_pos(animal))
            if ~isempty(axonRois_stacked{animal,pos})
                for axon = 1:size(axonRois_stacked{animal,pos},2)
                    axonRois_stacked{animal,pos}{1,axon} = axonRois_stacked{animal,pos}{1,axon}+nROIs_cumsum(pos-1);
                end
            else
                disp([num2str(animal) ',' num2str(pos)])
            end
        end
        axonRois_permice{animal} = cat(2,axonRois_stacked{animal,:});
    end
    
    data_perAnimal = cell(1,n_animals);
    Posterior_perAnimal = cell(n_animals,1);
    for animal = 1:n_animals
        data_perAnimal{animal} = cat(1,data{animal, :});
        Posterior_perAnimal{animal} = BayesianDecoder_v4_RotSPK(data_perAnimal{animal},training_fraction,0,0,axonRois_permice{animal} ,RespROIs.info,RespROIs.info.animalID{animal},az_vector,saveFolder,PlotAnimalData);
    end
    
    % % -- spatial aspect of decoding
    for animal = 1:n_animals
        Accuracy_deltaToActual(animal,:) = Posterior_perAnimal{animal}.delta_decod2actual;
        Accuracy_deltaToActual_sh(animal,:) = Posterior_perAnimal{animal}.delta_decodShuff2actual;
        Accuracy_deltaToActual_abs(animal,:) = Posterior_perAnimal{animal}.delta_decod2actual_abs;
        Accuracy_deltaToActual_sh_abs(animal,:) = Posterior_perAnimal{animal}.delta_decodShuff2actual_abs;
        [N,edges1] = histcounts(Accuracy_deltaToActual(animal,:),[-10.5:10.5]);
        data_N(animal,:) = N;
        [N,edges] = histcounts(Accuracy_deltaToActual_sh(animal,:),[-10.5:10.5]);
        shuffle_N(animal,:) = N;
        [N,edges] = histcounts(Accuracy_deltaToActual_abs(animal,:),[-0.5:12.5]);
        data_N_abs(animal,:) = N;
        [N,edges] = histcounts(Accuracy_deltaToActual_sh_abs(animal,:),[-0.5:12.5]);
        shuffle_N_abs(animal,:) = N;
    end
    
    % -- plot and stats
    figure;
    subplot(1,2,1); hold on
    p1=plot([-10:10],[data_N./780]');
    set(p1,'Color',[1 0 0 0.2]);
    p2= plot([-10:10],[shuffle_N./780]');
    set(p2,'Color',[0 0 0 0.2]);
    plot([-10:10],mean(data_N)./780,'r','linewidth',3)
    plot([-10:10],mean(shuffle_N)./780,'k','linewidth',3)
    xlim([-10 10]); xlabel('Dist. to actual sound source location')
    ylabel('fraction decoded')
    
    subplot(1,2,2); hold on
    p1=plot([data_N_abs./780]');
    set(p1,'Color',[1 0 0 0.2]);
    p2= plot([shuffle_N_abs./780]');
    set(p2,'Color',[0 0 0 0.2]);
    plot(mean(data_N_abs)./780,'r','linewidth',3)
    plot(mean(shuffle_N_abs)./780,'k','linewidth',3)
    xlim([0 12]); xlabel('Relative dist. to actual sound source location')
    ylabel('fraction decoded')
    
    % % - anova attempt
    % % RM ANOVA
    groups = [ones(n_animals,1);2*ones(n_animals,1)];
    groups = num2str(groups);
    varNames_az = cell(1,14);
    for i=2:14
        varNames_az{i} = ['az',num2str(i-1)];
    end
    varNames_az{1} = 'Group';
    data_anova = [data_N_abs;shuffle_N_abs];
    
    t = table(groups,data_anova(:,1),data_anova(:,2),data_anova(:,3),data_anova(:,4),data_anova(:,5),data_anova(:,6),data_anova(:,7),data_anova(:,8),data_anova(:,9),data_anova(:,10),data_anova(:,11),data_anova(:,12),data_anova(:,13),...
        'VariableNames',varNames_az);
    Meas = table([1:13]','VariableNames',{'Az'});
    rm = fitrm(t,'az1-az13~Group','WithinDesign',Meas);
    rmTable1 = ranova(rm);
    
    varNames_az = cell(1,22);
    for i=2:22
        varNames_az{i} = ['az',num2str(i-1)];
    end
    varNames_az{1} = 'Group';
    data_anova = [data_N;shuffle_N];
    
    t = table(groups,data_anova(:,1),data_anova(:,2),data_anova(:,3),data_anova(:,4),data_anova(:,5),data_anova(:,6),data_anova(:,7),data_anova(:,8),data_anova(:,9),data_anova(:,10),data_anova(:,11),data_anova(:,12),data_anova(:,13),...
        data_anova(:,14),data_anova(:,15),data_anova(:,16),data_anova(:,17),data_anova(:,18),data_anova(:,19),data_anova(:,20),data_anova(:,21),...
        'VariableNames',varNames_az);
    Meas = table([1:21]','VariableNames',{'Az'});
    rm = fitrm(t,'az1-az21~Group','WithinDesign',Meas);
    rmTable2 = ranova(rm);
    
    
    subplot(1,2,2); title(['RM ANOVA, int. p = ' num2str(rmTable1{2,6},3)]);
    subplot(1,2,1); title(['RM ANOVA, int. p = ' num2str(rmTable2{2,6},3)]);
    
    saveas(gcf,[saveFolder '_SpatialDecoding_hist.tif']);
    saveas(gcf,[saveFolder '_SpatialDecoding_hist.svg']);
    savefig([saveFolder '_SpatialDecoding_hist']);
    pause(0.1); close
    
    % % average the log likelihood across mice
    allPosteriors = NaN(nAzPos,nAzPos,n_animals);
    for i = 1:n_animals
        allPosteriors(:,:,i) = Posterior_perAnimal{i}.avg_PerStimType;
    end
    
    figure
    test2 = mean(allPosteriors,3);
    y_ax = 1;
    for i = 1:nAzPos
        test = test2(:,i);
        [~,I]=max(test,[],'all','linear');
        subtightplot(nAzPos,1,i);hold on;
        imagesc(test');colormap gray
        x_ax = i;
        plot([x_ax-.5 x_ax+.5 x_ax+.5 x_ax-.5 x_ax-.5],[y_ax-.5 y_ax-.5 y_ax+.5 y_ax+.5 y_ax-.5],'r-')
        decodedPos_x = I;
        plot(decodedPos_x,1,'o','MarkerSize',1)
        set(gca,'ydir','reverse')
        axis equal;   ylim([0.5 1.5]); xlim([0.5 16.5])
        xticks(''); yticks('')
    end
    set(gcf, 'units', 'normalized','position',[.3 .05 .4 .9])
    sgtitle(['Posterior probabilities, avged across mice'],'interpreter','none')
    
    saveas(gcf,[saveFolder '_SpatialDecoding.tif']);
    saveas(gcf,[saveFolder '_SpatialDecoding.svg']);
    savefig([saveFolder '_SpatialDecoding']);
    pause(0.1); close
end

%% mega virtual mouse, aka Mickey Mouse
if do_MickeyMouse
% - axonize the data
finalData = cell(n_animals,max(n_pos));
temp = cellfun('size',RespROIs.data,2);
temp(temp==0)=[];
min_TimeVectorLength = min(temp,[],'all');
for animal = 1:n_animals
    for pos = 1:n_pos(animal)
            clear RpzROI
            session_data = RespROIs.data{animal, pos}(:,1:min_TimeVectorLength,:,:) ;
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
    end % loop thru pos
end % loop thru mice
                   
% - concatenate all the data
data_all=cat(1,finalData{:});

ideal_nAxons = [10 25 50 100 150 200 300];
nAxons = 100;
c=1;
while ideal_nAxons(c+1) < size(data_all,1)
    c = c+1;
    nAxons = [nAxons ideal_nAxons(c)];
end
% if size(data_all,1)>500
%     nAxons = [nAxons ,500,1000];
% end
nDraws = 10; % 100 in Mazo et al., 2024
counter = 0;
for i = 1:length(nAxons)
    for ii = 1:nDraws
        counter=counter+1;
        y = datasample(data_all,nAxons(i),'Replace',false);
        Posterior_MickeyMouse{i,ii} = BayesianDecoder_v4_RotSPK(y,training_fraction,NaN,0,[],RespROIs.info,'MickeyMouse',az_vector,[],0);
    end
end


%% decoding error as a function of number of axons
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
ylabel('Decoding error'); ylim([0 100])
xlabel('Number of axons'); xticks(nAxons)
end

%% save
if do_MickeyMouse
Decoder_Output.MickeyMouse = Posterior_MickeyMouse;
Decoder_Output.delta = delta;
end
if doSingleSession
    Decoder_Output.AcrossMice = Accuracy;
    Decoder_Output.AcrossMice_sh = Accuracy_shuffle;
end
end
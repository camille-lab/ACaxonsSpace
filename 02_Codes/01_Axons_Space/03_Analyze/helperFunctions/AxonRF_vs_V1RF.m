function [RF_perpos,GOF_perpos,hasRF_perpos] = AxonRF_vs_V1RF(RespROIs,rsquare_th,Boutons,saveFolder,dataset)
% Boutons = 1, boutons; 0, axons
fprintf(1,'calculating RF center...')

%% some definition and variables initiations
n_animals = size(RespROIs.info.animalID ,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);

base_window = RespROIs.info.base_window;
try
    resp_window = RespROIs.info.resp_window{1,1};
catch
    resp_window = RespROIs.info.resp_window;
end

color_plot = distinguishable_colors(n_animals);

timeVect = 0:1/6.0962:7;

xData = -20:10:100;

mean_RF = NaN(n_animals,max(n_pos));
nRespBoutons = NaN(n_animals,max(n_pos));
RF_perpos = cell(n_animals,max(n_pos));
GOF_perpos = cell(n_animals,max(n_pos));
hasRF_perpos = cell(n_animals,max(n_pos));

RF_permice = cell(1,n_animals);
GOF_permice = cell(1,n_animals);

%% main loop - calculate RF az center
counter = 0;
for animal = 1:n_animals % loop through the animals
    for pos = 1:n_pos(animal) % loop through the sessions
        data = RespROIs.data{animal,pos};
        nROIs = size(data,1);
        
        tBase = false(size(timeVect));
        tBase(timeVect>=base_window(1) & timeVect<=base_window(2))=true;
        tBase = tBase(1:size(data,2));
        tAna = false(size(timeVect));
        tAna(timeVect>=resp_window(1) & timeVect<=resp_window(2))=true;
        tAna = tAna(1:size(data,2));
        
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
        
        RF = zeros(1,nROIs); GOF = zeros(1,nROIs);hasRF = false(1,nROIs);
        for ROI = 1:nROIs
            temp = squeeze(mean(mean(dFF(ROI,tAna,:,:),4),2));
            temp2 = reshape(temp,13,3);
            yData = mean(temp2,2);
            minResp = min(yData);
            if minResp < 0
                yData = yData-minResp;
            end
            
            [maxResp,~]=max(yData);
            yData = double(yData);
            [fitresult, gof]  = fit(xData',yData,'gauss1','Lower',[0 -20 0],'Upper',[maxResp+0.2*maxResp 100 Inf]);
            
            RF(ROI) = fitresult.b1;
            GOF(ROI) = gof.adjrsquare;
            
            if gof.adjrsquare > rsquare_th
                hasRF(ROI) = true;
            end
        end % end loop through ROIs
        
        if Boutons
            % - transform ROIs into boutons
            nBoutons_max = max(RespROIs.nBoutonsPerROI{animal,pos});
            for BoutonPerROI = 2:nBoutons_max
                multipleBoutons = find(RespROIs.nBoutonsPerROI{animal,pos}==BoutonPerROI);
                RF = [RF repmat(RF(multipleBoutons),1,BoutonPerROI-1)];
                GOF = [GOF repmat(RF(multipleBoutons),1,BoutonPerROI-1)];
                hasRF = [hasRF repmat(hasRF(multipleBoutons),1,BoutonPerROI-1)];
            end
            ROItype='boutons';
        elseif Boutons == 0
            corrThresh = 0.3;
            SameAxonROIs = getAxonID_CM_v2(RespROIs.pearsonR{animal,pos},corrThresh);
            temp = cellfun(@length,SameAxonROIs);
            axonWMultipleROI = find(temp>1);
            for i = axonWMultipleROI
                    meanF = mean(mean(mean(data(SameAxonROIs{1,i},:,:,:),2),3),4);
                    [~,idx_maxF] = sort(meanF,'descend');
                    bestWRF = find(hasRF(SameAxonROIs{1,i}(idx_maxF)),1);
                    
                    hasRF(SameAxonROIs{1,i}(idx_maxF(~bestWRF))) = false; % only keep the 'hasRF' of the best ROI
                    
                    if sum(hasRF(SameAxonROIs{1,i}))>1
                        counter = counter+1;
                        HasRF=hasRF(SameAxonROIs{1,i});
                        RF_sameAxon{counter} = RF(SameAxonROIs{1,i}(HasRF));
                    end
            end
            ROItype='axons';
        else
            ROItype='ROIs';
        end
        
        mean_RF(animal,pos) = mean(RF(hasRF));
        nRespBoutons(animal,pos) = sum(hasRF);
        
        RF_perpos{animal,pos}=RF(hasRF);
        GOF_perpos{animal,pos}=GOF(hasRF);
        hasRF_perpos{animal,pos} = hasRF;
        
    end % end loop through sessions
    RF_permice{animal} = cat(2,RF_perpos{animal,:});
    GOF_permice{animal} = cat(2,GOF_perpos{animal,:});
end %end loop through mice

%% plot RF from different ROIs of the same axon
if Boutons == 0
   for i = 1:length(RF_sameAxon)
       dRF{i,:} = diff(RF_sameAxon{i});       
   end

deltaRF = cat(2,dRF{:});
figure;histogram(abs(deltaRF),[0:5:120],'Normalization','probability');
xlabel('delta RF between ROIs of the same axon (degrees)')
ylabel('probability')

saveas(gcf,[saveFolder filesep dataset '_' ROItype '_RFsROIsSameAxon.tif']);
savefig([saveFolder filesep dataset '_' ROItype '_RFsROIsSameAxon']);
end

%%  plot bouton RF = f(V1RF)
width = 5;
isGoodPos = false(n_animals,max(n_pos));
ScatterFig = figure; subplot(1,4,[1 3]); hold on
indivAnimalFig = figure;
for animal = 1:n_animals % loop through the animals
    for pos = 1:n_pos(animal) % loop through the sessions
        nBoutons = length(RF_perpos{animal,pos});
        if nBoutons > 10
            isGoodPos(animal,pos) = true;
            jitter = (rand(nBoutons,1)-0.5)*width;
            figure(ScatterFig)
            scatter(RespROIs.V1az(animal,pos)*ones(length(nBoutons),1)+jitter',RF_perpos{animal,pos},...
                'MarkerEdgeColor','none','MarkerFaceColor',color_plot(animal,:),'MarkerFaceAlpha',0.2);
            plot([RespROIs.V1az(animal,pos)-2.5 RespROIs.V1az(animal,pos)+2.5],[mean_RF(animal,pos) mean_RF(animal,pos)],...
                '-','linewidth',3,'color',color_plot(animal,:))
            % -- same for the final figure
            figure(indivAnimalFig)
            subplot(1,n_animals,animal);hold on
            scatter(RespROIs.V1az(animal,pos)*ones(length(nBoutons),1)+jitter',RF_perpos{animal,pos},...
                'MarkerEdgeColor','none','MarkerFaceColor',color_plot(animal,:),'MarkerFaceAlpha',0.2);
            plot([RespROIs.V1az(animal,pos)-2.5 RespROIs.V1az(animal,pos)+2.5],[mean_RF(animal,pos) mean_RF(animal,pos)],...
                '-','linewidth',3,'color',[0 0 0])
        end
    end
end

figure(indivAnimalFig)
for animal = 1:n_animals
    forFit = mean_RF(animal,:);
    V1az_all =RespROIs.V1az(animal,:);
    izGoodPos = isGoodPos(animal,:);
    forFit = forFit(izGoodPos)';
    V1az_all = V1az_all(izGoodPos);
    X = [ones(length(V1az_all),1) V1az_all'];
    subplot(1,n_animals,animal)
    plot([-20 100],[-20 100],'k:')
    [b,bint,~,~,stats] = regress(forFit,X);
    yCalc1 = X*b;
    plot(V1az_all,yCalc1,'k--','linewidth',2)
    title([RespROIs.info.animalID{animal} newline ...
        'y=' num2str(b(2),3) 'x+' num2str(b(1),3) newline ...
        'R2=' num2str(stats(1),3) '|p=' num2str(stats(3),3)])
    if animal == 1
        ylabel(['Axon ' ROItype ' RF']);
    end
    xlabel('V1 RF'); xlim([-30 110])
    ylim([-30 110])
    axis square
end
set(gcf,'units','normalized','position',[.05 .4 .9 .3])
saveas(gcf,[saveFolder filesep dataset '_' ROItype '_BoutonsRFvsV1RF_IndivMouse.tif']);
savefig([saveFolder filesep dataset '_' ROItype '_BoutonsRFvsV1RF_IndivMouse']);
close

% --- fit the mean RF per session
figure(ScatterFig)
plot([-20 100],[-20 100],'k:')
xlabel('V1 RF'); xlim([-30 110])
ylabel(['Axon ' ROItype ' RF']); ylim([-30 110])

meanRF_all = reshape(mean_RF,size(mean_RF,1)*size(mean_RF,2),1);
isGoodPos_reshape = reshape(isGoodPos,size(mean_RF,1)*size(mean_RF,2),1);
meanRF_all = meanRF_all(isGoodPos_reshape);
V1az_all = reshape(RespROIs.V1az,size(RespROIs.V1az,1)*size(RespROIs.V1az,2),1);
V1az_all = V1az_all(isGoodPos_reshape);

X = [ones(length(V1az_all),1) V1az_all];
[b1,~,~,~,stats1] = regress(meanRF_all,X);
yCalc1 = X*b1;

plot(V1az_all,yCalc1,'k--','linewidth',2)
axis equal

equation1 = sprintf('y=%.2fx+%.1f',b1(2),b1(1));

% --- fit all the RF
V1RF_perpos = cell(n_animals,max(n_pos));
for animal = 1:n_animals
    for pos = 1:n_pos(animal)
        V1RF_perpos{animal,pos} = repmat(RespROIs.V1az(animal,pos),1,size(RF_perpos{animal,pos},2));
    end
end

V1RF_permice = cell(n_animals,1);
RF_permice = cell(n_animals,1);
for animal = 1:n_animals
    RF_permice{animal} = cat(2,RF_perpos{animal,:});
    V1RF_permice{animal} = cat(2,V1RF_perpos{animal,:});
end
RF_all = cat(2,RF_permice{:});
V1RF_all = cat(2,V1RF_permice{:});

X = [ones(length(V1RF_all),1) V1RF_all'];
[b2,~,~,~,stats2] = regress(RF_all',X);
yCalc2 = X*b2;
plot(V1RF_all,yCalc2,'k--','linewidth',2)
equation2 = sprintf('y=%.2fx+%.1f',b2(2),b2(1));

title('RF azimtuh center (deg, head coord.)')

subplot(1,4,4); hold on
for animal = 1:n_animals
    if sum(isGoodPos(animal,:))>0
        scatter(0,1-0.05*(animal-1),'MarkerEdgeColor','none','MarkerFaceColor',color_plot(animal,:),'MarkerFaceAlpha',0.2);
        text(0.1,1-0.05*(animal-1),num2str(RespROIs.info.animalID{animal}))
    end
end
text(0.1,1-0.05*n_animals-0.05,{'session mean',equation1,['R2=' num2str(stats1(1),2)],['p=' num2str(stats1(3))]})
text(0.1,1-0.05*n_animals-0.25,{'all boutons' equation2,['R2=' num2str(stats2(1),2)],['p=' num2str(stats2(3))]})

temp = reshape(nRespBoutons,n_animals*max(n_pos),1);
temp = temp(isGoodPos_reshape);
nRespBoutons_all = cumsum(temp,'omitnan');
nRespBoutons_all = nRespBoutons_all(end);
text(0.1,1-0.05*n_animals-0.45,{['n=' num2str(nRespBoutons_all) ' ' ROItype ','],[num2str(length(temp)) ' sessions']})

ylim([0 1]);xlim([0 1]);axis off

saveas(gcf,[saveFolder filesep dataset '_' ROItype '_BoutonsRFvsV1RF.tif']);
saveas(gcf,[saveFolder filesep dataset '_' ROItype '_BoutonsRFvsV1RF.svg']);
savefig([saveFolder filesep dataset '_' ROItype '_BoutonsRFvsV1RF']);
close

fprintf(1,'done \n')
end


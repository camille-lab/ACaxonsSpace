function Decoder_Output = BayesianDecoder_batch_SingleRow(RespROIs,MainFolder,az_vector)

%% a few parameters
corrThresh = 0.3;         % to determine which ROIs belong to the same axon
training_fraction = 0.8;
animalList = RespROIs.info.animalID;

%% a few definitions
n_animals = size(RespROIs.info.animalID,2);
ispos = ~cellfun(@isempty,RespROIs.info.data_details);
n_pos = sum(ispos,2);

%% mega virtual mouse, aka Mickey Mouse
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

nAxons = [10 25 50 100];
nDraws = 100; % 100 in Mazo et al., 2024
counter = 0;
for i = 1:length(nAxons)
    for ii = 1:nDraws
        counter=counter+1;
        y = datasample(data_all,nAxons(i),'Replace',false);
        Posterior_MickeyMouse{i,ii} = BayesianDecoder_v4_SingleRow(y,training_fraction,NaN,0,[],RespROIs.info,'MickeyMouse',az_vector,[]);
    end
end

%% (decoded-actual) position as a function of number of axons
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

%% save
Decoder_Output.MickeyMouse = Posterior_MickeyMouse;
Decoder_Output.delta = delta;
end
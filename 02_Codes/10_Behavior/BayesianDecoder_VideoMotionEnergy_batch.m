% function BayesianDecoder_VideoMotionEnergy_batch(resp_window,animals,eyedata,dataset,SaveFolder)
SaveFolder = 'D:\AVspace_final\03_Plots\10_BehaviorVideos\Face'; 

resp_window = [1.5 3];
dataset = 'CBA'; % 'BL6','CBA' -- Bl6 data is not presented in Mazo et al., 2024

%% Define some variables
training_fraction = 0.8;
nRep = 20;
az_vector = [-20:10:100];
nAzPos = length(az_vector);
nElePos = 3;
nTrialTypes = 39;
nTrials =nTrialTypes*nRep;
nPCs = 500;

%% mice to use
if strcmpi(dataset,'bl6')
    animals = {'CMad97','CMad98'};
elseif  strcmpi(dataset,'cba')
    animals = {'CMad103','CMad104','CMad106','CMad109','CMad110','CMad111'};
end

SaveFolder = [SaveFolder filesep filesep dataset]; 

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

%%
n_animals = length(animals);
collapsed = NaN(length(animals),size(faceData.videoME_sorted,2),13);
collapsed_shuffle = NaN(length(animals),size(faceData.videoME_sorted,2),13);

decodingError = NaN(length(animals),size(faceData.videoME_sorted,2),5,39,4);
decodingError_sh = NaN(length(animals),size(faceData.videoME_sorted,2),5,39,4);
for mice = 1:n_animals
    for pos = 1:size(faceData.videoME_sorted,2)
        if ~isempty(faceData.videoME_sorted{mice,pos})
            data = faceData.videoME_sorted{mice,pos}(:,:,:,1:nPCs);
            
            nFrames = size(data,3);
            fps = nFrames/5; % a trial is 5s
            x_vect = 0:1/fps:5; x_vect=x_vect(1:nFrames);

            tAna = false(size(x_vect));
            tAna(x_vect>=resp_window(1) & x_vect<=resp_window(2))=true;
            tAna = tAna(1:size(data,3));
            
            %% Built the decoder and the response matrix
            for iii = 1:nRuns
                respDist = squeeze(nanmean(data(:,:,tAna,:),3));
                
                sdMatrix = squeeze(std(respDist(:,trainingTrials(iii,:),:),[],2,'omitnan'));
                RespMatrix = squeeze(nanmean(respDist(:,trainingTrials(iii,:),:),2));
                
                meanResp = respDist(:,decodedTrials(iii,:),:);
                meanResp = reshape(meanResp,nDecodedTrials,nPCs);
                
                % - shuffle -
                respDist_shuffle = NaN(nTrialTypes,nTrainingTrials,nPCs);
                for i = 1:nTrialTypes
                    for j = 1:nTrainingTrials
                        respDist_shuffle(i,j,:) = respDist(i,j,randperm(nPCs));
                    end
                end
                sdMatrix_shuffled = squeeze(std(respDist_shuffle,[],2,'omitnan'));
                RespMatrix_shuffled = squeeze(nanmean(respDist_shuffle,2));
                shuffle_meanResp = meanResp;
                % ------------
                
                for ii = 1:nDecodedTrials
                    for j = 1:nTrialTypes
                        for i = 1:nPCs
                            Posterior_perROI(i) = -log(sdMatrix(j,i)*sqrt(2*pi)) - (meanResp(ii,i)-RespMatrix(j,i))^2./(2*sdMatrix(j,i)^2);
                            Posterior_perROI_shuffled(i) = -log(sdMatrix_shuffled(j,i)*sqrt(2*pi)) - (shuffle_meanResp(ii,i)-RespMatrix_shuffled(j,i))^2./(2*sdMatrix_shuffled(j,i)^2);
                        end
                        Posterior_perStimType(j) = nansum(Posterior_perROI);
                        Posterior_perStimType_shuffled(j) = nansum(Posterior_perROI_shuffled);
                    end
                    PosteriorMatrix{iii}(:,ii) = Posterior_perStimType;
                    PosteriorMatrix_shuffled{iii}(:,ii) = Posterior_perStimType_shuffled;
                end
            end
            
            %%
            for iii = 1:nRuns
                PosteriorMatrix_perStimType{iii} = reshape(PosteriorMatrix{iii},...
                    nTrialTypes,nTrialTypes,nDecodedTrials_perTrialtype);
                PosteriorMatrix_perStimType_shuffled{iii} = reshape(PosteriorMatrix_shuffled{iii},...
                    nTrialTypes,nTrialTypes,nDecodedTrials_perTrialtype);
            end
            
            c = cat(3,PosteriorMatrix_perStimType{:}); % concatenate across runs
            
            %% Decode
            % -- Spatial distribution of decoding
            [X,Y]= meshgrid([-20:10:100],[-20:20:20]);
            Xlin = reshape(X',39,1); Ylin = reshape(Y',39,1);
        
            deltaDecoded2Actual = NaN(nRuns,nTrialTypes,nDecodedTrials_perTrialtype);
            deltaDecoded2Actual_sh = NaN(nRuns,nTrialTypes,nDecodedTrials_perTrialtype);
            for iii = 1:nRuns
                for j = 1:nTrialTypes
                    for i = 1:nDecodedTrials_perTrialtype
                        [~, I] = max(PosteriorMatrix_perStimType{iii}(:,j,i));
                        [~, I2] = max(PosteriorMatrix_perStimType_shuffled{iii}(:,j,i));
%                         if  I==j || I==j+13 || I==j+26 || I==j-13 || I==j-26 || I==j+1 || I==j-1 || I==j+14 || I==j+12 || I==j+27 || I==j+25 || I==j-12 || I==j-14 || I==j-25 || I==j-27
%                             DecodedAz(iii,j,i) = 1;
%                         else
%                             DecodedAz(iii,j,i) = 0;
%                         end
%                         if I2==j || I2==j+13 || I2==j+26 || I2==j-13 || I2==j-26 || I2==j+1 || I2==j-1 || I2==j+14 || I2==j+12 || I2==j+27 || I2==j+25 || I2==j-12 || I2==j-14 || I2==j-25 || I2==j-27
%                             DecodedAz_shuffle(iii,j,i) = 1;
%                         else
%                             DecodedAz_shuffle(iii,j,i) = 0;
%                         end
                        
                        deltaDecoded2Actual(iii,j,i) = norm([Xlin(I),Ylin(I)]-[Xlin(j),Ylin(j)]);
                        deltaDecoded2Actual_sh(iii,j,i) = norm([Xlin(I2),Ylin(I2)]-[Xlin(j),Ylin(j)]);
%                      
                    end
                    
                end
            end
            
decodingError(mice,pos,:,:,:) = deltaDecoded2Actual;
decodingError_sh(mice,pos,:,:,:) = deltaDecoded2Actual_sh;
        end
    end % end loop through recording position
end % end loop through mice

%%
% % - calculate chance
for i = 1:39
    for ii = 1:39
        chance(i,ii) =  norm([Xlin(i),Ylin(i)]-[Xlin(ii),Ylin(ii)]);
    end
end
chance_avg = mean(chance,'all');

for i = 1:n_animals
mean_decodingError(i) = mean(decodingError(i,:,:,:,:),'all','omitnan');% squeeze(mean(mean(deltaDecoded2Actual,1),3));
mean_decodingError_sh(i) = mean(decodingError_sh(i,:,:,:,:),'all','omitnan');
end

figure; hold on
plot([.5 2.5],[chance_avg chance_avg],'k:')
plot(repmat([1 2]',1,n_animals),[mean_decodingError;mean_decodingError_sh],'k')
scatter(ones(n_animals,1),mean_decodingError,'r')
scatter(2*ones(n_animals,1),mean_decodingError_sh,'k')
plot([.9 1.9;1.1 2.1],[mean(mean_decodingError) mean(mean_decodingError_sh);mean(mean_decodingError) mean(mean_decodingError_sh)],...
    'k','linewidth',2)
xlim([.5 2.5]) ; xticks([1 2])
ylabel('decoding error (degree)'); ylim([40 60]); 
[~,p]=ttest(mean_decodingError,mean_decodingError_sh);
title(['p ttest, p=' num2str(p)])

% % - save
savefig([SaveFolder filesep dataset 'DecodingError'])
saveas(gcf,[SaveFolder filesep dataset 'DecodingError.png'])
saveas(gcf,[SaveFolder filesep dataset 'DecodingError.svg'])

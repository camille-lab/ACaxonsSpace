% written Sep 2021
% use directly R. otherwise same as v1 but for the new input (data_Sorted_v4)

function axonRois = getAxonID_CM_v2(R,corrThresh) %correlation_data
sameAxon = R>corrThresh;
sessionN = 1;
nROIs = size(R,1);

percentBiggerThresPerSession= mean(mean(sameAxon,2)*length(sameAxon)-1)/length(sameAxon);
twobincountsCorrelations=histcounts(R,[0 corrThresh 0.99]);
bincountsCorrelations=histcounts(R,0:0.05:0.999);

correlationsPercetBigThres=bincountsCorrelations(1,2)/ sum(bincountsCorrelations(1,:),2);
% correlationsPercetBigThres is the percentage of pairs that have highger
% correlation than threshold. was 0.2% in the paper.

% -- identify correlated ROIs
for ROI = 1:nROIs
    axonsListPerSession{1,ROI} = find(sameAxon(ROI,:));
end

axonsListPerSessionMat{sessionN} = inf(ROI);
for ROI = 1:nROIs
    axonsListPerSessionMat{sessionN}(ROI,1:size(axonsListPerSession{sessionN,ROI},2)) = axonsListPerSession{sessionN,ROI};
end

% -- matrix of unique correlated ROIs
axonsTEMP{sessionN} = unique(axonsListPerSessionMat{sessionN},'rows');

% ---
axonsID_f = {};
% for sessionN=1:length(sessions)
for ROI = 1:size(axonsTEMP{sessionN},1)
    axonsID_f{sessionN,ROI}=axonsTEMP{sessionN}(ROI,~isinf(axonsTEMP{sessionN}(ROI,:)) );
end
% end


%% ---- me scripting here: to seed the ROIs
clear axon
temp = axonsID_f;
soloROI = cell(1,length(temp));
for i =length(temp):-1:1 % loop through axons and remove the soloistes from 'temp'
    if size(temp{i})==1
        soloROI{i}=temp{i};
        temp{i}=[];
    end
end
clubbedROIs=temp(~cellfun(@isempty,temp));
soloROI = soloROI(~cellfun(@isempty,soloROI));

if isempty(clubbedROIs)      % no ROIs pass correlation threshold
    axonRois=soloROI;  
else                         % seed the actual axons
    c=1;
    while length(clubbedROIs)>0
        clear temp
        rdSeed= randperm(length(clubbedROIs),1);
        rdPair = randperm(length(clubbedROIs{rdSeed}),2);
        for i =1:2
            temp{i} = find(R(clubbedROIs{rdSeed}(rdPair(i)),:)>corrThresh); % find the correlated ROIs to either ROIs of the pair
        end
        temp2=cat(2,temp{:});
        axon{c} = unique(temp2);
        clubbedROIs{rdSeed}=[];                        % remove this axon
        clubbedROIs = clubbedROIs(~cellfun(@isempty,clubbedROIs));   % remove the field to resize 'test'
        R(:,axon{c})=NaN;                       % remove the clubbed ROIs from correlation matrix
        c=c+1;
    end
    
    axon=axon(~cellfun(@isempty,axon));
    axonRois = [soloROI,axon];
end

% test3 = cat(2,axonRois{:});
% temp=unique(test3);

% indexToDupes = find(not(ismember(1:numel(test3),i))) % check that no ROI is represented twice


% R2 = corrcoef(data(test3,:)');
% figure; imagesc(R2)
% test2=cat(2,axon{:});
% for i = 1:length(test2)
%     ROIsID = test2{i};
%     if any(ROIsID
% end

% Transform the data back into trials
% Sort the data according to trial types

function data_sorted = SortData_CM_v4_SPKalone_Batch(fTracesArray,animalID,position,path,saveDir_Var)
nTrialsPer2Pfile = 4;

%% Preliminary calculations, some checks

% - load trial timeline; if animalID does not match ID of the stim ops
% file, it will fail
timeline_dir = dir([path.folder filesep path.name filesep '*' char(position) filesep 'trial_timeline' '*.csv']);
trialtimeline = load([timeline_dir.folder filesep timeline_dir.name]);
disp(timeline_dir.name)

saveDir_variable = [saveDir_Var,filesep,char(animalID),filesep];

DataID = [char(animalID) '_' char(position)];
disp(DataID)
if  ~exist(saveDir_variable,'dir')
    mkdir(saveDir_variable)
end

extracting_ops.software = 2;
extracting_ops.SplitROIs = 0;
extracting_ops.plotFlag = 0;
extracting_ops.PureTone = 0;
extracting_ops.zoom = 8;
axonID_data = 1;

planesArray = 1:length(fTracesArray);
disp(['n planes=' num2str(planesArray(end))]);

dt = datestr(now,'yyyymmdd');
oldfolder = cd(saveDir_variable);
diaryName = ['SortDatav4_Log_' DataID '_' dt];
diary off; diary(diaryName)
cd(oldfolder)

nFTraces = length(fTracesArray);
nFiles = length(fTracesArray{1}.data);
nTrials = nFiles*nTrialsPer2Pfile;

if mod(nTrials,39) ~= 0 % i.e not the same number of stim for each position
    nTrials = floor(nTrials/39)*39;
    disp('trials types were not balanced');
    unbalancedTrials = 1;
else
    unbalancedTrials = 0;
end

data_sorted = fTracesArray{1,1};
try
    data_sorted.header = fTracesArray{1,1}.H;
catch
    data_sorted.header =  fTracesArray{1,1}.header;
end

for i=1:nFTraces
    ROI_planeNumber_temp{i} = planesArray(i)*ones(1,size(fTracesArray{1,i}.data{1,i},1));
    nROIs_perPlane(i) = size(fTracesArray{1,i}.data{1,i},1);
end
ROI_planeNumber = cat(2,ROI_planeNumber_temp{:});
nROIs = sum(nROIs_perPlane);

fprintf(['nROIs per plane = ' num2str(nROIs_perPlane)])
disp([' | Total = ' num2str(nROIs)])

minmalTrialLength = min(fTracesArray{1,1}.header.numberOfFrames);
if minmalTrialLength < 6*8
    fprintf(2,'Trials were chopped\n')
    chopped = 1;
    keyboard
else
    chopped = 0;
end
minmalTrialLength = minmalTrialLength - mod(minmalTrialLength,nTrialsPer2Pfile);
minmalTrialLength = double(minmalTrialLength);

%% Concatenate all the ROIs in single matrix
% this concatenate all the ROIs from different planes already, trial-by-trial
if nFTraces == 1 % if single plane
    for i = 1:nTrials
        data_sorted.Fval{:,i} = data_sorted.data{1,i}(:,1:minmalTrialLength);
        data_sorted.spikes{:,i} = data_sorted.spikes{1,i}(:,1:minmalTrialLength);
    end
end
for i = 1:nFiles
    for k = 2:nFTraces
        data_sorted.data{:,i} = [data_sorted.data{1,i}(:,1:minmalTrialLength); fTracesArray{1,k}.data{1,i}(:,1:minmalTrialLength)];
        data_sorted.spikes{:,i} = [data_sorted.spikes{1,i}(:,1:minmalTrialLength); fTracesArray{1,k}.spikes{1,i}(:,1:minmalTrialLength)];
    end
end

% stack trials in the 3rd dimension (not sorted)
for k = 1:nFiles
    Fvalues(:,:,k) = cat(1,data_sorted.data{:,k});
    spks(:,:,k) = cat(1,data_sorted.spikes{:,k});
end

Fvalues = reshape(Fvalues,nROIs,minmalTrialLength/nTrialsPer2Pfile,nTrials);
spks = reshape(spks,nROIs,minmalTrialLength/nTrialsPer2Pfile,nTrials);

dtCa = 1/fTracesArray{1}.header.frameFrequency;
disp(['sampling frequency (per plane) is: ', num2str(1/dtCa,3), 'Hz'])

%% Calculate the pearson's correlation
% -- F values
concat=cat(2,data_sorted.data{:});
R_fval_all = corrcoef(concat');

% -- spikes
concat2=cat(2,data_sorted.spikes{:});
R_spikes_all = corrcoef(concat2');

%% Sort the trials according to SPK stim ID
nTrialTypes_Spk = 39;

%  -- Get the SPK ID for each trial
SPK = trialtimeline(:,4:6);

for i = 1:3
    spkID{:,i} = find(SPK(:,i)~=-1);
end
for i = 1:3
    spkID2{:,i} = SPK(spkID{:,i},i);
end

for i = 1:3
    SPK_timing{i} = spkID{1,i}(1:2:end,1);
end
SPK_ID{1} = spkID2{1,1}(1:2:end,1)+15;     % Top board controls columns 6-9
SPK_ID{2} = spkID2{1,2}(1:2:end,1)+15+12;  % Middle board controls columns 10-13
SPK_ID{3} = spkID2{1,3}(1:2:end,1);

Timing_all = cat(1,SPK_timing{:});
[~,timing_sorted] = sort(Timing_all);
SpkID_all = cat(1,SPK_ID{:});
SpkIDall_sorted = SpkID_all(timing_sorted);
SpkIDall_sorted = SpkIDall_sorted +1;

tracesDist_SPK = NaN(nROIs,minmalTrialLength/nTrialsPer2Pfile,nTrialTypes_Spk,nTrials/nTrialTypes_Spk);
spikesDist_SPK = NaN(nROIs,minmalTrialLength/nTrialsPer2Pfile,nTrialTypes_Spk,nTrials/nTrialTypes_Spk);

matrix = [1 :13;14:26;27:39];
linearized = reshape(matrix,3*13,1); % or matrix(:)
for ii = 1:nTrialTypes_Spk
    trial_index(ii,:) = find(SpkIDall_sorted==ii);
    tracesDist_SPK(:,:,linearized(ii),:) = Fvalues(:,:,find(SpkIDall_sorted==ii));
    spikesDist_SPK(:,:,linearized(ii),:) = spks(:,:,find(SpkIDall_sorted==ii));
end

%% Save s2p info
for k = 1:nFTraces
    temp{k} = fTracesArray{k}.s2p;
end
data_sorted.s2p = cat(1,temp{:});

%% Save
fprintf('saving...')

data_sorted = rmfield(data_sorted,'data');
data_sorted = rmfield(data_sorted,'spikes');

data_sorted.Info.dtCa               = dtCa;
data_sorted.Info.nROIs             = nROIs;
data_sorted.Info.nROIs_perPlane    = nROIs_perPlane;
data_sorted.Info.ROI_planeNumber   = ROI_planeNumber;
data_sorted.Info.nTrialTypes_Spk   = nTrialTypes_Spk;
data_sorted.Info.nTrials           = nTrials;
data_sorted.Info.ActualnRep        = nTrials/nTrialTypes_Spk;
data_sorted.Info.saveName          = DataID;
data_sorted.Info.trialtimelineUsed = timeline_dir;
data_sorted.Info.extracting_ops    = extracting_ops;
data_sorted.Info.dateOfExtraction  = dt;
data_sorted.Info.trial_index       = trial_index;

data_sorted.pearsonR.dFF     = R_fval_all;
data_sorted.pearsonR.spikes  = R_spikes_all;
    
data_sorted.SPKwn.tracesDist = tracesDist_SPK;
data_sorted.SPKwn.spikesDist = spikesDist_SPK;
data_sorted.SPKwn.trials_SPK = SpkIDall_sorted;

saveName = ['dataSorted_' DataID];
save([saveDir_variable, filesep, saveName],'data_sorted', '-v7.3')

disp('done')
diary off
end
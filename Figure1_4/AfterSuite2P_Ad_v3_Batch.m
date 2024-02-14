% Chop the concatenated suite2P data back into 2P files (could be several
% trials)
% Use 'iscell'
% Stack data across planes and ignores flyback

function Ftraces_all = AfterSuite2P_Ad_v3_Batch(animal,pos,path,trial_length,varargin)
%Script is gonna import the different Fall values from suite2p folder

%% Preambule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isFlyBack = 1;

foldername = [path.folder filesep path.name filesep 'Axons'];

DataID = [char(animal) '_' char(pos)];
n_planes = 4;
frameFrequency =  30.4812/(n_planes+isFlyBack);
if nargin == 5
    nTrials_expected = cell2mat(varargin);
else
    nTrials_expected = 195;
end
expectedNumberOfFrames = 38; % only if need to chopp over- or under-sized files AFTER suite2p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find the suite2p folder
iscell =[];
s2p_folder = dir([foldername filesep '\suite2p']); 

match = {}; counter = 1;
for i = 1:length(s2p_folder)
    if ~isempty(regexp(s2p_folder(i).name,'^plane'))
        match{end+1} = s2p_folder(i).name;
        index(counter) = i;
        counter = counter +1;
    end
end
if length(index)>1
    index = index(1:n_planes);
    plane = 'multiples';
else
    plane = 'single';
end
dt = datestr(now,'yyyymmdd');

filetype = 'tif'; % type of files to be processed
saveFolder = [foldername '\analysis'];
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

oldfolder = cd(saveFolder);
diaryName = ['FtracesLog_' dt];
diary off; diary(diaryName)
cd(oldfolder)

saveName = ['Ftraces_' DataID '_' dt];

disp(saveFolder); disp(saveName);

%% Get the number of frames per 2P files to prep for chopping concatenated traces.
% NofFrames is the same for all planes. Use first plane to retrieve info
path = [s2p_folder(index(1)).folder filesep s2p_folder(index(1)).name];
load([path '\Fall.mat'])
numberOfFrames = ops.frames_per_file;

nTrials = length(numberOfFrames);

if nTrials ~= nTrials_expected
    fprintf(2,'wrong number of files')
    
    files = subdir(fullfile(foldername,['*',filetype]));   % list of filenames (will NOT search all subdirectories)
    [H, Data] = opentif(files(1).name);
    
    nTrials = length(files);
    disp([', found ',num2str(nTrials), ' tif files'])
    
    numberOfFrames = zeros(nTrials,1);
    FlyBackFrame = H.SI.hFastZ.discardFlybackFrames; % If flyback was automatically discarded, will set it to 1
    f =  waitbar(0);
    for i = 1:nTrials
        info = imfinfo(files(i).name);
        numFrames(i) = numel(info)/(size(Data,3)*(size(Data,5)+FlyBackFrame));
        waitbar(i/nTrials,f,[num2str(i/nTrials*100,2),'%'])
    end
    close(f)
    j = 1; shortTrialFrames = 0; short_trial =0;
    for i = 1:nTrials
        if numFrames(i) > 2*frameFrequency*n_planes*trial_length
            numberOfFrames(j+1) = numFrames(j)-floor(numFrames(i)/2);
            numberOfFrames(j) = floor(numFrames(i)/2);
            j = j+1;
            short_trial = 0;
        elseif numFrames(i) < frameFrequency*n_planes*trial_length-3

            numberOfFrames(j)= NaN;
            shortTrialFrames = shortTrialFrames+numFrames(i);
            short_trial = 1;
        else
            if ~short_trial
            numberOfFrames(j) = numFrames(i);
            else
                numberOfFrames(j) = numFrames(i)+shortTrialFrames;
                shortTrialFrames = 0;
            end
            short_trial = 0;
        end
        j = j+1;
    end
    numberOfFrames(isnan(numberOfFrames))=[];
    nTrials = length(numberOfFrames);
    if nTrials ~= nTrials_expected
        keyboard
    end
end

disp(['nTrials: ' num2str(nTrials)])

%% Actually chop the concatenated trials into stacked 2P files
for ii = 1:length(index)
    path = [s2p_folder(index(ii)).folder filesep s2p_folder(index(ii)).name];
    
    fprintf(s2p_folder(index(ii)).name)
    load([path '\Fall.mat'])
    
    counter = 1;
    is_cell = logical(iscell(:,1));
    for i = 1:nTrials
        Ftraces_all{ii}.data{i} = F(is_cell,counter:counter+numberOfFrames(i)-1);
        Ftraces_all{ii}.spikes{i} = spks(is_cell,counter:counter+numberOfFrames(i)-1);
        counter = counter + numberOfFrames(i);
    end

    fprintf(' | nROIs: %d \n',size(Ftraces_all{ii}.data{1},1))
    
    %% get header info from intact/splitted file
    Ftraces_all{ii}.header.frameFrequency = frameFrequency;
    Ftraces_all{ii}.header.DataID = DataID;
    Ftraces_all{ii}.header.mainDir = foldername;    
    
    if isfield(stat{1,1},'iplane')
        for i = 1:size(F,1)
            Ftraces_all{ii}.planeNumber(i) = stat{1,i}.iplane+1;
        end
    end
    
    Ftraces_all{ii}.Fneu = Fneu;
    
    Ftraces_all{ii}.s2p.iscell   = iscell(:,1);
    Ftraces_all{ii}.s2p.meanImg  = ops.meanImg;
    Ftraces_all{ii}.s2p.meanImgE = ops.meanImgE;
    Ftraces_all{ii}.s2p.yOffset  = ops.yoff;
    Ftraces_all{ii}.s2p.xOffset  = ops.xoff;
    Ftraces_all{ii}.s2p.spikes   = spks(logical(iscell(:,1)),:);
    for i = 1:length(is_cell)
        if is_cell(i)
            Ftraces_all{ii}.s2p.stat{1,i}  = stat{1,i} ;
        end
    end
    
    Ftraces_all{ii}.header.numberOfFrames = numberOfFrames;
    
    Ftraces_all{ii}.Info.saveName = DataID;
    Ftraces_all{ii}.Info.mainDir = foldername;
    Ftraces_all{ii}.Info.fileList = ops.filelist;
    Ftraces_all{ii}.Info.plane = plane;
end

%% save
fileName = [saveFolder,'\',saveName,'.mat'];
save(fileName,'Ftraces_all')
diary off
% Run AfterSuite2p and SortData_CM_v4.
% Works for both all data in principle (2-80 kHz Bl6 or CBA, 2-80 kHz TEA)
% trial timeline file in 'mainFolder' 'mouseID' '*posXX' filesep 'trial_timeline' '*.csv'
%%
mainFolder = 'D:\2P\'; % 'oragnized like 'D:\2P\' 'mouseID' '*posXX'
saveDir_Var = ''; % this was 'C:\Users\camil\Data\AdAxons_analysis'

dataset = 'CBA'; % 'Bl6','CBA'

trial_length = 5; % (s)
%%
switch dataset
    case 'CBA'
        % CMad103, pos8-13 curated and sorted
        % CMad104, pos1-9  curated and sorted, all done
        % CMad106, pos1-5  curated and sorted, all done
        % CMad107, pos1-7  curated and sorted, all done
        % CMad109, pos1-9  curated and sorted, all done
        % CMad110, pos1-10 curated and sorted, all done
        % CMad111, pos1-10 curated and sorted, all done
        
        animalID = {'CMad103','CMad104','CMad106','CMad107','CMad109','CMad110','CMad111'};
        positions = {
            'pos02','pos03','pos05','pos06','pos08','pos09','pos10','pos11','pos12','pos13';...      % CMad103
            'pos01','pos02','pos03','','','','','','','';...                                         % CMad104
            'pos01','pos02','pos03','pos04','pos05','','','','','';...                               % CMad106
            'pos01','pos02','pos03','pos04','pos05','pos06','','','','';...                          % CMad107
            'pos01','pos02','pos03','pos04','pos05','pos06','pos07','pos08','pos09','';...           % CMad109
            'pos01','pos02','pos03','pos04','pos05','pos06','pos07','pos08','pos09','pos10';...      % CMad110
            'pos01','pos02','pos03','pos04','pos05','pos06','pos07','pos08','pos09','pos10'};        % CMad111

    case 'Bl6'
        %         animalID = {'CMad97','CMad98','CMad99'};
        %         positions = {
        %             'pos01','pos02','pos03','pos04','pos05','pos06','pos07','pos08';...  % CMad97
        %             'pos01','pos02','pos03','pos04','pos05','pos06','pos07','pos08';...  % CMad98
        %             'pos01','pos02','pos03','pos04','pos05','pos06','pos07','pos08'};    % CMad99
        animalID = {'CMad99'};
        positions = {
            'pos01','pos02','pos03','pos04'};    % CMad99
end

%%
for mouse = 1:length(animalID)
    fprintf(2,[animalID{mouse} newline])
    actual_pos=~cellfun(@isempty,positions);
    for pos = 1:size(positions,2)
        if actual_pos(mouse,pos)
            path = dir([mainFolder char(animalID(mouse)) filesep '*' char(positions{mouse,pos})]);
            disp([newline path.name])
            
            Ftraces_all = AfterSuite2P_Ad_v3_Batch(animalID{mouse},positions(mouse,pos),path,trial_length);
            
            data_sorted = SortData_CM_v4_SPKalone_Batch(Ftraces_all,animalID{mouse},positions(mouse,pos),path,saveDir_Var);
        end
    end
    fprintf(newline)
end

disp('all done')
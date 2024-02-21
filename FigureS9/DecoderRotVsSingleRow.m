function DecoderRotVsSingleRow
nAxons = [10 25,50,100]; 
nDraws = 100;
nAzPos = 12;
myfit = fittype('a + b*log(x)',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','b'});

az_vector = [-10:10:100];

%% manage the input

    MD_SingleRow = load('C:\Users\camil\AudSpace\Data\Decoder\MD_SingleRow_Deg_2mice.mat');
    RotSpk = load('C:\Users\camil\AudSpace\Data\Decoder\RotSpk_deg_2mice.mat');
    
%% a few definitions
DecoderFig = figure;
SaveFolder = ['C:\Users\camil\AudSpace\Plots\Decoder' filesep 'MDvsROT'];
if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end

%% decoding error (i.e, distance (decoder-actual) position, in 2D, in degrees
decodingError_MD_SingleRow = MD_SingleRow.Decoder_Output_SingleRowMD.delta(1:length(nAxons),1:nDraws);
decodingError_RotSpk = RotSpk.Decoder_Output_RotSPK.delta(1:length(nAxons),1:nDraws);

% -- plot
Acc_avg = mean(decodingError_MD_SingleRow,2);
P = prctile(decodingError_MD_SingleRow,[5 95],2);
figure(DecoderFig); hold on
fill([nAxons flip(nAxons)],[P(:,1)' flip(P(:,2))'] ,[0 0 1],'edgecolor','none','facealpha',0.2)
plot(nAxons,Acc_avg,'bo-');

Acc_avg = mean(decodingError_RotSpk,2);
P = prctile(decodingError_RotSpk,[5 95],2);
fill([nAxons flip(nAxons)],[P(:,1)' flip(P(:,2))'] ,[1 0 0],'edgecolor','none','facealpha',0.2)
plot(nAxons, Acc_avg,'ro-');

temp = NaN(12,12);
for i=1:12
    for j = 1:12
    temp(i,j) = norm(az_vector(i)-az_vector(j));
    end
end
chance = mean(temp,'all');

plot([0 nAxons(end)],[chance chance],'k:')
yticks([0 25 50 75 100]); ylim([0 100]); ylabel('decoding error (degrees)')
xticks(nAxons);xlabel('number of axons')

% -- save
sgtitle('ROTvsMD')
set(gcf,'units','normalized','position',[.2 .1 .7 .7])
savefig(DecoderFig,[SaveFolder filesep 'ErrorVsNaxons_ROTvsMD'])
saveas(gcf,[SaveFolder filesep 'ErrorVsNaxons_ROTvsMD.tif']);
saveas(gcf,[SaveFolder filesep 'ErrorVsNaxons_ROTvsMD.svg']);    
    

end
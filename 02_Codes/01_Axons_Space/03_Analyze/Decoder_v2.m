% dataset = '2-20kHz' '2-80kHz'
% 'All_Error.mat' and 'All_Error.mat' contain the 'Decoder_Output' variables created running 'Batch_dataResponsive_v2.m'
% change the paths in line 12 and 15
function Decoder_v2(varargin)
MainFolder = 'D:\AVspace_final';
nAxons = [10 25,50,100,200,300];
nDraws = 100;

%% manage the input
[indx,tf] = listdlg('PromptString','choose the dataset','ListString',{'2-20kHz','2-80kHz'});
if indx == 1
    load([MainFolder filesep '01_Data\Decoder\All_Error.mat'])
    dataset = '2-20kHz';
else
    load([MainFolder filesep '01_Data\Decoder\All_Error_80kHz.mat')
    dataset = '2-80kHz';
end
disp(dataset)

%% save folder
SaveFolder = ['C:\Users\camil\AudSpace\Plots\Decoder' filesep dataset];
if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end

%% Data
decodingError_ACaud = mean(Decoder_Output_AC.delta,2);
decodingErrorCI_ACaud = prctile(Decoder_Output_AC.delta,[5 95],2);

if strcmp(dataset,'2-20kHz')
    decodingError_LMled = mean(Decoder_OutputLMled.delta,2);
    decodingErrorCI_LMled = prctile(Decoder_OutputLMled.delta,[5 95],2);
    
    decodingError_LMspk = mean(Decoder_OutpuLMspkt.delta,2);
    decodingErrorCI_LMspk = prctile(Decoder_OutpuLMspkt.delta,[5 95],2);
elseif strcmp(dataset,'2-80kHz')
    decodingError_LMled = mean(Decoder_Output_CBA_80kHz.delta,2);
    decodingErrorCI_LMled = prctile(Decoder_Output_CBA_80kHz.delta,[5 95],2);
    
    decodingError_LMspk = mean(Decoder_Output_Bl6_80kHz.delta,2);
    decodingErrorCI_LMspk = prctile(Decoder_Output_Bl6_80kHz.delta,[5 95],2);
end
x=-20:10:100; y = -20:20:20;
[X,Y] = meshgrid(x,y);
Xlin = reshape(X',39,1); Ylin = reshape(Y',39,1);
for i = 1:39
    for ii = 1:39
        chance(i,ii) =  norm([Xlin(i),Ylin(i)]-[Xlin(ii),Ylin(ii)]);
    end
end
chance_avg = mean(chance,'all');

% -- plot
figure; hold on
plot([nAxons(1) nAxons(end)],[chance_avg chance_avg],'k:')
fill([nAxons flip(nAxons)],[decodingErrorCI_ACaud(:,1)' flip(decodingErrorCI_ACaud(:,2))'] ,[0 0 .5],'edgecolor','none','facealpha',0.2)
plot(nAxons,decodingError_ACaud,'b');
fill([nAxons flip(nAxons)],[decodingErrorCI_LMled(:,1)' flip(decodingErrorCI_LMled(:,2))'] ,[.5 0 0],'edgecolor','none','facealpha',0.2)
plot(nAxons,decodingError_LMled,'r');
fill([nAxons flip(nAxons)],[decodingErrorCI_LMspk(:,1)' flip(decodingErrorCI_LMspk(:,2))'] ,[0 0 0],'edgecolor','none','facealpha',0.2)
plot(nAxons,decodingError_LMspk,'k');

ylim([0 75]); yticks([0 25 50 75]);ylabel('error')
xticks(nAxons); xlabel('number of axons')


% % - save
sgtitle(dataset)
set(gcf,'units','normalized','position',[.2 .1 .7 .7])
savefig(gcf,[SaveFolder filesep 'ErrorVsNaxons'])
saveas(gcf,[SaveFolder filesep 'ErrorVsNaxons.tif']);
saveas(gcf,[SaveFolder filesep 'ErrorVsNaxons.svg']);
end
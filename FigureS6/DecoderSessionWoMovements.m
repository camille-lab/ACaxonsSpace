load('C:\Users\camil\AudSpace\Data\Decoder\ACaud_2-80kHz_allSessionsWithoutMovement.mat')
Decoder_Output2 = Decoder_Output;
mean_deltaAz = mean(Decoder_Output2.delta,2);
CI_deltaAz = prctile(Decoder_Output2.delta,[5 95],2);
nAxons = [10 25 50 100 200 252];
figure;hold on
fill([nAxons flip(nAxons)],[CI_deltaAz(:,1); flip(CI_deltaAz(:,2))],'b','edgecolor','none','facealpha',.2)
plot(nAxons,mean_deltaAz,'r')
ylim([0 100])

x=-20:10:100; y = -20:20:20;
[X,Y] = meshgrid(x,y);
Xlin = reshape(X',39,1); Ylin = reshape(Y',39,1);
for i = 1:39
    for ii = 1:39
        chance(i,ii) =  norm([Xlin(i),Ylin(i)]-[Xlin(ii),Ylin(ii)]);
    end
end
chance_avg = mean(chance,'all');
plot([nAxons(1) nAxons(end)],[chance_avg chance_avg],'k:')

xlabel('Number of axons');xticks(nAxons)
ylabel('Decoding error');ylim([20 60])
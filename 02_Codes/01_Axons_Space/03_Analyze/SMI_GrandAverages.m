load('C:\Users\camil\AudSpace\Data\SMI\GrandAverages\Wilcoxon_dot01_dot15ampTh_0sd\SMI_grandAverages.mat');
figure;
subplot(1,3,2); hold on
scatter(ones(size(smi_ACaud,1),1),nanmean(smi_ACaud,2),'filled')
scatter(3*ones(size(smi_LMaud,1),1),nanmean(smi_LMaud,2),'filled')
scatter(2*ones(size(smi_LMvis,1),1),nanmean(smi_LMvis,2),'filled')
plot([0.7 0.9],[nanmean(nanmean(smi_ACaud)) nanmean(nanmean(smi_ACaud))],'k','linewidth',3)
sem = std(nanmean(smi_ACaud,2))./sqrt(size(smi_ACaud,1));
plot([.8 .8],[nanmean(nanmean(smi_ACaud))-sem nanmean(nanmean(smi_ACaud))+sem],'k','linewidth',2)
plot([2.1 2.3],[nanmean(nanmean(smi_LMvis)) nanmean(nanmean(smi_LMvis))],'k','linewidth',3)
sem = std(nanmean(smi_LMvis,2))./sqrt(size(smi_LMvis,1));
plot([2.2 2.2],[nanmean(nanmean(smi_LMvis))-sem nanmean(nanmean(smi_LMvis))+sem],'k','linewidth',2)
plot([3.1 3.3],[nanmean(nanmean(smi_LMaud)) nanmean(nanmean(smi_LMaud))],'k','linewidth',3)
sem = std(nanmean(smi_LMaud,2))./sqrt(size(smi_LMaud,1));
plot([3.2 3.2],[nanmean(nanmean(smi_LMaud))-sem nanmean(nanmean(smi_LMaud))+sem],'k','linewidth',2)

ylim([0 1]);yticks([0 0.5 1]);
xlim([.5 3.5]); xticks([1 2 3])
set(gca,'Tickdir','out')
% smi_LMvis2 = NaN(8*7,1);smi_LMvis2(1:18) = smi_LMvis(:);
% smi_LMaud2 = NaN(8*7,1);smi_LMaud2(1:18) = smi_LMaud(:);
% [p,tbl,stats]=anova1([smi_ACaud(:),smi_LMvis2,smi_LMaud2],[],'off');
smi_LMvis2 = NaN(8,7);smi_LMvis2([1:3],[1:6]) = smi_LMvis;
smi_LMaud2 = NaN(8,7);smi_LMaud2([1:3],[1:6]) = smi_LMaud;
[p,tbl,stats]=anova1([nanmean(smi_ACaud,2),nanmean(smi_LMvis2,2),nanmean(smi_LMaud2,2)],[],'off');
c=multcompare(stats,'Display','off');
zetitle = sprintf('1W ANOVA, F(%.3g,%.3g)=%.5g, p=%.3g',tbl{2,3},tbl{3,3},tbl{2,5},p);
zetitle2= sprintf('post-hoc:1vs2=%.4g, 1vs3=%.4g, 2vs3=%.4g',c(1,6),c(2,6),c(3,6));
title([zetitle newline zetitle2])

subplot(1,3,1); hold on
scatter(ones(size(smi_ACaud,1)*size(smi_ACaud,2),1),smi_ACaud(:),'filled')
scatter(2*ones(size(smi_LMaud,1),1),nanmean(smi_LMaud,2),'filled')
scatter(2*ones(size(smi_LMvis,1)*size(smi_LMvis,2),1),smi_LMvis(:),'filled')
ylim([0 1]);yticks([0 0.5 1]); ylabel('SMI - grand average')
xlim([.5 2.5]);xticks([1 2])
set(gca,'Tickdir','out')
[~,p]=ttest2(smi_ACaud(:),smi_LMvis(:));
title(['ttest,p=' num2str(p,4)]);

subplot(1,3,3); hold on
scatter(ones(size(smi_ACaud,1),1),nanmean(smi_ACaud,2),'filled')
scatter(2*ones(size(smi_LMvis,1),1),nanmean(smi_LMvis,2),'filled')
plot([0.7 0.9],[nanmean(nanmean(smi_ACaud)) nanmean(nanmean(smi_ACaud))],'k','linewidth',3)
sem = std(nanmean(smi_ACaud,2))./sqrt(size(smi_ACaud,1));
plot([.8 .8],[nanmean(nanmean(smi_ACaud))-sem nanmean(nanmean(smi_ACaud))+sem],'k','linewidth',2)
plot([2.1 2.3],[nanmean(nanmean(smi_LMvis)) nanmean(nanmean(smi_LMvis))],'k','linewidth',3)
sem = std(nanmean(smi_LMvis,2))./sqrt(size(smi_LMvis,1));
plot([2.2 2.2],[nanmean(nanmean(smi_LMvis))-sem nanmean(nanmean(smi_LMvis))+sem],'k','linewidth',2)

ylim([0 1]);yticks([0 0.5 1]);
xlim([.5 2.5]); xticks([1 2 3])
set(gca,'Tickdir','out')
[~,p]=ttest2(nanmean(smi_ACaud,2),nanmean(smi_LMvis,2));
title(['ttest,p=' num2str(p,4)])

%%
SaveFolder_GrandAverage = ['C:\Users\camil\AudSpace\Plots\SMI' filesep 'GrandAverage' filesep];
if ~exist(SaveFolder_GrandAverage,'dir')
    mkdir(SaveFolder_GrandAverage)
end
saveas(gcf,[SaveFolder_GrandAverage 'SMI_GA.tif']);
saveas(gcf,[SaveFolder_GrandAverage 'SMI_GA.svg']);
savefig([SaveFolder_GrandAverage 'SMI_GA']);


load('C:\Users\camil\AudSpace\Data\SMI\Wilcoxon_dot01_0sd_0dot15ampTh\SMI_ACaud_LMaud_LMvis.mat');
%%
clearvars -except SMI_LMvis SMI_LMaud SMI_ACvis SMI_ACaud 

smi_bl6_ACaud = cell(size(SMI_ACaud,1),1);
for i = 1:size(SMI_ACaud,1)
    smi_bl6_ACaud{i} = cat(1,SMI_ACaud{i,:});
end


for i = 1:size(SMI_LMaud,1)
    smi_bl6_LM20kHz{i} = cat(1,SMI_LMaud{i,:});
end

for i = 1:size(SMI_LMvis,1)
    smi_bl6_LMvis{i} = cat(1,SMI_LMvis{i,:});
end

%% cumsum
figure; hold on

for i = 1:length(smi_bl6_ACaud)
p1 = cdfplot(smi_bl6_ACaud{i});
p1.Color = [0.6350 0.0780 0.1840];
p1.LineWidth = .05;
end
smi_bl6_20kHz_all = cat(1,smi_bl6_ACaud{:});
p2 = cdfplot(smi_bl6_20kHz_all);
p2.Color = [0.6350 0.0780 0.1840];
p2.LineWidth = 5;

for i = 1:length(smi_bl6_LM20kHz)
p1 = cdfplot(smi_bl6_LM20kHz{i});
p1.Color = [0.3010 0.7450 0.9330];
p1.LineWidth = .05;
end
smi_bl6_LM20kHz_all = cat(1,smi_bl6_LM20kHz{:});
p3 = cdfplot(smi_bl6_LM20kHz_all);
p3.Color = [0.3010 0.7450 0.9330];
p3.LineWidth = 5;

for i = 1:length(smi_bl6_LMvis)
p1 = cdfplot(smi_bl6_LMvis{i});
p1.Color = [0 0.4470 0.7410];
p1.LineWidth = .05;
end
smi_bl6_LMvis_all = cat(1,smi_bl6_LMvis{:});
p4 = cdfplot(smi_bl6_LMvis_all);
p4.Color = [0 0.4470 0.7410];
p4.LineWidth = 5;

legend([p2,p3,p4],{'Bl6 2-20kHz','LM 2-20kHz','LM Vis','AC Vis'},'location','northwest')
xlabel('SMI')
ylabel('Cum. fraction')


%% plot mean per mouse
for i = 1:8
meanACaud(i) = nanmean(smi_bl6_ACaud{i});
end
for i= 1:3
meanLMaud(i) = nanmean(smi_bl6_LM20kHz{i});
 meanLMvis(i) = nanmean(smi_bl6_LMvis{i});   
end

% % - anova
y = NaN(8,3);
y(:,1) = meanACaud; y(1:3,2) = meanLMaud; y(1:3,3) = meanLMvis;
[p,tbl,stats]= anova1(y,[],'off');
if p<0.05
    c = multcompare(stats,'displayopt','off');
end

% %  - plot mean across mice
figure;
subplot(1,3,[1 2]);hold on
for i = 1:3
scatter(i*ones(8,1),y(:,i),'filled')
end
nMice = sum(~isnan(y));
sem = std(y,[],'omitnan')./sqrt(nMice);
plot([1:3;1:3],[nanmean(y)-sem;nanmean(y)+sem],'k')
plot([.9:2.9;1.1:3.1],[nanmean(y);nanmean(y)],'k')
xlim([.5 3.5]); xticks([1 2 3])
ylim([0 1]);yticks([0 .5 1])
%
title1 = sprintf('1W ANOVA F(%.3g,%.3g)=%.3g, p=%.3g',...
    cell2mat(tbl(2,3)),cell2mat(tbl(3,3)),cell2mat(tbl(2,5)),p);
title(title1)
subplot(1,3,3)
text(0,1,num2str(c(:,[1 2 6]),3),'horizontalalignment','left','verticalalignment','top')
axis off

%% k-sample anderson darling test
% docu: https://cran.r-project.org/web/packages/kSamples/kSamples.pdf
Y= [smi_bl6_LMvis_all,ones(5964,1)];
Y= [Y;smi_bl6_LM20kHz_all,2*ones(1271,1)];
Y= [Y;smi_bl6_20kHz_all,3*ones(5632,1)];
AnDarksamtest(Y)
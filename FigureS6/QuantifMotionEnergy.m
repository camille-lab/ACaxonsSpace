saveFolder = 'C:\Users\camil\AudSpace\Plots\Face';

trialLength = 5; % (s)
stimOnset = 1;   % (s)
tAna = [0.2 1.2];
tBase = [-1 0];

motionBase_all = NaN(8,10);
motionResp_all = NaN(8,10);
isBehavior = false(8,10);
for dataset = 1:2
    disp('loading data...')
    if dataset == 1
        load('C:\Users\camil\Data\Face\SpeakerSpace\BL6\faceData.mat')
    else
        load('C:\Users\camil\Data\Face\SpeakerSpace\CBA\faceData.mat')
    end    
    [n_mice,n_pos] = size(faceData.video_motion);      

    for mouse = 1:n_mice
        for positions = 1:n_pos
            if ~ isempty(faceData.video_motion{mouse,positions})
                data = faceData.video_motion{mouse,positions};
                
                fps = faceData.fps{mouse,positions};
                x_vect = -stimOnset:1/fps:trialLength;
                analysis_window = false(size(data,3),1);
                analysis_window(x_vect>tAna(1) & x_vect<tAna(2)) = true;
                base_window = false(size(data,3),1);
                base_window(x_vect>tBase(1) & x_vect<tBase(2)) = true;
                
                temp1 = mean(data(:,:,base_window),3);
                temp2 = mean(data(:,:,analysis_window),3);
                
                baseline = mean(mean(temp1));
                response = mean(mean(temp2));
                
                base_test = temp1(:);
                resp_test = temp2(:);
                [~,p] = ttest(base_test,resp_test);
                if p<0.05 && baseline<response
                    if dataset == 1
                        isBehavior(mouse,positions) = true;
                    else
                        isBehavior(mouse+2,positions) = true;
                    end
                end
                
                if dataset == 1
                    motionBase_all(mouse,positions) = baseline;
                    motionResp_all(mouse,positions) = response;
                else
                    motionBase_all(mouse+2,positions) = baseline;
                    motionResp_all(mouse+2,positions) = response;
                end
            end
        end
    end
end

%% Figure
base = motionBase_all(:);
resp = motionResp_all(:);
isBehav = isBehavior(:);

figure;
subplot(1,2,1); hold on
isBehav_data = [base(isBehav) resp(isBehav)]';
plot([ones(1,sum(isBehav));2*ones(1,sum(isBehav))],[isBehav_data],'k')
noBehav_data = [base(~isBehav) resp(~isBehav)]';
plot([3*ones(1,sum(~isBehav));4*ones(1,sum(~isBehav))],[noBehav_data],'k')
plot([.9 1.9 2.9 3.9;1.1 2.1 3.1 4.1],[[mean(isBehav_data,2) mean(isBehav_data,2)]' [nanmean(noBehav_data,2) nanmean(noBehav_data,2)]'],'r','linewidth',2)
[~,p1] = ttest(base(isBehav) , resp(isBehav));
[~,p2] = ttest(base(~isBehav) ,  resp(~isBehav));
nNoBehav = sum(~isnan(noBehav_data(1,:)));
nBehav = sum(isBehav);
for i=1:100
   temp = randsample(nBehav,nNoBehav,true);
   y = isBehav_data(:,temp);
   relChange(i) = mean(y(2,:)-y(1,:));
   [~,p_temp] = ttest(y(1,:),y(2,:));
   p_resamp(i) = p_temp;
end
title(['Behav, ' num2str(min(p_resamp),3) '<p<' num2str(max(p_resamp),3) '| NoBehav,p=' num2str(p2,3)])
xlim([.5 4.5]); xticks([1:4]); xticklabels({'base','sound','base','sound'});
ylabel('motion energy (z)')

subplot(1,2,2); hold on
histogram(relChange);
plot([nanmean(resp(~isBehav) - base(~isBehav)) nanmean(resp(~isBehav) - base(~isBehav))],[0 35],'r')
xlabel('Behav change w/ sound')
title('mean noBehav (red) vs distribution from resampled Behav (blue)')

% - save
set(gcf,'units','normalized','position',[.2 .1 .7 .7])
savefig(gcf,[saveFolder filesep 'MotionEnergy'])
saveas(gcf,[saveFolder filesep 'MotionEnergy.tif']);
saveas(gcf,[saveFolder filesep 'MotionEnergy.svg']);
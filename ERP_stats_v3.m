function [] = ERP_stats_v3(EEG,ERP1,ERP2,label1,label2,channels,n_permutes)

% plots two ERPs with their conf interval
% computes permutations t test 
% uses cluster correction for multiple comparison correction 

% Inputs:
% EEG           - an EEGLab EEG structure file where the data come from for
%                 channal locations and times
% ERP1,ERP2     - 3D data (channals, times, participants) of conditions 1/2
%                 or hits/CRs
% label1,label2 - names of the two conditions 
% channels      - channels to take the ERP data of 

% update: 20191016
% it only does permutation test on a specified region, instead of on the
% entire signal

% update: 20191017
% has small edits compared to v2
% titles are A and B align to left 
% lines are stacked so that the ERP is not behind the basline line 

%% (not so) stupid way to make a string for title from the string cells
for t = 1:numel(channels)
    if t==1
        tit = channels{t};
    else
        tit = strcat(tit, {' '}, channels{t});
    end
end
tit = char(tit);


%% Read the electrodes in the channel group
tic
for i = 1:size(channels,2)
    
    if (i==1)
        Cond_1 = zeros(size(channels,2),EEG.pnts,size(ERP1,3));
        Cond_2 = zeros(size(channels,2),EEG.pnts,size(ERP2,3));
    end
    
    % and find the index (channel number) of that channel
    % first input channal to plot second input all channels
    channel_index = strcmpi(channels(i),{EEG.chanlocs.labels});
    
    % put the error for not finding the channal here
    if (sum(channel_index)==0); error(['Electrode ' channels{i} ' not found']); end
    
    Cond_1(i,:,:) = ERP1(channel_index,:,:); %Hits
    Cond_2(i,:,:) = ERP2(channel_index,:,:); %CRs
end

% average per electrode groups and reshape it into P x T (participants by times)
Cond_1 = squeeze(mean(Cond_1,1))';
Cond_2 = squeeze(mean(Cond_2,1))';
disp(['* take data from defined channels [' num2str(toc) ' secs]'])



%% keep relevant times of the ERP (instead of -1000 to 2000) (this is most relevant for permutation analysis)


% the part of the signal for permuation test. 
t1      = dsearchn(EEG.times',100);
t2      = dsearchn(EEG.times',700);

t_window       = EEG.times; % times
t_window       = t_window(t1:t2);
Cond_1_window  = Cond_1(:,t1:t2);
Cond_2_window  = Cond_2(:,t1:t2);


% still dont even need the entire signal for ploting, only use times -200 to 1000 
t1      = dsearchn(EEG.times',-200);
t2      = dsearchn(EEG.times',1000);

t       = EEG.times; % times
t       = t(t1:t2);
Cond_1  = Cond_1(:,t1:t2);
Cond_2  = Cond_2(:,t1:t2);

%% Plot the two conditions
tic

subplot(211)


line([0,0],[-1000 2000],'color',[.5 .5 .5],'LineWidth',.5)
line([-200 1200],[0,0],'color',[.5 .5 .5],'LineStyle','--','LineWidth',1)

x  = t;
y1 = mean(Cond_1);
y2 = mean(Cond_2);

% standard error of the mean 
y1_sem = std(Cond_1)/sqrt(size(Cond_1,1));
y2_sem = std(Cond_2)/sqrt(size(Cond_1,1));

X_plot   = [x, fliplr(x)];
Y1_plot  = [y1-1.96.*y1_sem, fliplr(y1+1.96.*y1_sem)]; % compute 95% CI
Y2_plot  = [y2-1.96.*y2_sem, fliplr(y2+1.96.*y2_sem)];

hold on
set(gca,'Ydir','reverse') % since other plots of reversed (still not sure why we plot in reverse)

f1 = fill(X_plot, Y1_plot , 1,....
    'facecolor','red', ...
    'edgecolor','none', ...
    'facealpha', 0.1);
p1 = plot(x,y1,'LineWidth',1);

f2 = fill(X_plot, Y2_plot , 1,....
    'facecolor','blue', ...
    'edgecolor','none', ...
    'facealpha', 0.1);
p2 = plot(x,y2,'LineWidth',1);

% t = title(['A) ERPs averaged across ' ,tit]);
% set(t, 'units', 'normalized'); % instead of units based on data, this way if the axis/data are different you still put the godamn text in the same place
% set(t, 'position', [.05 1.08]); % first value x then y, so first value left right, the other up down

t = title('A)');
set(t, 'units', 'normalized'); % instead of units based on data, this way if the axis/data are different you still put the godamn text in the same place
set(t, 'position', [-.1 1.08]); % first value x then y, so first value left right, the other up down


% axis limits
x_axis_limit = [-200 1000];  % in ms
y_axis_limit(1) = floor(min([Y1_plot Y2_plot].*2))/2; % y min lim micro volts
y_axis_limit(2) = ceil(max([Y1_plot Y2_plot].*2))/2; % y max lim

set(gca,'xlim',x_axis_limit,'ylim',y_axis_limit);
box off
xlabel('Time (ms)')
ylabel('Voltage (\muV)')

% put legend on southeast or northeast depending on which side must of the ERP falls
if mean(Y1_plot)>0, position='southeast'; else, position='northeast'; end

legend([f1 p1 f2 p2], [label1 ' 95% CI' ],[label1 ' grnd avg'],[label2 ' 95% CI' ],[label2 ' grnd avg'],...
    'Location',position); 
legend boxoff

disp(['* plot two conditions [' num2str(toc) ' secs]'])


%% 
%% set some colours because they look nice 
% TODO for 
colour0 = [.0, .28, .73];  % absolute Zero 
colour1 = [.49, .73, .91]; % Aero 
colour2 = 'r'; % red

%% Do permutation test

% Compute the null distribution (i think aka monte carlo approximation)
% obs_diff = Cond_2 - Cond_1;
% obs_diff =  Cond_1 - Cond_2;
obs_diff =  Cond_1_window - Cond_2_window; % only inlude window of the signal that we want to do permutation test on 

% figure;clf;plot(t,mean(obs_diff)) % plot the difference if you want 

% p-value
pval = 0.05;

% convert p-value to Z value
zval = abs(norminv(pval));

if isempty(n_permutes)
n_permutes = 500;
end 

subjN = size(Cond_1,1);

for permi = 1:n_permutes
    
    if permi==1
        % initialize null hypothesis maps
%         permmaps = zeros(n_permutes,size(Cond_1,2));
        permmaps = zeros(n_permutes,size(Cond_1_window,2));
    end
    %     RandSubjs = 0.5<rand(1,12);
    
    RandSubjs = randi(2,1,subjN)-1; % make an array containing random 0s and 1s that is used to multiply random particiapts data by -1
    
    Null_diff = obs_diff;
%     Null_diff(:,:,RandSubjs==1) = - Null_diff(:,:,logical(RandSubjs)); % fro frequency maps
    Null_diff(RandSubjs==1,:) = - Null_diff(logical(RandSubjs),:);
    
    permmaps(permi,:) = mean(Null_diff,1);
    
end
% 
diffmap = mean(obs_diff); % average the difference map over participants 

%% Show non-corrected thresholded maps

% compute mean and standard deviation maps
mean_h0 = mean(permmaps);
std_h0  = std(permmaps);

% now threshold real data...
% first Z-score
zmap = (diffmap-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;

% now some plotting...
% figure(2), clf

% plot the difference wave 
% subplot(212)
% plot(t,diffmap,'color',colour0)
% title('the observed difference')
% box off

% plot the significant parts in red 
subplot(212)
Y = diffmap;

plot(t_window, Y, 'color', [.55 .57 .67],'LineWidth',1) % plot the difference signal 
hold on 

Y(~logical(zmap)) = NaN; % this is to plot only the significant parts 
% plot the signifcatn part in a different colour 
plot(t_window, Y, 'color', [1 .65 .0],'LineWidth',1) % chrome yellow 

% t = title(['B) Difference wave (' label1 ' minus ' label2 ')']);
% set(t, 'units', 'normalized'); % instead of units based on data, this way if the axis/data are different you still put the godamn text in the same place
% set(t, 'position', [ -130 1.2]); 

t = title('B)');
set(t, 'units', 'normalized'); % instead of units based on data, this way if the axis/data are different you still put the godamn text in the same place
set(t, 'position', [-.1 1.08]); % first value x then y, so first value left right, the other up down


set(gcf,'color','w');
box off


%% Cluster based correction - two types 
% 1) based on the size of the cluser 
% 2) based on the maximum value of the cluster

max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max

% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
    
    % get extreme values (smallest and largest)
    temp = sort( reshape(permmaps(permi,:,:),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
    
end


%% 1) based on the size of the cluser

%% Show histograph of maximum cluster sizes

% TODO useful might want to include this in figure later 
% % figure(3), clf
% % hist(max_cluster_sizes,20);
% % xlabel('Maximum cluster sizes'), ylabel('Number of observations')
% % title('Expected cluster sizes under the null hypothesis')

% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

%% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end

% plot the significant parts in red
subplot(212)
Y = diffmap;
Y(~logical(zmap)) = NaN; % this is to plot only the significant parts 
plot(t_window, Y, 'color', 'r') % plot the signifcatn part in a different colour 

% features of the second plot 
line([0,0],[-1 1],'color',[.5 .5 .5],'LineWidth',.5)
line([-200 1000],[0,0],'color',[.5 .5 .5],'LineStyle','--','LineWidth',1)
xlabel('Time (ms)')
ylabel('Voltage (\muV)')

% put legend on southeast or northeast depending on which side must of the ERP falls
if mean(diffmap)>0, position='southeast'; else, position='northeast'; end

legend('Difference wave','time points p<.05','clusters p<.05',...
    'Location',position);
legend boxoff



end


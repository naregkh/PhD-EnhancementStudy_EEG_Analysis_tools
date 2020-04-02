function [] = TF_Cluster_Based_Permutation_v2(Struct)

% Clusterd-based permutation test 
% for a within participants effect 

% % Input
% TF struct     

% read the data from the struct
NSubj       = size(Struct.Hits,3);
times       = Struct.times;
frex        = Struct.frex;

%% cut the data and do permutation on useful areas only (no reason to run
% it pre 0, since we cant have any reason to expect anything there,
% obviously try it for sanity check)

t1      = dsearchn(times',0);
% t2      = dsearchn(times',500);
t2      = dsearchn(times',1000);
% t2      = dsearchn(times',1200);

times         = times(t1:t2);
Struct.Hits   = Struct.Hits(:,t1:t2,:);
Struct.CRs    = Struct.CRs(:,t1:t2,:);

%% Part that runs the permutations of null hypothesis - creating the distribution of data given null hypothesis

% Compute the difference between conditions (for all participants)

% obs_diff = Struct.CRs - Struct.Hits; 
% positive diff -> desynchronization for hits 
% negative diff -> synchronization for hits 

obs_diff = Struct.Hits - Struct.CRs;
% negative diff -> synchronization for hits 
% positive diff -> desynchronization for hits 

% % plot the difference for sanity check
% % figure;
% % Struct.clim = [-1 1]; % set colour limit of the plot
% % TF_MikeXC_plot(mean(Obs_diff,3),Struct,['Diff ',Struct.chanName]) %
% check code for ERP to get channel names in a string 

% p-value
pval = 0.05;
% pval = 0.1;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 1000;


for permi = 1:n_permutes
    
    if permi==1
        % initialize null hypothesis maps
        permmaps = zeros(n_permutes,size(frex,2),size(times,2));
    end
    
    RandSubjs = randi(2,1,NSubj)-1; % make an array containing random 0s and 1s that is used to multiply random particiapts data by -1
    
    Null_diff = obs_diff;
    Null_diff(:,:,RandSubjs==1) = - Null_diff(:,:,logical(RandSubjs));
    
    permmaps(permi,:,:) = mean(Null_diff,3);
    
end

diffmap = mean(obs_diff,3); % average the difference map over participants

%% Show non-corrected thresholded maps

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

% now threshold real data...
% first Z-score
zmap = (diffmap-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;


% top plot - signigficant regions 

clim = [0 1];

figN = get(gcf,'Number');

subplot(311) % For ploting two types of cluster correction 
imagesc(times,frex,diffmap);
hold on
contour(times,frex,logical(zmap),1,'linecolor','k');
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
title('Power values and outlined significance regions')

colorbar
h = colorbar;
ylabel(h, 'Power (dB)')


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
%% Cluster based correction - two types
% 1) Based on the size of the cluser
% 2) Based on the maximum value of the cluster

% permutes clusters - collects cluster sizes and clusters max values 

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

% Now we have the cluster sizes and max cluster value distribution 

%% 1) Based on the size of the cluser
% i think this is testing a one sided hypothesis, 
% if the size of the cluster is larger than it would be expected if Null is
% true

%% Show histograph of maximum cluster sizes

% % % figure(401), clf
% % % hist(max_cluster_sizes,20);
% % % xlabel('Maximum cluster sizes'), ylabel('Number of observations')
% % % title('distribution of cluster sizes under the null hypothesis')

% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));
% % % 
% % % line([cluster_thresh,cluster_thresh],[0 200],'color',[.5 .5 .5],'LineWidth',.5)


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

% middle figure - cluster size corrected 

figure(figN)

subplot(312)
imagesc(times,frex,diffmap)
hold on
contour(times,frex,logical(zmap),1,'linecolor','k')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('After cluster-size correction')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

colorbar
h = colorbar;
ylabel(h, 'Power (dB)')


%% 2) based on the maximum value of the cluster

%% Now with max-pixel-based thresholding

thresh_lo = prctile(max_val(:,1),100-100*(1-pval/2)); % for a two tailed test 
thresh_hi = prctile(max_val(:,2),    100*(1-pval/2));

% % % % find the threshold for lower and upper values
% % % thresh_lo = prctile(max_val(:,1),    100*pval); % what is the
% % % thresh_hi = prctile(max_val(:,2),100-100*pval); % true p-value?
% % % % note about the above code: a 2-tailed test actually requires pval/2 on each tail;
% % % % thus, the above is actually testing for p<.1 !

% % plot the min max value distribution
% % figure(400), clf
% % subplot(211)
% % hist(max_val(:,1),30);
% % xlabel('Minimum pixel difference'), ylabel('Number of observations')
% % title('Distribution of Max pixel value under the null hypothesis')
% % line([thresh_lo,thresh_lo],[0 100],'color','r','LineWidth',.5)
% % 
% % subplot(212)
% % hist(max_val(:,2),30);
% % xlabel('Maximum pixel difference'), ylabel('Number of observations')
% % title('Distribution of Min pixel value under the null hypothesis')
% % line([thresh_hi,thresh_hi],[0 100],'color','r','LineWidth',.5)


% bottom plot - max pixel corrected 

% threshold real data
zmap = diffmap;
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;

figure(figN)
subplot(313)
imagesc(times,frex,diffmap)
hold on
contour(times,frex,logical(zmap),1,'linecolor','k')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('After maximum-pixel-value correction)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','normal')

set(gcf,'color','w');
colorbar
h = colorbar;
ylabel(h, 'Power (dB)')


end


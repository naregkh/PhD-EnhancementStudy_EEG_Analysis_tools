function [FT] = Extract_Power(FT,t,f)
%% function to extract the power
% struct is nareg's tf structures
% t is the times to extract
% f is the frequencies to extract
% it extracts for Hits and CRs

% inputs
% t = [180 250]; % time window
% f = [45  76]; % frequency band

frex = FT.frex;
times = FT.times;

%% need to make sure the f and t are within the range
if ((t(1) <= min(times)) || (t(2) >= max(times)))
    error('time is not in the range (NK)')
elseif ((f(1) <= min(frex)) || (f(2) >= max(frex)))
    error('frequency is not in the range (NK)')
end

% find the time and frequency indices of the tf
f_index = dsearchn(frex',f');
t_index = dsearchn(times',t');

%% loop for participants


FT.stats.CR   = Extract(FT.CRs);
FT.stats.Hits = Extract(FT.Hits);
FT.stats.ALL  = Extract(FT.ALL);

% % % for debugging take this function out and use matlab debugging 

    function struct = Extract(tf)
        
%         power_time = zeros(size(tf,3),numel(times));
        power_time = zeros(size(tf,3),size(tf,2));
        power_mp = zeros(size(tf,3),1);
        
        for subji = 1:size(tf,3)
            power2average = tf(f_index(1):f_index(2),:,subji); % index the power values
            
            power_time(subji,:) = mean(power2average,1); % average the power values across the power window defined 
            
            power_mp(subji,1) = mean(power_time(subji,t_index(1):t_index(2))); % now average the power across the time window
        end
        struct.power_time = power_time;
        struct.power = power_mp;
        struct.times = t;
        struct.frex  = f;
       
    end


end
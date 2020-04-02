function [times2save,frex,tf] = TF_Decomposition(EEG,channel2use,min_freq,max_freq,num_frex,range_cycles,type)
% code mustly either taken or inspired from MikeXCohen loglinTF.m with  added by myself for debuggin and more functions (NK)
% mikexcohen@gmail.com

%% initial parameters

% vector of time points to save in post-analysis downsampling
times2save = -300:20:1200; % in ms

basetime = [-500 -200];

% frequency parameters
% min_freq =  2;
% max_freq = 50; 
% num_frex = 40;
% other wavelet parameters
% range_cycles = [ 4 10 ]; % original parameters 

% % which channel to plot
% channel2use = 'pz';

%% parameter conversion and other initializations

% % load sampleEEGdata.mat

% time vector converted to indices
times2saveidx = dsearchn(EEG.times',times2save');
basetimeidx   = dsearchn(EEG.times',basetime');

% frequencies vector
% frex = logspace(log10(min_freq),log10(max_freq),num_frex);
frex = linspace(min_freq,max_freq,num_frex);

% wavelet parameters
s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex) ./ (2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;


% FFT parameters
nWave = length(wavtime);
nData = EEG.pnts * EEG.trials; % This line is different from above!!
nConv = nWave + nData - 1;

% initialize output time-frequency data
tf = zeros(length(frex),length(times2save));

%% (non-)phase locked: (NK)

switch type
    case 'all'
    case 'PhaseLocked' % simply average the signal over the trials 
        PhaseLocked      = EEG; % creating a new EEG file for the phase locked acitivity
        PhaseLocked.data = mean(EEG.data,3); % Averaging over trials to get the ERP/phase locked activity
        EEG = PhaseLocked;
    case 'NonPhaseLocked' % main signal minus the phase locked = non phase locked
        PhaseLocked      = EEG; % creating a new EEG file for the phase locked acitivity
        PhaseLocked.data = mean(EEG.data,3); % Averaging over trials to get the ERP/phase locked activity
        NonPhaseLocked = EEG;
        for i=1:size(EEG.data,3)
            NonPhaseLocked.data(:,:,i) = EEG.data(:,:,i) - PhaseLocked.data;
        end
        EEG = NonPhaseLocked;
    otherwise
        error('Incorrect input for type')
end

%% make sure inputed electrode exist
if sum((strcmpi(channel2use,{EEG.chanlocs.labels})))==0
    error(['Electrode [ ' channel2use ' ] does not exits'])
end


%% run convolution

% now compute the FFT of all trials concatenated
alldata = reshape( EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:) ,1,[]);
dataX   = fft( alldata ,nConv );

% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet and get its FFT
    % the wavelet doesn't change on each trial...
    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    % now run convolution in one step
    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);
    
    % and reshape back to time X trials
    as = reshape( as, EEG.pnts, EEG.trials );
    
    % compute power and average over trials
    temppow  = mean( abs(as).^2 ,2);
    tf(fi,:) = 10*log10( temppow(times2saveidx) / mean(temppow(basetimeidx)) );
end



end


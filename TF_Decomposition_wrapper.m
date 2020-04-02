function [TF] = TF_Decomposition_wrapper(TF,FileNames,path,params)
% Wrapper for TF_Decomposition
% Reads all participants and per all conditions
t1 = tic;

%%
min_freq     = params.min_freq;
max_freq     = params.max_freq;
num_frex     = params.num_frex;
range_cycles = params.cycles;
type         = params.type;
channels     = TF.channels;

%%

Response = {'All', 'Hits', 'CRs'};
% Response = {'Hits', 'CRs'};
tf_subject = zeros(num_frex,76,size(FileNames,1));

for respi = 1:numel(Response) % per response
    
    for subji = 1:size(FileNames,1) % per participant
        t2 = tic;
        
        % in case empty spaces at the are a problem use strrep
        toload = [path,char(Response(respi)),'\',(FileNames(subji,:))];
        load(toload); % matlab Load is faster than EEGlab's load
        
        % tf per channel
        for ch = 1:numel(channels)
            thischan = char(channels(ch));
%             [times2save,frex,tf_elect(:,:,ch)] = TF_MikeXC_Params(EEG,thischan,min_freq,max_freq,num_frex,range_cycles,type);
            [times2save,frex,tf_elect(:,:,ch)] = TF_Decomposition(EEG,thischan,min_freq,max_freq,num_frex,range_cycles,type);
        end
        
        % averages over electrodes if more than one
        tf_electMean = mean(tf_elect,3);
        tf_subject(:,:,subji) = tf_electMean;

        disp(['subject ' num2str(subji) ' resp ' num2str(respi) ' completed in ' num2str(toc(t2))])

    end
    
    % save to the TF file 
    switch char(Response(respi))
        case 'All'
            TF.ALL  = tf_subject;
        case 'Hits'
            TF.Hits = tf_subject;
        case 'CRs'
            TF.CRs  = tf_subject;
    end
    
end


TF.times = times2save;
TF.frex = frex;
TF.params = params;

disp(['Tf_wrapper took ', num2str(toc(t1)),' seconds']);



end


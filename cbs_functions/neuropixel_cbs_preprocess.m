function bank_model ...
    = neuropixel_cbs_preprocess(myKsDir, chans_to_use, ops)

r = getOr(ops, 'r', 3); % channel-wise pca dimensionality
nspikes_per_unit = getOr(ops, 'nspikes_per_unit', 100); % number of spikes to analyze per unit (puts a cap on the 

sp = loadKSdir(myKsDir);

bank_model.channel_info = neuropixel_get_channel_info(fullfile(myKsDir, sp.dat_path));

%for now hardcode cluster 
cluster = unique(sp.spikeTemplates);

clu_spike_count = zeros(1, length(cluster));
for iuniq = 1:length(cluster) % loop through unique clusters
    ispike_this = find(sp.spikeTemplates == ...
                       cluster(iuniq));
    spike_times_this = sp.st(ispike_this); 
    clu_spike_count(iuniq) = length(spike_times_this);
end

l_low_num_spikes = clu_spike_count < nspikes_per_unit;

cluster(l_low_num_spikes) = [];

wf_spike_times = ceil(sp.st(ismember(sp.spikeTemplates, cluster)) * 30000);
wf_spike_clusters = sp.spikeTemplates(ismember(sp.spikeTemplates, cluster));


nchan_all = sp.n_channels_dat - 1;

%initialize gwfparams
gwfparams.dataDir = myKsDir;
apD = dir(fullfile(myKsDir, '*ap*.bin'));
gwfparams.fileName = apD(1).name;
gwfparams.dataType = 'int16';
gwfparams.nCh = sp.n_channels_dat;
gwfparams.wfWin = [-30 30];
gwfparams.nWf = nspikes_per_unit;
gwfparams.spikeTimes = wf_spike_times;
gwfparams.spikeClusters = wf_spike_clusters;


disp('getting waveforms');
wf = getWaveFormsModified2(gwfparams);


whitening_mat = readNPY(fullfile(gwfparams.dataDir, 'whitening_mat.npy'));

%% preprocessing
wf_size = size(wf.waveForms);

nunits = wf_size(1);

nchan = wf_size(3);
ntimes = wf_size(4);

working_wf = reshape(wf.waveForms, prod(wf_size(1:3)), ntimes)';

% do median subtraction of a prespike window
disp('subtracting medians');
wf_meds = median(working_wf, 1);
working_wf = bsxfun(@minus, working_wf, wf_meds);

% do high pass filtering
disp('filtering waveforms')
hp_cutoff = 150;
Fs = 30e3; 
[b1, a1] = butter(3, hp_cutoff/Fs*2, 'high'); 
filtered_wf = filter(b1, a1, working_wf);
filtered_wf = reshape(filtered_wf', wf_size);

% do whitening filter
disp('applying whitening filter');
for iunit = 1:nunits
    for ispike = 1:nspikes_per_unit
        filtered_wf(iunit, ispike, :, :) = ...
            whitening_mat'*sq(filtered_wf(iunit, ispike, :, :));        
    end
end

clear working_wf;


%% training channel-wise pca models

chan_map = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               
% 1-indexed

% each channel should have a pca model associated with it.

used_temps = sp.temps(cluster + 1, :, :);

channel_models = struct();

use_all_spikes = false; % use waveforms from all spikes from all units, not just the ones whose templates intersected with the unit

disp('performing pca on each channel');
pca_wf = zeros(nunits, nspikes_per_unit, nchan, r);
for ich = 1:nchan
    if ismember(chan_map(ich), chans_to_use)                    
        
        if use_all_spikes            
            tpls = true(1, nunits);
        else
            tpls = sq(sum(abs(used_temps(:, :, ich)), 2)>0)';
        end
        
        ntpl = sum(tpls);      
        channel_models(ich).tpls = tpls;

        if sum(tpls)            
            try
            selected_wfs = filtered_wf(tpls, :, ich, :);              
            catch e
               keyboard; 
            end
            selected_wfs = reshape(selected_wfs, ntpl*nspikes_per_unit, ntimes); 
                        
            [coeff, score, latent] = pca(selected_wfs);
                                                
            channel_models(ich).coeff = coeff;
            channel_models(ich).latent = latent;
                
            pca_wf(tpls, :, ich, :) = reshape(score(:, 1:r), [ntpl, nspikes_per_unit,1, r]);
            
        end
    
    else
%         fprintf('ignoring ch %d when computing pca\n', chan_map(ich));
   
    end
    
end

pca_wf = permute(pca_wf, [1,2,4,3]);
Xred = reshape(pca_wf, [nunits*nspikes_per_unit, nchan*r]);

ilabels = repmat((1:nunits)', nspikes_per_unit, 1);
N = nunits*nspikes_per_unit;
n = size(Xred, 2);
labels = false(N, nunits);
for i = 1:N
    labels(i, ilabels(i)) = true;
end



%% scatter matrices

[S_Wx, S_bx] = neuropixel_cbs_form_scatter_matrices(Xred, labels, nunits);

% sparsify the within class covariance
% S_Wx = sparse(S_Wx);


%% populate struct ouptut

bank_model.S_Wx = S_Wx;
bank_model.S_bx = S_bx;

bank_model.filtered_wf = filtered_wf;

bank_model.channel_models = channel_models;    
bank_model.Xred = Xred;
bank_model.labels = labels;
bank_model.ilabels = ilabels;    

bank_model.nunits = nunits;
bank_model.nspikes_per_unit = nspikes_per_unit;
bank_model.nchan = nchan;
bank_model.nchan_all = nchan_all;
bank_model.ntimes = ntimes;
bank_model.nspikes_total = N;
bank_model.r = r;
bank_model.chan_map = chan_map;
bank_model.chans_to_use = chans_to_use;
bank_model.cluster = cluster;
bank_model.sp = sp;

bank_model.wf_unit_ids = wf.unitIDs;

bank_model.raw_dir = myKsDir;
    



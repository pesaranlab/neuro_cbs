function [snrs_by_method, med_amps_by_method] = get_channel_neuron_snrs(data_dir, methods, cluster_timing_mat, cluster_mat, ops2)


avg_spike_rate_cutoff = getOr(ops2, 'avg_spike_rate_cutoff', 0);
banks_to_analyze = getOr(ops2, 'banks_to_analyze', 0:2);

nspikes_per_unit = 100;
snrs_by_method = {};
med_amps_by_method = {};
for imethod = 1:length(methods)
    method = methods{imethod};
    
    snrs_all_banks = {};
    med_amps_all_banks = {};
    
    for ibank = 1:length(banks_to_analyze)
        
        
        bank = banks_to_analyze(ibank);
        
        if bank == 2
            chans_to_use = 1:192;
        else
            chans_to_use = 1:384;
        end
        
        rec_path_dir = fullfile(data_dir, sprintf('rec_bank%d_%s', bank, method), sprintf('rec_bank%d_%s_imec0', bank, method));        
        rec_path = fullfile(rec_path_dir, sprintf('rec_bank%d_%s_t0.imec0.ap.bin', bank, method));

        cluster_timing = cluster_timing_mat{imethod, ibank};
        clusters = cluster_mat{imethod, ibank};
                      
        sp = loadKSdir(rec_path_dir);

        
        %% screening neurons by rate
        cluster = unique(sp.spikeTemplates);

        avg_spike_rates = zeros(1, length(cluster));
        last_spike_time = max(sp.st);
        for iuniq = 1:length(cluster) % loop through unique clusters
            ispike_this = find(sp.clu == ...
                               cluster(iuniq));
            spike_times_this = sp.st(ispike_this);
            avg_spike_rates(iuniq) = length(spike_times_this)./last_spike_time;  
        end

        l_low_spike_rate = avg_spike_rates < avg_spike_rate_cutoff;

        cluster(l_low_spike_rate) = [];                
        
        %% analyzing spike amplitudes
                    
        nclus = length(cluster);

        % for each cluster, get cross-trial median amplitudes
        med_amps = zeros(1, nclus);
        for icl = 1:nclus
            clid = cluster(icl);
            i_this_cl = (sp.spikeTemplates == clid);
            amplitudes = sp.tempScalingAmps(i_this_cl);
            med_amps(icl) = median(amplitudes);
        end
        med_amps = med_amps(~isnan(med_amps));
        med_amps_all_banks = [med_amps_all_banks, med_amps];

        
        %% getting waveforms

        wf_spike_times = ceil(sp.st((ismember(int32(sp.spikeTemplates), int32(cluster)))) * 30000);
        wf_spike_clusters = sp.spikeTemplates((ismember(int32(sp.spikeTemplates), int32(cluster))));
        
        %initialize gwfparams
        gwfparams.dataDir = rec_path_dir;
        apD = dir(fullfile(rec_path_dir, '*ap*.bin'));
        gwfparams.fileName = apD(1).name;
        gwfparams.dataType = 'int16';
        gwfparams.nCh = 385;
        gwfparams.wfWin = [-30 30];        
        gwfparams.nWf = nspikes_per_unit;        
        gwfparams.spikeTimes = wf_spike_times;
        gwfparams.spikeClusters = wf_spike_clusters;        
        
        disp('getting waveforms');
        wf = getWaveFormsModified(gwfparams);       
        
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
                
        %% SNR

        chan_map = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               
        % 1-indexed

        % each channel should have a pca model associated with it.

        used_temps = sp.temps(cluster + 1, :, :);

        channel_models = struct();

        use_all_spikes = false; % use waveforms from all spikes from all units, not just the ones whose templates intersected with the unit

        % disp('performing pca on each channel');
        % pca_wf = zeros(nunits, nspikes_per_unit, nchan, r);

        snrs_by_chan = {};
        for ich = 1:nchan
            if ismember(chan_map(ich), chans_to_use)

                % find all templates that this channel is involved in
                if use_all_spikes            
                    tpls = true(1, nunits);
                else
                    tpls = sq(sum(abs(used_temps(:, :, ich)), 2)>0)';
                end
                ntpl = sum(tpls);

                channel_models(ich).tpls = tpls;

                if sum(tpls)
                    % collect all the waveforms for this channel from the different units
                    ftpls = find(tpls);
                    % (ntpl*nspikes_per_unit) x ntimes    

                    selected_wfs = filtered_wf(tpls, :, ich, :);  
                    mean_wf = mean(selected_wfs, 2);

                    centered_selected_wfs = bsxfun(@minus, selected_wfs, mean_wf);
                    signal_rms = sq(sqrt(sum(mean_wf.^2, 4))/sqrt(size(mean_wf, 4)));                       
                    noise_rms = sq(sqrt(sum(sum(centered_selected_wfs.^2, 2), 4)))/sqrt(size(selected_wfs, 2)*size(selected_wfs, 4));

                    snr_db = 20*log10(signal_rms./noise_rms);

                    snrs_by_chan{ich} = snr_db;

                end
            else
                fprintf('ignoring ch %d when computing snr\n', chan_map(ich));
            end
        end       
        
        snrs_all_banks = [snrs_all_banks, snrs_by_chan];
    end
    
    snrs_by_method{imethod} = unfold_cell(snrs_all_banks);
    med_amps_by_method{imethod} = unfold_cell(med_amps_all_banks);
end
    

end


function v = unfold_cell(cell_of_vectors)
nelem = 0;

for i = 1:length(cell_of_vectors)
   nelem = nelem + length(cell_of_vectors{i}); 
end

v = zeros(nelem, 1);
offset = 1;
for i = 1:length(cell_of_vectors)
    this_len = length( cell_of_vectors{i});
    
    v(offset : offset +  this_len-1) = cell_of_vectors{i}(:);

   offset = offset + this_len;
end

end

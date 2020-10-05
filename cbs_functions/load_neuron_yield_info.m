function [neuron_yields, cluster_timing_mat, cluster_mat] = load_neuron_yield_info(data_dir, methods, ops2)

banks_to_analyze = getOr(ops2, 'banks_to_analyze', 0:1) ;  % restricting analysis to a subset of banks that
                         % were determined be at least partially in the
                         % brain.
avg_spike_rate_cutoff = getOr(ops2, 'avg_spike_rate_cutoff', 0);
                         
bin_width = 60; % 2 ms
nCh = 385;
avg_spike_rate_cutoff = 100/120; % Hz

cluster_timing_mat = {};
cluster_mat = {};
neuron_yields = zeros(length(methods), length(banks_to_analyze));
for imethod = 1:length(methods)
    method = methods{imethod};
    for ibank = 1:length(banks_to_analyze)
        
        bank = banks_to_analyze(ibank);
        
        rec_path_dir = fullfile(data_dir, sprintf('rec_bank%d_%s', bank, method), sprintf('rec_bank%d_%s_imec0', bank, method));
        
        rec_path = fullfile(rec_path_dir, sprintf('rec_bank%d_%s_t0.imec0.ap.bin', bank, method));

        clusters =  readNPY(fullfile(rec_path_dir, 'spike_templates.npy')) + 1;       
        sp_times = double(readNPY(fullfile(rec_path_dir, 'spike_times.npy')));
        
        dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
        fileStruct = dir(rec_path);
       
        nSamps = fileStruct.bytes/(nCh*dataTypeNBytes);
        num_bins = ceil(nSamps / bin_width);
        
        num_spikes_cutoff = floor((nSamps / 30000) * avg_spike_rate_cutoff);

        cluster_timing = {};
        cluster_count = 0;
        for k = 1:max(clusters)            
            sp_k = clusters == k;
            num_spikes = sum(sp_k);                        
            if num_spikes > num_spikes_cutoff       
                cluster_count = cluster_count + 1;
                spt = sp_times(sp_k);
                
                spt_bin = zeros(num_bins, 1);                        
                spt_bin(ceil(spt / bin_width)) = 1;                
                spt_bin = sparse(spt_bin); 
                
                cluster_timing{k}.spt = spt;
                cluster_timing{k}.spt_bin = spt_bin;                
            end                        
        end
        
        cluster_mat{imethod, ibank} = clusters;
        
        cluster_timing_mat{imethod, ibank} = cluster_timing;
        neuron_yields(imethod, ibank) = cluster_count;        
    end
end
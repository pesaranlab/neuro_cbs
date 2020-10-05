function bank_models = get_bank_model_neuropixel(raw_files, ops)

% GET_BANK_MODEL_NEUROPIXEL: Performs waveform extraction and filtering for the cbs and mucbs methods.
%   inputs:
%      raw_files: cell array of nbanks recording paths pointing to the `*.ap.bin` file for each bank.
%      ops: options struct
%
%   outputs:
%      bank_models: struct containing filtered waveform data and metadata needed for CBS
%
%   relevant options: 
%      ops.probe_version: currently only the Neuropixel Phase 3b2, aka Phase 1 probe (default: 'phase1')
%      ops.spikesorter: currently only Kilosort2 is supported (default: 'kilosort2')
%      ops.selection_method: available methods: 'cbs', 'mucbs'. mucbs is an approximation to cbs
%                            that uses a diagonal approximation to the within-class covariance.


probe_version = getOr(ops, 'probe_version', 'phase1');
spikesorter = getOr(ops, 'spikesorter', 'kilosort2');
selection_method = getOr(ops, 'selection_method', 'cbs');

switch probe_version
    case 'phase1'
        nbanks = 3;        
        
        assert(length(raw_files)==nbanks, ...
            'for a phase1 neuropixel, raw_files must be a cell array of length of 3');
        
        valid_chans = {};
        for ibank = 1:nbanks
            if ibank == 3
                valid_chans{ibank} = 1:192;
            else
                valid_chans{ibank} = 1:384;
            end        
        end
        
        nchan_all = 384;
        
    otherwise
        error('probe version not supported');
end

switch spikesorter
    case 'kilosort2'                
        
        for ibank = 1:nbanks        
            ks_dir = fileparts(raw_files{ibank});  
            
            switch selection_method
                case 'cbs'
                    ops.do_prewhitening = getOr(ops, 'do_prewhitening',...
                        true);
                    ops.do_channel_pca = getOr(ops, 'do_channel_pca',...
                        true);                   
                                                    
                case 'mucbs'                                       
                    ops.do_prewhitening = getOr(ops, 'do_prewhitening', ...
                        false);
                    ops.do_channel_pca = getOr(ops, 'do_channel_pca', ...
                        false);
                    
                otherwise
                    error('selection method %d not supported', ...
                        selection_method);
            end            
            
            bank_models(ibank) = ...
                neuropixel_cbs_preprocess(ks_dir, valid_chans{ibank}, ops);
            
        end        
    otherwise
        error('spikesorter not supported', spikesorter);        
end


function [is_kilosorted, is_all_one_bank, bank] = check_neuropix_bank_rec(ap_bin_file)

bank = -1;
[imec_dir, imec_file, imec_ext] = fileparts(ap_bin_file);

is_kilosorted = exist(fullfile(imec_dir, 'spike_templates.npy'), 'file') ~= 0;

channel_info = neuropixel_get_channel_info(ap_bin_file);

banks = [channel_info.bank];
is_all_one_bank = all(banks == channel_info(1).bank);
if is_all_one_bank
    bank = banks(end);
    
end

% neuropixel specific?
if any(banks == 2)
    % special treatment of bank 2. 
    if sum(banks ==  2) == 192
        is_all_one_bank = true;
        bank = 2;
    else
        is_all_one_bank = false;        
    end    
end




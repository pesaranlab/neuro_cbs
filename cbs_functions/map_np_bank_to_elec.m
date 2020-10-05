function elec_mapped_scores = map_np_bank_to_elec(scores, bank_model)

nbanks = length(bank_model);
nchan = 384;

elec_mapped_scores = zeros(1, 960);

for ibank = 1:nbanks
        
    chans_to_use = bank_model(ibank).chans_to_use;
    i_ch_select = ismember(bank_model(ibank).chan_map, chans_to_use);
    
    bank_scores = scores{ibank}(i_ch_select);
    
    
    chmap = bank_model(ibank).chan_map(i_ch_select);
    ch_mapped_scores = zeros(1, nchan);    
    ch_mapped_scores(chmap) = bank_scores;
    
    elecs_in_bank = unique(min(960, (ibank-1)*384 + (1:384)));
    elec_mapped_scores(elecs_in_bank) = ch_mapped_scores(1:length(elecs_in_bank));
end
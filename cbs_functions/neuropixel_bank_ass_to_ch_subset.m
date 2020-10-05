function ch_enabled = neuropixel_bank_ass_to_ch_subset(bank_assignment)

nchan = 384;
ch_enabled = {};

for ibank = 1:3
    if ibank == 3
        ch_enabled{ibank} = false(1, 192);
    else        
        ch_enabled{ibank} = false(1, nchan);
    end    
end

for ichan = 1:nchan    
    ch_enabled{bank_assignment(ichan)+1}(ichan) = true;    
end
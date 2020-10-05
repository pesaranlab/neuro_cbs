function h = plot_bank_assignments_img(bank_assignments, bank_models)

channel_infos = {};
for ibank = 1:length(bank_models)
    channel_infos{ibank} = bank_models(ibank).channel_info;
end
nchan = length(bank_assignments);

nbanks = length(channel_infos);

nelecs = 960;

electrodes_by_bank = {};
for ibank = 1:nbanks
    electrodes_by_bank{ibank} = unique( ...
                                min(960, 384*(ibank-1) + (1:384)) );
end

electrode_on = zeros(1, nelecs);
electrode_img = zeros(480, 2);
rows = zeros(1, nchan);
cols = zeros(1, nchan);
for ichan = 1:nchan
    ibank = bank_assignments(ichan) + 1;
    row = channel_infos{ibank}(ichan).row;
    col = channel_infos{ibank}(ichan).col;
    
    electrode_img(row, col) = 1;
    
    elecs_in_bank = electrodes_by_bank{ibank};
    elec = elecs_in_bank(ichan);
    electrode_on(elec) = 1;
end

% electrode_img = flipud(electrode_img);

h = tvimage(electrode_img');
daspect([1 15 1]);
end


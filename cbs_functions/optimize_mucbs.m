function [ampscore_bank_assignments, ampscore_select_img] = ...
    optimize_mucbs(ampscores, bank_models, ops)

probe_version = getOr(ops, 'probe_version', 'phase1');

switch probe_version    
    case 'phase1'
        elec_mapped_ampscores = map_np_bank_to_elec(ampscores, bank_models);
        elec_mapped_ampscore_selections = zeros(1, 960);

        ampscore_bank_assignments = zeros(1, 384);
        for ich = 1:384
            avail_ass = neuropixel_get_avail_bank_assignments(ich);

            cand_scores = zeros(1, length(avail_ass));
            for iass = 1:length(avail_ass)
                ielec = avail_ass(iass)*384 + ich;
                cand_scores(iass) = elec_mapped_ampscores(ielec);
            end

            [~, imax] = max(cand_scores);

            ielec = avail_ass(imax)*384 + ich;
            ampscore_bank_assignments(ich) = avail_ass(imax);

            elec_mapped_ampscore_selections(ielec) = 1;

        end
        
        ampscore_select_img = get_img_from_elec_mapped_scores(elec_mapped_ampscore_selections);
        
    otherwise
        error('probe_version not supported');
end


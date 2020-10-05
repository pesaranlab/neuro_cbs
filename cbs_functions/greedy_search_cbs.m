function [bank_assignments, search_output] = greedy_search_cbs(bank_models, bank_assignments_init, ops)

do_classifier = getOr(ops, 'do_classifier', true);
condense_validation_scores = getOr(ops, 'condense_validation_scores', true); 
banks_to_validate = getOr(ops, 'banks_to_validate', 1:length(bank_models));
use_banded_solver = getOr(ops, 'use_banded_solver', true); 

%% creating sparse representation for banded solver, if used

if use_banded_solver
    for ibank = 1:length(bank_models)    
        S_Wx = bank_models(ibank).S_Wx;        
        bank_models(ibank).S_Wx_sparse = sparse(S_Wx);    
    end
end

%% search procedure iterations

bank_assignments = bank_assignments_init;

if do_classifier
    % initial classification performance
    ch_enabled = neuropixel_bank_ass_to_ch_subset(bank_assignments);
    [~, init_prop_correct] = ...
        neuropixel_cbs_eval_criterion(bank_models, ch_enabled, true, ops);        
end

max_iter = getOr(ops, 'max_iter', 100); % max number of passes through the array.

nchan_all = bank_models(1).nchan_all;

nbanks = length(bank_models);

Jmax = 0;
Jtrace = [];
Jmax_prev = 0;
bank_assignment_record = [];
n_units_per_bank = [bank_models.nunits];

prop_correct_record = zeros(max_iter, nbanks);
for iiter = 1:max_iter    
    
    Jtraj = zeros(1, nchan_all);
    ind=1;
    bank_assignment_chunk = zeros(nchan_all, nchan_all);    
        
    for ch = randperm(nchan_all) % 1 based                
        avail_ass = neuropixel_get_avail_bank_assignments(ch);               
        J = zeros(1, length(avail_ass));        
        current_ass = bank_assignments(ch);
        
        for iass = 1:length(avail_ass)
            if  current_ass == avail_ass(iass)
                bank_assignments(ch) = current_ass;
                J(iass) = Jmax;
            else
                bank_assignments(ch) = avail_ass(iass);
                ch_enabled = neuropixel_bank_ass_to_ch_subset(bank_assignments);     
                
                ops.do_classifier = false;
                J(iass) = neuropixel_cbs_eval_criterion(bank_models, ...
                    ch_enabled, false, ops);        
            end            
        end
        
        [~, imax] = max(J);
                
        Jmax = J(imax);
        Jtraj(ind) = Jmax;
        
        bank_assignments(ch) = avail_ass(imax);            
        bank_assignment_chunk( ind, :) = bank_assignments;
        
        ind = ind + 1;
    end    
                
    bank_assignment_record = [bank_assignment_record; bank_assignment_chunk];        
    fprintf('pass %d  Jmax = %f\n', iiter, Jmax);
    if do_classifier    
        
        ch_enabled = neuropixel_bank_ass_to_ch_subset(bank_assignments);
        [~, prop_correct, prop_correct_global] = ...
            neuropixel_cbs_eval_criterion(bank_models, ch_enabled, true, ...
            ops);        
                        
        fprintf('  proportion correct = %0.1f%%\n', ...
            prop_correct_global*100);
                
        prop_correct_record(iiter, :) = prop_correct;
    end
            
    if Jmax_prev == Jmax        
        break;
    else
        Jmax_prev = Jmax;
    end
    
    Jtrace = [Jtrace, Jtraj];
        
end

search_output.Jtrace = Jtrace;
search_output.bank_assignment_record = bank_assignment_record;
   
%% classifier validation

if do_classifier
    last_prop_correct = prop_correct_record(iiter, :);
    a = bsxfun(@minus, prop_correct_record(1:iiter, :), last_prop_correct);
    % all(a >= 0, 2) 
    
    if condense_validation_scores
        n_units_per_bank = [bank_models.nunits];
        a = a*n_units_per_bank'./sum(n_units_per_bank);
    end

    is_better_than_last = all(a >= 0, 2);

    if any(is_better_than_last(1:iiter-1))
        [~, imax] = max(sum(a(is_better_than_last, :), 2));

        b = find(is_better_than_last);

        best_iter = b(imax);
        best_assignments = bank_assignment_record(best_iter*nchan_all, :);
        bank_assignments = best_assignments;

        best_prop_correct = prop_correct_record(best_iter, :);
    else
        best_prop_correct = last_prop_correct;
    end

    if condense_validation_scores
        if best_prop_correct*n_units_per_bank' < init_prop_correct*n_units_per_bank'
            disp('defaulting to initial condition');
            best_prop_correct = init_prop_correct;
            bank_assignments = init_bank_assignments;
        end
    else
        if all(best_prop_correct < init_prop_correct)        
            disp('defaulting to initial condition');
            best_prop_correct = init_prop_correct;
            bank_assignments = init_bank_assignments;  
        end
    end

    classif_results.init_prop_correct = init_prop_correct;
    classif_results.prop_correct_record = prop_correct_record(1:iiter, :);
    classif_results.best_prop_correct = best_prop_correct;

    search_output.classif_results = classif_results;    
    
    best_prop_correct_global = sum(best_prop_correct(banks_to_validate).*n_units_per_bank(banks_to_validate)/sum(n_units_per_bank(banks_to_validate)));
    
    fprintf('  best proportion correct = %0.1f%%\n', best_prop_correct_global*100);
end


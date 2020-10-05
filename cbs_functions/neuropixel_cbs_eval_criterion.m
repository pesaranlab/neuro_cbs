function [J1, prop_correct_Y, prop_correct_global] = ...
    neuropixel_cbs_eval_criterion(bank_model, ch_enabled, do_classifier, ops)

banks_to_validate = getOr(ops, 'banks_to_validate', 1:length(bank_model));
use_banded_solver = getOr(ops, 'use_banded_solver', true);

prop_correct_Y = 0;

if nargin < 3
   do_classifier = false; 
end

J1 = 0;
for ibank = 1:length(bank_model)

    i_enabled = ch_enabled{ibank};
    
    if any(i_enabled)

        N = size(bank_model(ibank).Xred, 1);                
        
        if use_banded_solver
            S_Wx_full = bank_model(ibank).S_Wx_sparse;           
        else
            S_Wx_full = bank_model(ibank).S_Wx;
        end
        
        S_bx_full = bank_model(ibank).S_bx;
                
        labels = bank_model(ibank).labels;
        ilabels = bank_model(ibank).ilabels;

        Expander = ones(bank_model(ibank).r, 1);
                
        i_possible = bank_model(ibank).chan_map <= length(i_enabled);        
        chan_map = bank_model(ibank).chan_map(i_possible);
                        
        i_valid_enabled = i_enabled(chan_map);
        i_expanded = logical(vec(Expander*i_valid_enabled)');
                
        S_Wx = S_Wx_full(i_expanded, i_expanded);    
        S_bx = S_bx_full(i_expanded, i_expanded);
                
        m = bank_model(ibank).nunits - 1;           
        n = size(S_Wx, 1);        
        nunits = bank_model(ibank).nunits;
        eps = 1e-12;                
        
        S2 = (S_Wx + eps*eye(n));
        
        Sol = S2\S_bx;
                        
        if do_classifier
            Xred = bank_model(ibank).Xred(:, i_expanded);    
            
            % classification statistics
            nfolds = 4;
            test_prob = 0.25;
            accur = zeros(1, nfolds);
            for ifold = 1:nfolds
                
                idx = randperm(N);                
                itest = idx(1:round(N*test_prob));                                
                itrain = setdiff(1:N, itest);
                            
                [A, V, D, trSol, class_means] = ...
                    neuropixel_cbs_train_model(Xred(itrain, :), ...
                    labels(itrain, :), nunits);

                [~, test_ilabels_Y] = ...
                    neuropixel_cbs_test_model(Xred(itest, :), ...
                    class_means, A);

                accur(ifold) = sum(ilabels(itest) == ...
                    test_ilabels_Y)/length(itest);
            end

            prop_correct_Y(ibank) = mean(accur);
        end        
        
        x_space_crit = real(trace(Sol));        
        J1 = J1 + x_space_crit;
                        
    else
        prop_correct_Y(ibank) = 0;
    end
end

n_units_per_bank = [bank_model.nunits];

if do_classifier
    prop_correct_global = ...
        sum(prop_correct_Y(banks_to_validate).*...
        n_units_per_bank(banks_to_validate)/...
        sum(n_units_per_bank(banks_to_validate)));
end

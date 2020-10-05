function selection_map = optimize_selection_map(bank_models, ops)

% OPTIMIZE_SELECTION_MAP: Finds a selection map that optimizes a separability criterion.
%   inputs:
%       bank_models: preprocessed output of get_bank_model_XYZ
%       ops: options
%   outputs:
%       selection_map: final selection map. selection_map is a vector of length
%           nchan_all, the total number of recording channels (not just ones used
%           in the analysis). Each element is the 1-based electrode index for each 
%           channel
%
%   relevant options:
%       ops.init_method: initialization method for cbs.
%           Available methods: 'mucbs' (default), 'checker', 'line'
%       ops.probe_version: currently only 'phase1' probes are supported.
%       ops.do_plot: boolean of whether to plot optimization results.
%       ops.plot_dir: if do_plot is true, this is the directory to
%           save plots to.

init_method = getOr(ops, 'init_method', 'mucbs');
probe_version = getOr(ops, 'probe_version', 'phase1');
do_plot = getOr(ops, 'do_plot', false);
save_plot = getOr(ops, 'save_plot', false);
plot_dir = getOr(ops, 'plot_dir', '.');
banks_to_validate = getOr(ops, 'banks_to_validate', 1:3);
do_classifier = getOr(ops, 'do_classifier', false);
% compute AmpScores

nbanks = length(bank_models);
ampscores = {};
tic
for ibank = 1:nbanks
    
    wf_size = size(bank_models(ibank).filtered_wf);
    nunits = wf_size(1);
    nspikes_per_unit = wf_size(2);
    nchan = wf_size(3);
    ntimes = wf_size(4);
       
    wfMeans = sq(mean(bank_models(ibank).filtered_wf, 2));
    global_mean = sq(mean(mean(wfMeans, 1), 2));    
        
    fprintf('computing AmpScores for bank %d of %d\n', ibank, nbanks);
    
    chMap = readNPY(fullfile(bank_models(ibank).raw_dir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
    nChInMap = numel(chMap);
    
    %calculate amplitudes
    numclus = length(bank_models(ibank).wf_unit_ids);
    
    n = ntimes*nChInMap;
    betweens = zeros(numclus, n);
    
    for i = 1:numclus
        betweens(i, :) = vec(bsxfun(@minus, sq(wfMeans(i, :, :)), global_mean').^2);
    end
    
    withins = zeros(numclus, n);
    for i = 1:numclus
        rs_filtered_wf = reshape(sq(bank_models(ibank).filtered_wf(i, :, :, :)), nspikes_per_unit, n);
        rs_wf_means = reshape(wfMeans(i, :, :), 1, n);
        withins(i, :) = mean(bsxfun(@minus, rs_filtered_wf, rs_wf_means).^2, 1);
    end
    
    betweens_neursum = sum(betweens, 1);
    withins_neursum = sum(withins, 1);
    
    scores = zeros(1, nChInMap);
    for ich = 1:nChInMap
        iselect = ich:nChInMap:n;
        
        scores(ich) = sum(betweens_neursum(iselect)./withins_neursum(iselect));
    end    
    ampscores{ibank} = scores;    
end
fprintf('Ampscore Computation '); toc;

switch ops.selection_method
    case 'cbs'        
        tic
        % initialize bank_assignments        
        switch init_method
            case 'line'
                bank_assignments_init = repmat([1, 0], 1, 384/2);% line                
            case 'checker'
                bank_assignments_init = repmat([1, 0, 0, 1], 1, 384/4);% checker                
            case 'mucbs'
                bank_assignments_init = optimize_mucbs(ampscores, bank_models, ops);%mucbs              
        end                        
   
        % greedy search        
        [bank_assignments, search_output] = greedy_search_cbs(bank_models, bank_assignments_init, ops);                
        
        if do_plot
            
            f = figure;
            
            p = get(f, 'Position');
            set(f, 'Position', [p(1:2), [920, 600]]);
            
            h1 = subplot(2,5,[1, 6]);
            plot_np_channels(ampscores, bank_models);
            xlabel('AmpScores');
        
            h2 = subplot(2,5,[2, 7]);
            plot_bank_assignments_img(bank_assignments_init, bank_models);
            switch init_method
                case 'line'
                    init_method_str = 'Line';
                case 'checker'
                    init_method_str = 'Checker';
                case 'mucbs'
                    init_method_str = '\muCBS';
            end
            xlabel(sprintf('%s selection', init_method_str));
            
            h3 = subplot(2,5,[3, 8]);
            plot_bank_assignments_img(bank_assignments, bank_models);
            xlabel('CBS selection');
            linkaxes([h1, h2, h3]);
            colormap hot;
            
            subplot(2, 5, [4, 5])
            plot(search_output.Jtrace, 'Color', [0.1, 0.8, 0.1]);
            ylabel('selection objective J(\theta)')
            xlabel('iteration');
            box off;
            
            
           if do_classifier
            nunits_vec = [bank_models.nunits];
            prop_correct = [search_output.classif_results.init_prop_correct;
                            search_output.classif_results.prop_correct_record];
                    
            avg_prop_correct = ...
                (prop_correct(:, banks_to_validate)*nunits_vec(banks_to_validate)')...
                ./sum(nunits_vec(banks_to_validate));
            
            subplot(2, 5, [9, 10]);
            plot(100*avg_prop_correct, 'o-', 'Color', [0.1, 0.8, 0.1]);
            xlabel('pass #');
            ylabel('classif. accuracy %')
            box off;
           end
           
           if save_plot
               saveas(f, fullfile(plot_dir, 'cbs_results.png'));
           end
        end
                          
        fprintf('CBS Greedy Search Completed. '); toc
    case 'mucbs'                                
        
        bank_assignments = optimize_mucbs(ampscores, bank_models, ops);   
        if do_plot
            f = figure(1);            
            
            h1 = subplot(1,3,1);
            plot_np_channels(ampscores, bank_models);
            xlabel('AmpScores');

            h2 = subplot(1,3,2);
            plot_bank_assignments_img(bank_assignments, bank_models);
         
            xlabel('\muCBS selection');

            linkaxes([h1, h2]);
            colormap hot;    
            if save_plot
                saveas(f, fullfile(plot_dir, 'mucbs_results.png'));
            end
        end
                
    otherwise
        error('unsupported selection method');       
end

%% converting bank assignment to selection map
% The selection map is a vector of size (1, nchan_all), where the ith
% element is the 1-based electrode assignment for the ith channel.
switch probe_version
    case 'phase1'
        selection_map = zeros(1, length(bank_assignments));
        for ichan = 1:length(bank_assignments)
            selection_map(ichan) = 384*bank_assignments(ichan) + ichan;
        end
end


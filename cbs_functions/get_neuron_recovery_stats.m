function [all_units_recovered, recovery_percentages] = ...
get_neuron_recovery_stats(methods, neuron_yields, ...
    cluster_timing_mat, ops2)

thresh = getOr(ops2, 'thresh', 0.8);
banks_to_analyze = getOr(ops2, 'banks_to_analyze', 0:2);

idense = find(ismember(methods, 'dense'));

all_units_recovered = zeros(size(neuron_yields));
for imethod = 1:length(methods)
        
    for ibank = 1:length(banks_to_analyze)
        
        binned_sp_sparse_cell = {};
        sparse_num_clus = neuron_yields(imethod, ibank);
        nsp_sparse = zeros(1, sparse_num_clus);
        
        ct = cluster_timing_mat{imethod, ibank};
        supra_thresh_clusters = find(cellfun(@(x)~isempty(x), ct));

        for n = 1:neuron_yields(imethod, ibank)
                        
            binned_sp_sparse_cell{n} = ...
                sparse(ct{supra_thresh_clusters(n)}.spt_bin);
            nsp_sparse(n) = sum(binned_sp_sparse_cell{n});
        end
                
        dense_timing = cluster_timing_mat{idense, ibank};                
        nclus_dense = neuron_yields(idense, ibank);
        
        sim_matrix = zeros(neuron_yields(idense, ibank), ...
                        neuron_yields(imethod, ibank));
        
        count_dense_clus_recov = 0;
        dense_clus_recov = [];
        supra_thresh_dense_clusters = find(cellfun(@(x)~isempty(x), dense_timing));
        for m = 1:nclus_dense
            
            binned_sp_dense = sparse(dense_timing{supra_thresh_dense_clusters(m)}.spt_bin);
            nsp_dense = sum(binned_sp_dense);
            
            for n = 1:neuron_yields(imethod, ibank)
                denom = min(nsp_dense, nsp_sparse(n));
                binned_sp_sparse = binned_sp_sparse_cell{n};
                
                sim_matrix(m, n) = sum(binned_sp_dense & binned_sp_sparse)/denom;                
                if sim_matrix(m, n) > thresh                    
                    if ~ismember(m, dense_clus_recov)
                        count_dense_clus_recov = count_dense_clus_recov + 1;
                        dense_clus_recov(count_dense_clus_recov) = m;                    
                    end
                end
                 
            end
        end
        
        all_units_recovered(imethod, ibank) = count_dense_clus_recov;    
    end
end

recovery_percentages = sum(all_units_recovered, 2)./sum(neuron_yields(idense, :), 2)*100;


clear

%%%%%
%%%%% This script shows process of validating electrode selection methods.
%%%%% This script loads full-length subsampled recordings from each
%%%%% selection method. Subsampled recordings are composed of dense
%%%%% recordings whose channels corresponding to disabled electrodes have
%%%%% been zeroed. It is assumed that Kilosort2 has been run on these
%%%%% subsampled recordings prior to running any validations.
%%%%% 

%%%%% 
%%%%% First, we specify all of the selection methods we are comparing. For
%%%%% the purposes of this example script, the data corresponding with
%%%%% these subsampling methods follow an appropriate naming convention. 
%%%%% 
methods = {'dense', 'line', 'checker', 'mucbs', 'cbs'};
data_dir = '../data';
plot_dir = '../results';


%%%%% 
%%%%% Then, we specify options and run the analyses.
%%%%% 
ops2.thresh = 0.8;            % in range (0, 1). Similarity threshold for a 
                         % subsampled neuron to "recover" 
                         % a dense neuron. 0.8 was used in the paper.
ops2.banks_to_analyze = 0:1;  % restricting analysis to a subset of banks that
                         % were determined to be at least partially in the
                         % brain.
ops2.avg_spike_rate_cutoff = 100/120; % restricting analysis to the subset
                                     % of neurons whose average rate
                                     % was above the minimum rate accepted
                                     % during CBS. 100 spikes/2 minutes = 
                                     % 0.8 Hz was used in the paper.                                     
nmethods = length(methods);

% neuron yield: number of neurons discovered by each algorithm
[neuron_yields, cluster_timing_mat, cluster_mat] = ...
    load_neuron_yield_info(data_dir, methods, ops2);

% neuron recovery: number of available neurons that were "recovered" by
% each selection method
[all_units_recovered, recovery_percentages] = ...
    get_neuron_recovery_stats(methods, neuron_yields, ...
    cluster_timing_mat, ops2);

% SNR and Amplitude analysis: subsampled waveform properties like SNR and
% amplitude.
[snrs_by_method, med_amps_by_method] = get_channel_neuron_snrs(data_dir, methods,...
    cluster_timing_mat, cluster_mat, ops2);

% save('did_analyses.mat')

%% Plotting results

colors = {[0, 0, 1], [55, 171, 200]./255, [255, 201, 0]./255, [128, 0, 128]./255, [0 1 0]};     

f = figure(3);
p = get(f, 'Position');
set(f, 'Position', [p(1), p(2), 900, 475]);
bwcolors = [0.1, 0.1, 0.1; 0.7, 0.7, 0.7];
subplot(1, 2, 1);
b = bar([neuron_yields(1, :); zeros(1, size(neuron_yields, 2))], 'stacked');
for ibank = 1:length(ops2.banks_to_analyze)
    b(ibank).FaceColor = bwcolors(ibank, :);
end

ylim([0, max(sum(neuron_yields, 2))]);
ylabel('available neurons');
set(gca, 'xtick', 1);
set(gca, 'xticklabel', {'All'});

subplot(1,2,2);
method_labels = {'Line', 'Checker', '\muCBS', 'CBS'};
b = bar(neuron_yields(2:end, :), 'stacked');
for ibank = 1:length(ops2.banks_to_analyze)
    b(ibank).FaceColor = bwcolors(ibank, :);
end
set(gca, 'Xticklabel', method_labels);
legend('bank 0', 'bank 1', 'Location', 'southeast');
ylim([0, max(sum(neuron_yields, 2))]);
ylabel('neuron yield');

saveas(f, fullfile(plot_dir, 'neuron_yield.png'));

f = figure(4);
b = bar(recovery_percentages(2:end));
set(gca, 'Xticklabel', method_labels);
b.FaceColor = 'flat';
b.CData = cell2mat(colors(2:end)');
ylabel('% of All neurons recovered');

saveas(f, fullfile(plot_dir, 'neuron_recovery.png'));

f = figure(5);
subplot(1,2,1);
for imethod = 1:length(snrs_by_method)
    [fi, xi] = ecdf(snrs_by_method{imethod});    
    plot(xi, fi, 'Color', colors{imethod});
    hold on;    
end
legend(['All', method_labels], 'Location', 'southeast');
hold off;
xlabel('channel SNR (dB)');
ylabel('cumulative prob.');

subplot(1,2,2);
x = cell2mat(med_amps_by_method');
g = [];
for i = 1:nmethods   
    g = [g, i*ones(1, length(med_amps_by_method{i}))];
end

boxplot(x, g, 'Whisker', 0, 'Symbol', '', 'BoxStyle', 'outline', 'Color', cell2mat(colors'), 'Notch', 'on');
ylim([prctile(x, 15), prctile(x, 80)])
set(gca, 'xtick', []);

ylabel('spike amplitude (a.u.)')

saveas(f, fullfile(plot_dir, 'channel_snr_and_spike_ampl.png'));
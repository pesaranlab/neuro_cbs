function channel_info = neuropixel_get_channel_info(rec_path)

%% imec metadata

rec_path_split = split(rec_path, filesep);
rec_dir_path = strjoin(rec_path_split(1:length(rec_path_split) - 1), filesep);

ap_bin_name = char(rec_path_split(length(rec_path_split)));

ap_meta = sglx_util.ReadMeta(ap_bin_name, rec_dir_path);
%lf_meta = sglx_util.ReadMeta(lf_bin_name, imec_rec_dir);

%% getting channel mapping and populating drmap experiment stuff

a = regexp(ap_meta.acqApLfSy, '[0-9]*,[0-9]*,[0-9]*', 'match');
a = strsplit(a{1}, ',');
n_chan_ap = str2double(a{1});

subscripts_all = regexp(ap_meta.snsShankMap, ...
                        '([0-9]*:[0-9]*:[0-9]*:[0-9]*)', 'match');

imro_fields = regexp(ap_meta.imroTbl, ...
                     '([0-9]* [0-9]* [0-9]* [0-9]* [0-9]* [0-9]*)', 'match');

a = cellfun(@(x)strsplit(x, ':'), subscripts_all, 'UniformOutput', false);
b = cellfun(@(x)strsplit(x, ' '), imro_fields, 'UniformOutput', false);

channel_info = struct([]);
for ichan = 1:n_chan_ap
    shank = str2num(a{ichan}{1}) + 1; % make it 1 - based
    row   = str2num(a{ichan}{3}) + 1; 
    col   = str2num(a{ichan}{2}) + 1;        

    bank = str2num(b{ichan}{2});
    
    ref_type = str2num(b{ichan}{3});
    ap_gain = str2num(b{ichan}{4});
    lf_gain = str2num(b{ichan}{5});
    is_hipassed = str2num(b{ichan}{6});
    
    
    channel_info(ichan).shank = shank;
    channel_info(ichan).bank = bank;
    channel_info(ichan).electrode = bank*384 + ichan;   
    channel_info(ichan).row = row;
    channel_info(ichan).col = col;
    channel_info(ichan).ref_type = ref_type;
    channel_info(ichan).ap_gain = ap_gain;
    channel_info(ichan).lf_gain = lf_gain;
    channel_info(ichan).is_hipassed = is_hipassed;
    
end


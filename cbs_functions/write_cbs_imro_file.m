function write_cbs_imro_file(imro_file_path, selection_map, bank_models)

% WRITE_CBS_IMRO_FILE: Writes a selection map to a Neuropixel IMRO file. 
%   inputs:
%       imro_file_path: Absolute or relative path to *.imro file to create/write to.
%       selection_map: selection map to write. This is a vector of length nchan_all, the 
%           total number of recording channels on the probe (including ones excluded 
%           from analysis). The elements are the 1-based electrode indices that each channel
%           is assigned to.

ref_type = bank_models(1).channel_info(1).ref_type;
ap_gain = bank_models(1).channel_info(1).ap_gain;
lf_gain = bank_models(1).channel_info(1).lf_gain;
is_hipassed = bank_models(1).channel_info(1).is_hipassed;

write_imro_table_from_config(imro_file_path, selection_map, ref_type, ap_gain, lf_gain, is_hipassed);


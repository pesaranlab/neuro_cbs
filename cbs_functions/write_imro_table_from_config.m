function write_imro_table_from_config(filename, config, ref_id, ap_gain, lf_gain, ap_apply_hipass)

fid = fopen(filename, 'w');
if nargin < 3
   ref_id = 0;
   ap_gain = 500;
   lf_gain = 250;
   ap_apply_hipass = 1;
end

imroTbl = [];
nchan = 384;
banks = floor(config./nchan);

imroTbl = fprintf(fid, '(0,%d)', nchan);
for ichan = 1:nchan
    fprintf(fid, '(%d %d %d %d %d %d)', ...
            ichan - 1, ...
            banks(ichan), ...
            ref_id, ...
            ap_gain, ...
            lf_gain, ...
            ap_apply_hipass);            
end
fclose(fid);

end
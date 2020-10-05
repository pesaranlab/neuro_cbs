function img = get_img_from_elec_mapped_scores(elec_mapped_scores)

i_img = reshape(1:960,2, 480)';

img = elec_mapped_scores(i_img);
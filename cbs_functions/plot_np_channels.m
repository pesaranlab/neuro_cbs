function h = plot_np_channels(scores, bank_model)

elec_mapped_scores = map_np_bank_to_elec(scores, bank_model);

img = get_img_from_elec_mapped_scores(elec_mapped_scores);

h = tvimage(img');
daspect([1 15 1]);

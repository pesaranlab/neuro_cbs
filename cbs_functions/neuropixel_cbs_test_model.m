function [test_labels_Y, test_ilabels_Y] = neuropixel_cbs_test_model(Xred, class_means, A)
N = size(Xred, 1);
nclasses = size(class_means, 2);
Y = Xred*A;
r = min(100000000, size(Y, 2));

distances_Y = zeros(N, nclasses);
for iunit = 1:nclasses
    class_mean_Y = class_means(:, iunit);    
    eY = sum(bsxfun(@minus, Y(:, 1:r)', class_mean_Y(1:r)).^2, 1);
    distances_Y(:, iunit) = eY;
end

test_labels_Y = bsxfun(@minus, distances_Y, min(distances_Y, [], 2)) == 0;
[~, test_ilabels_Y] = min(distances_Y, [], 2);

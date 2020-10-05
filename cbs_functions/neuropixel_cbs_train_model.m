function [A, V, D, Sol, class_means] = neuropixel_cbs_train_model(Xred, labels, n_classes)
n = size(Xred, 2);
        
[S_Wx, S_bx] = neuropixel_cbs_form_scatter_matrices(Xred, labels, n_classes);

eps = 1e-12;        
S2 = (S_Wx + eps*eye(n));
Sol = S2\S_bx;

m = n_classes - 1;
[V, D] = eigs(Sol, min(size(Sol, 1), m));

m = min(size(Sol, 1), m);
A = V(:, 1:m);  

Y = Xred*A;

class_means = zeros(m, n_classes);


for iunit = 1:n_classes
    ithis = find(labels(:, iunit));
    %                 class_mean_X = mean(Xred(ithis, :), 1)';
    %                  eX = sum(bsxfun(@minus, Xred', class_mean_X).^2, 1);
    %                  distances_X(:, iunit) = eX;
    class_mean_Y = mean(Y(ithis, :), 1)';
   

    class_means(:, iunit) = class_mean_Y;

end


function [S_Wx, S_bx] = neuropixel_cbs_form_scatter_matrices(Xred, labels, n_classes)
n = size(Xred, 2);
N = size(Xred, 1);
label_probs = sum(labels, 1)./N;

S_Wx = zeros(n);
S_bx = zeros(n);
Mu = mean(Xred, 1)';

for iunit = 1:n_classes
    
    Xi = Xred(labels(:,iunit), :);
    Ni = size(Xi, 1);
    Mu_i = mean(Xi, 1)';
    
    centered_Xi = bsxfun(@minus, Xi, Mu_i');
    
    Sigma_i = 1/(Ni-1)*(centered_Xi'*centered_Xi );        
    
    S_Wx = S_Wx + label_probs(iunit)*Sigma_i;
    S_bx = S_bx + label_probs(iunit)*((Mu_i - Mu)*(Mu_i - Mu)');
end    
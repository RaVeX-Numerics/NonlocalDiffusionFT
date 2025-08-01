function Y = uq_ft_rk_wrapper(X, P)
% Wrapper function for FF to be used with UQLab.
% X: Matrix where each row is a sample of uncertain parameters:
%    [mu, alpha2, r, k, sigma, h0]
% P: Structure containing fixed parameters.
% Y: Vector of the Quantity of Interest (QoI) results.
M = P.M; % Number of spatial points factor
T = P.T; % Final time

Ns = size(X, 1); % Number of samples
Y = zeros(Ns, 1);

for i = 1:Ns
    h0_i     = min(max(0.5, X(i, 1)),3);
    mu_i     = min(max(0.1, X(i, 2)),2);
    alpha2_i = min(max(1, X(i, 3)),4);
    r_i      = min(max(0.5, X(i, 4)),2);
    k_i      = min(max(0.25, X(i, 5)),4);
    sigma_i  = min(max(0.1, X(i, 6)),2);

    f_handle = @(u) r_i .* u .* (1 - u / k_i);

    J_handle = @(x) exp(-x.^2 / sigma_i^2) / (sqrt(pi) * sigma_i);
    K_handle = @(z) 0.5 * (1 + erf(z / sigma_i));
    u0_handle = @(x) (h0_i^2 - x.^2)./ h0_i^2*3/(4*h0_i); % to guarantee that P0 = 1


    [X,~,H,G,~] = FF_lax( M,T, mu_i, h0_i, u0_handle, alpha2_i, f_handle, J_handle, K_handle);

    % Y(i) = H(end) - G(end); % Final domain width
    Y(i) = (H(end)-G(end))/(2*h0_i); 



end
end
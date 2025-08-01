function Y = uq_ft_rk_wrapper(X, P)
% Wrapper function for FT_RK to be used with UQLab.
% X: Matrix where each row is a sample of uncertain parameters:
%    [mu, alpha2, r, k, sigma, h0]
% P: Structure containing fixed parameters.
% Y: Vector of the Quantity of Interest (QoI) results.
    M = P.M; % Number of spatial points factor
    T = P.T; % Final time

    Ns = size(X, 1); % Number of samples
    Y = zeros(Ns, 1);

    for i = 1:Ns
        h0_i     = max(0.1, X(i, 1));
        mu_i     = max(0.1, X(i, 2));
        alpha2_i = max(0.1, X(i, 3));
        r_i      = max(0.1, X(i, 4));
        k_i      = max(0.1, X(i, 5));
        sigma_i  = max(0.1, X(i, 6)); 

        f_handle = @(u) r_i .* u .* (1 - u / k_i);

       J_handle = @(x) exp(-x.^2 / sigma_i^2) / (sqrt(pi) * sigma_i);
        K_handle = @(z) 0.5 * (1 + erf(z / sigma_i));
        u0_handle = @(x) (h0_i^2 - x.^2)./ h0_i^2; 


            [X, U, ht, gt, ~] = FT_RK(M, T, mu_i, h0_i, u0_handle, alpha2_i, f_handle, J_handle, K_handle);

            %Y(i) = ht(end);

            % Y(i) = gt(end); % Final position of left boundary g(T)
             Y(i) = ht(end) - gt(end); % Final domain width
      
            %  total_pop = trapz(X,U);
            % Y(i) = total_pop;

       
    end 
end 
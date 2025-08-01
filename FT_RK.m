% Front-tracking method with Runge-Kutta integration for non local Fisher-KPP problem
% Vera Egorova. March, 2025
function [x,u,ht,gt,t] = FT_RK(M,T, mu, h0, u0, alpha2, f, J, K,dt)
% Inputs:
% M: Number of spatial grid points in positive direction (total number of nodes is 2M+1).
% T: Final simulation time.
% mu: Parameter controlling the boundary dynamics.
% h0: Initial length of the domain.
% u0: Function handle for the initial population distribution \( u_0(x) \).
% alpha2: Diffusion coefficient.
% f: Function handle for the growth term \( f(u) \).
% J: Kernel function for nonlocal diffusion.
% K: Cumulative distribution function of the kernel.
% dt(optional): Time step size (default is computed based on stability criteria).
%
% Outputs:
% x: Spatial grid at final time step.
% u: Solution vector
% ht: Upper boundary values.
% gt: Lower boundary values.
% t: Time vector.

dx = h0/M;  % spatial step-size
if nargin < 10
    P0 = integral(u0,-h0,h0);
    dt = 1*dx/(mu*P0);
end
t = 0:dt:T;
N = length(t);
t = linspace(0,T,N);
dt = t(2)-t(1);
x0 = -h0:dx:h0; % initial spatial grid
x = x0;
u = u0(x);
ht = zeros(1, N);  gt = zeros(1, N);
r = zeros(1, N);  l = zeros(1, N);
ht(1) = h0;  gt(1) = -h0;

Start = cputime;
for n = 1:N-1
    delta_l = x(1)-gt(n);
    delta_r = ht(n) - x(end);
    Nx = length(x);
    k = floor((Nx-1)/2);
    Nxi = 2*k+1; % Nxi = Nx, if Nx is odd, otherwise  Nxi = Nx-1
    % if(Nxi<Nx)
    %     fprintf('Nx = %d, Nxi = %d\n',Nx,Nxi)
    % end
    fx = u.*K(x-ht(n));
    J1 = dx/3*( fx(1)+  4*sum(fx(2:2:Nxi-1)) + 2*sum(fx(3:2:Nxi-2)) +fx(Nxi)) +...
        0.5*(Nx-Nxi)*(fx(Nx) + fx(Nxi))*dx+0.5*(fx(1)*delta_l + fx(Nx)*delta_r);
    fx = u.*K(gt(n)-x);
    J2 = dx/3*( fx(1)+  4*sum(fx(2:2:Nxi-1)) + 2*sum(fx(3:2:Nxi-2)) +fx(Nxi)) +...
        0.5*(Nx-Nxi)*(fx(Nx) + fx(Nxi))*dx+0.5*(fx(1)*delta_l + fx(Nx)*delta_r);


    % RK-4:

    fx = J(x' - x) .* u;
    In= dx/3*( fx(:,1)+  4*sum(fx(:,2:2:Nxi-1),2) + 2*sum(fx(:,3:2:Nxi-2),2) +fx(:,Nxi))...
        +0.5*(Nx-Nxi)*(fx(:,Nx) + fx(:,Nxi))*dx+ 0.5*fx(:,1)*delta_l + 0.5*fx(:,Nx)*delta_r;
    k1u = alpha2*(In'-u) + f(u);
    k1g = -mu*J2;
    k1h = mu*J1;
    % ----- k2-----
    u2 = (u+0.5*dt*k1u);
    g2 = gt(n)+0.5*dt*k1g;
    h2 = ht(n)+0.5*dt*k1h;
    fx = J(x' - x) .*u2;

    In= dx/3*(fx(:,1)+  4*sum(fx(:,2:2:Nxi-1),2) + 2*sum(fx(:,3:2:Nxi-2),2) +fx(:,Nxi)) ...
        +0.5*(Nx-Nxi)*(fx(:,Nx) + fx(:,Nxi))*dx+0.5*(fx(:,1)*(x(1)-g2) + fx(:,Nx)*(h2-x(Nx)));
    fx = u2.*K(x-h2);
    J1 = dx/3*( fx(1)+  4*sum(fx(2:2:Nxi-1)) + 2*sum(fx(3:2:Nxi-2)) +fx(Nxi)) +...
        0.5*(Nx-Nxi)*(fx(Nx) + fx(Nxi))*dx+0.5*(fx(1)*(x(1)-g2) + fx(Nx)*(h2-x(Nx)));
    fx = u2.*K(g2-x);
    J2 = dx/3*( fx(1)+  4*sum(fx(2:2:Nxi-1)) + 2*sum(fx(3:2:Nxi-2)) +fx(Nxi)) +...
        0.5*(Nx-Nxi)*(fx(Nx) + fx(Nxi))*dx+0.5*(fx(1)*(x(1)-g2) + fx(Nx)*(h2-x(Nx)));


    k2u = alpha2*(In'-u2)+ f(u2);
    k2g = -mu*J2;
    k2h = mu*J1;

    % ----- k3-----
    u3 = (u+0.5*dt*k2u);
    g3 = gt(n)+0.5*dt*k2g;
    h3 = ht(n)+0.5*dt*k2h;
    fx = J(x' - x) .*u3;

    In= dx/3*(fx(:,1)+  4*sum(fx(:,2:2:Nxi-1),2) + 2*sum(fx(:,3:2:Nxi-2),2) +fx(:,Nxi)) ...
        +0.5*(Nx-Nxi)*(fx(:,Nx) + fx(:,Nxi))*dx+0.5*(fx(:,1)*(x(1)-g3) + fx(:,Nx)*(h3-x(Nx)));
    fx = u3.*K(x-h3);
    J1 = dx/3*( fx(1)+  4*sum(fx(2:2:Nxi-1)) + 2*sum(fx(3:2:Nxi-2)) +fx(Nxi)) +...
        0.5*(Nx-Nxi)*(fx(Nx) + fx(Nxi))*dx+0.5*(fx(1)*(x(1)-g3) + fx(Nx)*(h3-x(Nx)));
    fx = u3.*K(g3-x);
    J2 = dx/3*( fx(1)+  4*sum(fx(2:2:Nxi-1)) + 2*sum(fx(3:2:Nxi-2)) +fx(Nxi)) +...
        0.5*(Nx-Nxi)*(fx(Nx) + fx(Nxi))*dx+0.5*(fx(1)*(x(1)-g3) + fx(Nx)*(h3-x(Nx)));


    k3u = alpha2*(In'-u3)+ f(u3);
    k3g = -mu*J2;
    k3h = mu*J1;

    % ----- k4-----
    u4 = (u+dt*k3u);
    g4 = gt(n)+dt*k3g;
    h4 = ht(n)+dt*k3h;
    fx = J(x' - x) .*u4;

    In= dx/3*(fx(:,1)+  4*sum(fx(:,2:2:Nxi-1),2) + 2*sum(fx(:,3:2:Nxi-2),2) +fx(:,Nxi)) ...
        +0.5*(Nx-Nxi)*(fx(:,Nx) + fx(:,Nxi))*dx+0.5*(fx(:,1)*(x(1)-g4) + fx(:,Nx)*(h4-x(Nx)));
    fx = u4.*K(x-h4);
    J1 = dx/3*( fx(1)+  4*sum(fx(2:2:Nxi-1)) + 2*sum(fx(3:2:Nxi-2)) +fx(Nxi)) +...
        0.5*(Nx-Nxi)*(fx(Nx) + fx(Nxi))*dx+0.5*(fx(1)*(x(1)-g4) + fx(Nx)*(h4-x(Nx)));
    fx = u4.*K(g4-x);
    J2 = dx/3*( fx(1)+  4*sum(fx(2:2:Nxi-1)) + 2*sum(fx(3:2:Nxi-2)) +fx(Nxi)) +...
        0.5*(Nx-Nxi)*(fx(Nx) + fx(Nxi))*dx+0.5*(fx(1)*(x(1)-g4) + fx(Nx)*(h4-x(Nx)));


    k4u = alpha2*(In'-u4)+ f(u4);
    k4g = -mu*J2;
    k4h = mu*J1;
    u = u + dt/6*(k1u+2*k2u+2*k3u+k4u);
    gt(n+1) = gt(n)+dt/6*(k1g+2*k2g+2*k3g+k4g);
    ht(n+1) = ht(n)+dt/6*(k1h+2*k2h+2*k3h+k4h);

    % --- Grid Adaptation (Front Tracking) ---
    delta_l_new = x(1) - gt(n+1);
    delta_r_new = ht(n+1) - x(end);

    num_nodes_add_left = 0;
    if delta_l_new >= dx
        num_nodes_add_left = floor(delta_l_new / dx);
        x_new_left = zeros(1, num_nodes_add_left);
        u_new_left = zeros(1, num_nodes_add_left);

        u_curr = u;
        x_curr = x;
        g_curr = gt(n+1);

        for k = 1:num_nodes_add_left
            x_new_pt = x_curr(1) - dx;
            delta_l_k = x_curr(1) - g_curr;

            if length(u_curr) >= 2
                u_extrap = 2*delta_l_k/(dx+delta_l_k)*u_curr(1) - delta_l_k/(2*dx+delta_l_k)*u_curr(2);
            else
                u_extrap = u_curr(1) * max(0, (delta_l_k - dx)) / delta_l_k;
            end

            % Positivity check
            if u_extrap < 0 || isnan(u_extrap)
                if delta_l_k > 1e-10
                    u_new_val = u_curr(1) * max(0, (delta_l_k - dx)) / delta_l_k;
                else
                    u_new_val = 0; % If boundary is exactly at the point
                end
            else
                u_new_val = u_extrap;
            end

            x_new_left(num_nodes_add_left - k + 1) = x_new_pt;
            u_new_left(num_nodes_add_left - k + 1) = max(0, u_new_val);

            x_curr = [x_new_pt, x_curr];
            u_curr = [max(0, u_new_val), u_curr];
        end

        x = [x_new_left, x];
        u = [u_new_left, u];
    end

    num_nodes_add_right = 0;

    delta_r_new = ht(n+1) - x(end);
    if delta_r_new >= dx
        num_nodes_add_right = floor(delta_r_new / dx);
        r_nodes_added(n) = num_nodes_add_right; % Record

        x_new_right = zeros(1, num_nodes_add_right);
        u_new_right = zeros(1, num_nodes_add_right);

        u_curr = u;
        x_curr = x;
        h_curr = ht(n+1);

        for k = 1:num_nodes_add_right
            x_new_pt = x_curr(end) + dx;
            delta_r_k = h_curr - x_curr(end);

            if length(u_curr) >= 2
                u_extrap = 2*delta_r_k/(dx+delta_r_k)*u_curr(end) - delta_r_k/(2*dx+delta_r_k)*u_curr(end-1);
            else
                u_extrap = u_curr(end) * max(0, (delta_r_k - dx)) / delta_r_k;
            end

            if u_extrap < 0 || isnan(u_extrap)
                if delta_r_k > 1e-10
                    u_new_val = u_curr(end) * max(0, (delta_r_k - dx)) / delta_r_k;
                else
                    u_new_val = 0;
                end
            else
                u_new_val = u_extrap;
            end

            x_new_right(k) = x_new_pt;
            u_new_right(k) = max(0, u_new_val);

            x_curr = [x_curr, x_new_pt];
            u_curr = [u_curr, max(0, u_new_val)];
        end
        x = [x, x_new_right];
        u = [u, u_new_right];
    end
end
x = [gt(n+1), x, ht(n+1)];
u  = [0, u, 0];
Finish = cputime - Start;
fprintf("FT (M = %d) done!, T = %.2f, CPU time: %.2f s\n",M, t(n+1), Finish)
end




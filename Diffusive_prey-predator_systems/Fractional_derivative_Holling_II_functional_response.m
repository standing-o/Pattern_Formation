% fractional-derivative Holling ll functional response
% Caputo's derivative by Grunwald-Letnikov in time and Central differences in space

clear; clc; close all;

% spatial and time discretization
nx = 50; 
h = 1;  % space step
x = linspace(-0.5*h,50+0.5*h,nx+2);
dt = 10^(-10);  % time step
r = dt/h^2;
maxit = 100000;

% parameters which govern equation
alpha = 0.175;   % prey's density at which the predator has the maximum kill rate (0.175)
beta = 0.95;   % maximum birth (0.95)
gam = 0.5;   % dead rate of the predator
d = 20;   % diffusion rate of the predator
eta = 0.8;  sigma = 0.1;
ubar = 0.2; vbar = 0.1;

for pit = 1:1
    rng(pit);
    
    % initial condition
    u = ubar + sigma*(2*rand(nx, maxit)-1);  
    v = vbar + sigma*(2*rand(nx, maxit)-1);  
    
   % no-flux boundary condition 
    u(1, :) = u(2, :);
    u(end, :) = u(end-1, :);
    v(1, :) = v(2, :);
    v(end, :) = v(end-1, :); 
        
    nu = u; nv = v;
    
    % numerical scheme
    for it = 1:maxit   % time loop
        for ix = 2:nx-1   % space loop
            
            % set the source terms
            F = (u(ix, it)*v(ix, it))/(u(ix, it) + alpha);
            f = u(ix, it)*(1-u(ix, it)) - F ;
            g = -gam*v(ix, it) + beta*F ;
            
            nu(ix, it) = u(ix, it) + r*(u(ix-1, it) + u(ix+1, it) - 2*u(ix, it)) + dt*f;
            nv(ix, it) = v(ix, it) + r*d*(v(ix-1, it) + v(ix+1, it) - 2*v(ix, it)) + dt*g;
        end
        
        % reset the variables for next step
        u = nu;
        v = nv;
       
        % visualization
        if mod(it, 1000) == 0
            imagesc(nu(:, it));
%             set(gca, 'xtick', [], 'ytick', []);
            title([num2str(it)])
            colorbar;
            hold on;
            drawnow;
%            print('-djpeg', sprintf('storage path', pit);
        end
    end
 
end
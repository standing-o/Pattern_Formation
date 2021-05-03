clear; clc; close all;

% spatial discretization
xright = 10; nx = 50;
yright = 10; ny = floor(nx*yright/xright);
h = xright/nx;
x = linspace(-0.5*h, xright + 0.5*h, nx+2)';
y = linspace(-0.5*h, yright + 0.5*h, ny+2);

% parameters which govern equation
alpha = 0.175;  % prey's density at which the predator has the maximum kill rate
beta = 0.95;  % maximum birth
gamma = 0.5;  % dead rate of the predator
d = 20;  % diffusion rate of the predator
eta = 0.8;  sigma = 0.1;
ubar = 0.2; vbar = 0.1;


% time discretization
dt = 0.1*h^2;
maxit = 10000;
nn = maxit;

for pit = 1:1
    rng(pit);
    u = ubar + 0.1*(2*rand(nx+2, ny+2)-1);
    v = vbar + 0.1*(2*rand(nx+2, ny+2)-1);
    nu = u; nv = v;
    
    % numerical scheme
    for it = 1:maxit
        % no-flux boundary condition (in progress)
        u(1, :) = u(2, :);
        u(end, :) = u(end-1, :);
        v(1, :) = v(2, :);
        v(end, :) = v(end-1, :);        
        
        % set the source terms
        F = u(2:end-1, 2:end-1).*v(2:end-1, 2:end-1)./(u(2:end-1, 2:end-1) + alpha);
        f = u(2:end-1, 2:end-1).*(1-u(2:end-1, 2:end-1)) - F ;
        g = -gamma*v(2:end-1, 2:end-1) + beta*F ;
        
        % set the source terms
        nu(2:end-1, 2:end-1) = u(2:end-1, 2:end-1) + dt*(lap(u,h) + f);
        nv(2:end-1, 2:end-1) = v(2:end-1, 2:end-1) + dt*(d*lap(u,h) + g);
       
        % reset the variables for next step
        u = nu;
        v = nv;
       
        % visualization
        if mod(it, 100) == 0
            surf(x(2:end-1), y(2:end-1), v(2:end-1, 2:end-1)', 'linestyle', 'none');
            axis image;
            view(2);
            set(gca, 'xtick', [], 'ytick', []);
            title([num2str(it)])
            box on;
            shading interp;
            hold on;
            drawnow;
%            print('-djpeg', sprintf('storage path', pit);
        end
    end
 
end

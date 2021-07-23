% fractional-derivative Holling ll functional response
% Caputo's derivative by Grunwald-Letnikov in time and Central differences in space

clear; clc; close all;

% spatial and time discretization
nx=50; h=1;
dt = 5*10^(-3);  % time step
maxit = 10000;
T = dt*maxit;

x = linspace(-0.5*h,50+0.5*h,nx+2)';
t = linspace(0, T, maxit+1);

% parameters which govern equation
alpha = 0.175;   % prey's density at which the predator has the maximum kill rate (0.175)
beta = 0.95;   % maximum birth (0.95)
gam = 0.5;   % dead rate of the predator
d = 20;   % diffusion rate of the predator (0.1)
sigma = 0.1;
ubar = 0.2; vbar = 0.1;

eta = 0.8;  % fractional order
ddt = dt^eta;
r = eta/ddt - 2/h^2;

uu = zeros(50,maxit/100);

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
            
            % Grunwald-Letnikov coefficients
            c = zeros(1,it); cu = 0; cv = 0;
            if it>2
                for j = 2:it-1
                    c(j) = ((1-(1+eta)/j)^j)*(dt)^(-eta);
                    cu = cu + c(j).*u(ix,it-j);
                    cv = cv + c(j).*v(ix,it-j);
                end
            end

            nu(ix, it) = ddt*(r*u(ix, it) + (1/h^2)*(u(ix-1, it) + u(ix+1, it)) - cu + f);
            nv(ix, it) = ddt*(r*v(ix, it) + (d/h^2)*(v(ix-1, it) + v(ix+1, it)) - cv + g);
            
        end
        
        % reset the variables for next step
        u = nu;
        v = nv;
        
      
%         if mod(it, 100)==0
%             uu(:,it/100)=nu(:,it);
%         end
        
        
        % visualization
        if mod(it, 100) == 0
            surf(t(2:end)', x(2:end-1), u,'linestyle','none');
            xlabel('iterations'), ylabel('X')
            view(2);
            title([num2str(it)]);
            colorbar;
            drawnow;
        end
    end
 
end
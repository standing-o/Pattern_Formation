% first-derivative Holling ll functional response
% FTCS methods

clear; clc; close all;

% spatial and time discretization
nx = 50; 
h = 1;  % space step
x = linspace(-0.5*h,50+0.5*h,nx+2);
dt = 5*10^(-3);  % time step
r = dt/h^2;
maxit = 10^(5);

t = linspace(0, maxit*dt, maxit+1);

% parameters which govern equation
alpha = 0.175;   % prey's density at which the predator has the maximum kill rate (0.175)
beta = 0.95;   % maximum birth (0.95)
gamma = 0.5;   % dead rate of the predator
d = 20;   % diffusion rate of the predator
sigma = 0.1;

ubar = 0.2; vbar = 0.1;

for pit = 1:1
    rng(pit);

    u = zeros(nx,maxit+1); v = zeros(nx,maxit+1);
    
    % initial condition
    u(:,1) = ubar + sigma*(2*rand(nx, 1)-1);  
    v(:,1) = vbar + sigma*(2*rand(nx, 1)-1);  

    % numerical scheme
    for it = 1:maxit  % time loop
%         filename = 'orig.gif';
        
        % no-flux boundary condition
        u(1, it) = u(2, it);
        u(end, it) = u(end-1, it);
        v(1, it) = v(2, it);
        v(end, it) = v(end-1, it);
        
        for ix = 2:nx-1   % space loop
            % set the source terms
            F = (u(ix, it)*v(ix, it))/(u(ix, it) + alpha);
            f = u(ix, it)*(1-u(ix, it)) - F ;
            g = -gamma*v(ix, it) + beta*F ;
            
            u(ix, it+1) = u(ix, it) + r*(u(ix-1, it) + u(ix+1, it) - 2*u(ix, it)) + dt*f;
            v(ix, it+1) = v(ix, it) + r*d*(v(ix-1, it) + v(ix+1, it) - 2*v(ix, it)) + dt*g;
            
        end

%         if it == 1 
%             hh = figure;
%             frame = getframe(hh); 
%             im = frame2im(frame); 
%             [imind,cm] = rgb2ind(im,256); 
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%         end

        if mod(it, 5000) == 0
            surf(t(2:it), x(2:end-1), u(:,2:it), 'linestyle','none');
            view(2);
            xlabel('time'), ylabel('X')
            title([num2str(it)])
            colorbar;
            drawnow;
            print('-djpeg',sprintf('change_L',pit));
            
%             hh = figure;
%             frame = getframe(hh); 
%             im = frame2im(frame); 
%             [imind,cm] = rgb2ind(im,256); 
%             imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        end

    end
 
end

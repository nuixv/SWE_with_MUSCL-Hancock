%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		This code uses The finite volume method for
% 	solving The shallow water equations (SWE) with Topoghaphy.
%				        (1D Hump and Exact solution)
% 	Approximate the numerical flux by Haten-Lax-van Leer contact.
% The MUSCL–Hancock method is adopted to achieve over all second-order accuracy.
% 	The bed slope is estimated using the central-differencing scheme.
%
% 		coded by Narong Batsuwan, Narong.ba@hotmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
global grav;
grav = 9.806;

xmin = 0; xmax = 1000;
n = 41;   %number of cells across box
L = xmax - xmin; %box size (m)
dx = L/n;   %Cell size (m)

xx = xmin-1.5*dx:dx:xmax+1.5*dx;
x = xmin-2*dx:dx:xmax+2*dx; %x coordinate of cell center
nt = 2522; %number of time steps
time = zeros(nt+1,1);

%Initial condition
t = 0;  %start time
h = zeros(n+4,1);
hlx = zeros(n+4,1);
hrx = zeros(n+4,1);
eta = zeros(n+4,1);
u = zeros(n+4,1);
hu = h.*u;  %dicharge for x direction
f1 = zeros(n+4,1);  %fluxes in the $x$ direction
f2 = zeros(n+4,1);  %fluxes in the $x$ direction
hh = zeros(n+4,nt);
z = zeros(n+4,1);
zbx = zeros(n+4,1);

%Source term
for i=1:n+4
    if (x(i)>=125&&x(i)<=875)
        z(i) = 4.75*((1-cos(2*((x(i)-125)/750)*pi))/2);
    end
end

    eta(:) =6;
    %eta(52:104) = 1.0e-06;
    for i=1:n+4
        if (z(i)>=eta(i))
            h(i) = 0;
            eta(i) = z(i);
        else
            h(i) = eta(i) - z(i);
        end
    end

%--------------------------------------------------------------------------
for tstep=1:nt
    %boundary condition
      % x direction
       %left side
        h(1) = h(3);
        h(2) = h(3);
        eta(1) = 10;%eta(3);
        eta(2) = 10;%eta(3);
        hu(1) = hu(3);
        hu(2) = hu(3);
       % right side
        h(n+4) = h(n+2);
        h(n+3) = h(n+2);
        eta(n+4) = eta(n+2);
        eta(n+3) = eta(n+2);
        hu(n+4) = hu(n+2);
        hu(n+3) = hu(n+2);
        
     % return vilocity u,v
       for i = 1:n+4
            if (h(i)<1.0e-06)
               u(i) = 0.0;
            else
               u(i) = hu(i)/h(i);
            end
       end

     % computation of time step
        d = max(abs(u)+sqrt(grav*h));
        dt = (0.5*dx)/(max(max(d)));
        time(tstep+1) = time(tstep) + dt;
        
% computation in x direction  
% solution of Riemann ploblem
for i=3:n+3 %start at interface i-1/2
    % slope limiter
        % left side
            limhl = minmod((h(i-1)-h(i-2))/dx,(h(i)-h(i-1))/dx);
            limetal = minmod((eta(i-1)-eta(i-2))/dx,(eta(i)-eta(i-1))/dx);
            limhul = minmod((hu(i-1)-hu(i-2))/dx,(hu(i)-hu(i-1))/dx);
        % right side
            limhr = minmod((h(i)-h(i-1))/dx,(h(i+1)-h(i))/dx);
            limetar = minmod((eta(i)-eta(i-1))/dx,(eta(i+1)-eta(i))/dx);
            limhur = minmod((hu(i)-hu(i-1))/dx,(hu(i+1)-hu(i))/dx);
    % Data reconstruction
        % left side
            hl = h(i-1)+0.5*dx*limhl ;
            etal = eta(i-1)+0.5*dx*limetal;
            hul = hu(i-1)+0.5*dx*limhul ;
            zl = etal - hl;
        % right side
            hr = h(i)-0.5*dx*limhr;
            etar = eta(i)-0.5*dx*limetar;
            hur = hu(i)-0.5*dx*limhur;
            zr = etar - hr;
%--------------------------------------------------------------------------            
    % Single topography value
            zbx(i) = max(zl,zr);         
    % compute h well-balanced ;
          hl = max(0,etal-zbx(i));
          hr = max(0,etar-zbx(i));
          hlx(i) = hl;
          hrx(i) = hr;
%--------------------------------------------------------------------------            
    % compute vilocity u,v
         % dry bet
           if (hl<1.0e-6)
             ul = 0.0;
           % wet bed
           else
             ul = hul/hl;
           end
           % dry bed
           if (hr<1.0e-6)
             ur = 0.0;
           else
             ur = hur/hr;
           end

    % wave speed
            um = 0.5*(ul+ur)+sqrt(grav*hl)-sqrt(grav*hr);
            hm = (1/grav)*((0.5*(sqrt(grav*hl)+sqrt(grav*hr))+0.25*(ul-ur))^2);
        if (hl==0)
            Sl = ur-2*sqrt(grav*hr);
        elseif (hl>0)
            Sl = min(ul-sqrt(grav*hl),um-sqrt(grav*hm));
        end
        if (hr==0)
            Sr = ul+2*sqrt(grav*hl);
        elseif (hr>0)
            Sr = max(ur+sqrt(grav*hr),um+sqrt(grav*hm));
        end
            Sm = ((Sl*hr*(ur-Sr))-(Sr*hl*(ul-Sl)))/((hr*(ur-Sr))-(hl*(ul-Sl)));

    %compute HLLC flux
       % left side
            f1l = hl*ul;
            f2l = hl*ul^2+0.5*grav*hl^2;
       % right side
            f1r = hr*ur;
            f2r = hr*ur^2+0.5*grav*hr^2;
       % middle side
            f1m = (Sr*f1l-Sl*f1r+Sl*Sr*(hr-hl))/(Sr-Sl);
            f2m = (Sr*f2l-Sl*f2r+Sl*Sr*(hr*ur-hl*ul))/(Sr-Sl);

       % compute F flux
        if (Sl==0)&&(Sr==0)
            f1(i) = 0;
            f2(i) = 0;
        elseif (Sl>=0)
            f1(i) = f1l;
            f2(i) = f2l;
        elseif (Sl<=0)&&(Sm>=0)
            f1(i) = f1m;
            f2(i) = f2m;
        elseif (Sm<=0)&&(Sr>=0)
            f1(i) = f1m;
            f2(i) = f2m;
        else %(Sr<=0)
            f1(i) = f1r;
            f2(i) = f2r;
        end      
end

%Updating of solution
    for i=3:n+2
        %source term
           sox = ((hlx(i+1)+hrx(i))/2)*(((zbx(i+1)-zbx(i))/dx));
           
           h(i) = h(i)-((0.5)*dt/dx)*(f1(i+1)-f1(i));
           hu(i) =hu(i)-((0.5)*dt/dx)*(f2(i+1)-f2(i))-...
               (0.5)*dt*grav*sox; %source
           if (h(i)<0)
               h(i) = 1.0e-06;
               hu(i) = 1.0e-06;
           end
    end
    eta(:) = h(:) + z(:);    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% computation in x direction  
% solution of Riemann ploblem
for i=3:n+3 %start at interface i-1/2
    % slope limiter
        % left side
            limhl = minmod((h(i-1)-h(i-2))/dx,(h(i)-h(i-1))/dx);
            limetal = minmod((eta(i-1)-eta(i-2))/dx,(eta(i)-eta(i-1))/dx);
            limhul = minmod((hu(i-1)-hu(i-2))/dx,(hu(i)-hu(i-1))/dx);
        % right side
            limhr = minmod((h(i)-h(i-1))/dx,(h(i+1)-h(i))/dx);
            limetar = minmod((eta(i)-eta(i-1))/dx,(eta(i+1)-eta(i))/dx);
            limhur = minmod((hu(i)-hu(i-1))/dx,(hu(i+1)-hu(i))/dx);
    % Data reconstruction
        % left side
            hl = h(i-1)+0.5*dx*limhl ;
            etal = eta(i-1)+0.5*dx*limetal;
            hul = hu(i-1)+0.5*dx*limhul ;
            zl = etal - hl;
        % right side
            hr = h(i)-0.5*dx*limhr;
            etar = eta(i)-0.5*dx*limetar;
            hur = hu(i)-0.5*dx*limhur;
            zr = etar - hr;
%--------------------------------------------------------------------------            
    % Single topography value
            zbx(i) = max(zl,zr);         
    % compute h well-balanced ;
          hl = max(0,etal-zbx(i));
          hr = max(0,etar-zbx(i));
          hlx(i) = hl;
          hrx(i) = hr;
%--------------------------------------------------------------------------            
    % compute vilocity u,v
         % dry bet
           if (hl<1.0e-6)
             ul = 0.0;
           % wet bed
           else
             ul = hul/hl;
           end
           % dry bed
           if (hr<1.0e-6)
             ur = 0.0;
           else
             ur = hur/hr;
           end

    % wave speed
            um = 0.5*(ul+ur)+sqrt(grav*hl)-sqrt(grav*hr);
            hm = (1/grav)*((0.5*(sqrt(grav*hl)+sqrt(grav*hr))+0.25*(ul-ur))^2);
        if (hl==0)
            Sl = ur-2*sqrt(grav*hr);
        elseif (hl>0)
            Sl = min(ul-sqrt(grav*hl),um-sqrt(grav*hm));
        end
        if (hr==0)
            Sr = ul+2*sqrt(grav*hl);
        elseif (hr>0)
            Sr = max(ur+sqrt(grav*hr),um+sqrt(grav*hm));
        end
            Sm = ((Sl*hr*(ur-Sr))-(Sr*hl*(ul-Sl)))/((hr*(ur-Sr))-(hl*(ul-Sl)));

    %compute HLLC flux
       % left side
            f1l = hl*ul;
            f2l = hl*ul^2+0.5*grav*hl^2;
       % right side
            f1r = hr*ur;
            f2r = hr*ur^2+0.5*grav*hr^2;
       % middle side
            f1m = (Sr*f1l-Sl*f1r+Sl*Sr*(hr-hl))/(Sr-Sl);
            f2m = (Sr*f2l-Sl*f2r+Sl*Sr*(hr*ur-hl*ul))/(Sr-Sl);

       % compute F flux
        if (Sl==0)&&(Sr==0)
            f1(i) = 0;
            f2(i) = 0;
        elseif (Sl>=0)
            f1(i) = f1l;
            f2(i) = f2l;
        elseif (Sl<=0)&&(Sm>=0)
            f1(i) = f1m;
            f2(i) = f2m;
        elseif (Sm<=0)&&(Sr>=0)
            f1(i) = f1m;
            f2(i) = f2m;
        else %(Sr<=0)
            f1(i) = f1r;
            f2(i) = f2r;
        end      
end

%Updating of solution
    for i=3:n+2
        %source term
           sox = ((hlx(i+1)+hrx(i))/2)*(((zbx(i+1)-zbx(i))/dx));
           
           h(i) = h(i)-((0.5)*dt/dx)*(f1(i+1)-f1(i));
           hu(i) =hu(i)-((0.5)*dt/dx)*(f2(i+1)-f2(i))-...
               (0.5)*dt*grav*sox; %source
           if (h(i)<0)
               h(i) = 1.0e-06;
               hu(i) = 1.0e-06;
           end
    end
    eta(:) = h(:) + z(:);    
%--------------------------------------------------------------------------
    hh(:,tstep) = eta(:);
end
    %set(gca,'nextplot','replacechildren')
    %v=VideoWriter('HLLC_1D_Z_from_H_and_eta.avi');
    %v.FrameRate = 10;
    %open(v);
%for n=1:nt
    plot(xx,hh(:,nt),'-ok');
    hold on;
    plot(xx,z(:),'k');
    hold off;
    axis([xmin xmax 0 12])
    a = floor(time(nt));
    title(['Time = ',num2str(a)])
    xlabel('x(m)')
    ylabel('h+z(m)')
    pause(0.01);
    %frame=getframe(gcf);
    %writeVideo(v,frame);
%end
    %close(v);

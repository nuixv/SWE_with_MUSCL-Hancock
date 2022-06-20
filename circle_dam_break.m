%This code uses a first-order accurate Godunov-type finite volume
%method with HLLc method 
%slop limiter : minmod
%MUSCL-Hancock time integration is used so that the solver is second-order accurate.
%source term with SGM

clear
close all
clc
global grav;
grav = 9.806;

xmin = 0; xmax = 200;
ymin = 0; ymax = 200;
n = 85; %number of cells across box the x-axis
m = 85; %number of cells across box the y-axis
Lx = xmax - xmin; %box size (m)
Ly = ymax - ymin; %box size (m)
dx = Lx/n; %Cell size (m)
dy = Ly/m; %Cell size (m)

xx = xmin-1.5*dx:dx:xmax+1.5*dx;
yy = ymin-1.5*dy:dy:ymax+1.5*dy;
x = xmin-2*dx:dx:xmax+2*dx; %x coordinate of cell center
y = ymin-2*dy:dy:ymax+2*dy; %y coordinate of cell center
nt = 200; %number of time steps
time = zeros(nt+1,1);

%Initial condition
t = 0; %start time
h = zeros(n+4,m+4);
hlx = zeros(n+4,m+4);
hrx = zeros(n+4,m+4);
hly = zeros(n+4,m+4);
hry = zeros(n+4,m+4);
%eta = zeros(n+4,m+4);
u = zeros(n+4,m+4);
v = zeros(n+4,m+4);
hu = h.*u;  %dicharge for x direction
hv = h.*v;  %dicharge for y direction
f1 = zeros(n+3,m+2); %fluxes in the $x$ direction
f2 = zeros(n+3,m+2); %fluxes in the $x$ direction
f3 = zeros(n+3,m+2); %fluxes in the $x$ direction
g1 = zeros(n+2,m+3); %fluxes in the $y$ direction
g2 = zeros(n+2,m+3); %fluxes in the $y$ direction
g3 = zeros(n+2,m+3); %fluxes in the $y$ direction
hh = zeros(n+4,m+4,nt);
z = zeros(n+4,m+4);
zbx = zeros(n+4,m+4);
zby = zeros(n+4,m+4);

for i=1:n+4 % x<=0
    for k=1:n+4
        r = sqrt((x(i)-Lx/2)^2+(y(k)-Ly/2)^2);
        if (r <= 50)
            h(i,k) = 10.0;
        else
            h(i,k) = 0.0;
        end
    end
end
% for i=1:n+4
%     for j=1:m+4
%         if(z(i,j) <= 1)
%             h(i,j) = 1;
%         else
%             h(i,j) = 0;
%         end
%     end
% end
%     for i=1:n+4
%         for j=1:m+4
%             if (z(i,j)>=h(i,j))
%                 h(i,j) = 0.0;
%             end
%         end
%     end
   eta = h + z;
%--------------------------------------------------------------------------
for tstep=1:nt
    %boundary condition
      % x direction
       %left side
        h(1,:) = h(3,:);
        h(2,:) = h(3,:);              
        eta(1,:) = eta(3,:);
        eta(2,:) = eta(3,:);
        hu(1,:) = hu(3,:);
        hu(2,:) = hu(3,:);
        hv(1,:) = hv(3,:);
        hv(2,:) = hv(3,:);              
       % right side
        h(n+4,:) = h(n+2,:);
        h(n+3,:) = h(n+2,:);
        eta(n+4,:) = eta(n+2,:);
        eta(n+3,:) = eta(n+2,:);
        hu(n+4,:) = hu(n+2,:);
        hu(n+3,:) = hu(n+2,:);
        hv(n+4,:) = hv(n+2,:);
        hv(n+3,:) = hv(n+2,:);              
      % y direction
       %bottom side
        h(:,1) = h(:,3);
        h(:,2) = h(:,3);
        eta(:,1) = eta(:,3);
        eta(:,2) = eta(:,3);
        hu(:,1) = hu(:,3);
        hu(:,2) = hu(:,3);
        hv(:,1) = hv(:,3);
        hv(:,2) = hv(:,3);
       % Top side
        h(:,m+4) = h(:,m+2);
        h(:,m+3) = h(:,m+2);
        eta(:,m+4) = eta(:,m+2);
        eta(:,m+3) = eta(:,m+2);
        hu(:,m+4) = hu(:,m+2);
        hu(:,m+3) = hu(:,m+2);
        hv(:,m+4) = hv(:,m+2);
        hv(:,m+3) = hv(:,m+2);
%--------------------------------------------------------------------------      
     % return vilocity u,v
        for i=1:n+4
            for j=1:m+4
                if (h(i,j)<=1.0e-06)
                    u(i,j)=0.0;
                    v(i,j)=0.0;
                else
                    u(i,j)=hu(i,j)/h(i,j);
                    v(i,j)=hv(i,j)/h(i,j);
                end
            end
        end 
%--------------------------------------------------------------------------
    %computation of time step
     dt = 10; %set for big step
        for j=1:m+4
            for i=1:n+4
              if (h(i,j)==0)
                d = 10;
              elseif (h(i,j)>0)
                d = 0.5*min(dx/(abs(u(i,j))+sqrt(grav*h(i,j))),dy/(abs(v(i,j))+sqrt(grav*h(i,j))));
              end
              if (d<dt)
                dt = d;
              end
            end
        end
        time(tstep+1) = time(tstep) + dt;
%--------------------------------------------------------------------------
% computation in x direction  
% solution of Riemann ploblem
for j=3:m+2 %start at interface i-1/2
    for i=3:n+3
    % slope limiter
        % left side
            limhl = minmod((h(i-1,j)-h(i-2,j))/dx,(h(i,j)-h(i-1,j))/dx);
            limetal = minmod((eta(i-1,j)-eta(i-2,j))/dx,(eta(i,j)-eta(i-1,j))/dx);
            limhul = minmod((hu(i-1,j)-hu(i-2,j))/dx,(hu(i,j)-hu(i-1,j))/dx);
            limhvl = minmod((hv(i-1,j)-hv(i-2,j))/dx,(hv(i,j)-hv(i-1,j))/dx);
        % right side
            limhr = minmod((h(i,j)-h(i-1,j))/dx,(h(i+1,j)-h(i,j))/dx);
            limetar = minmod((eta(i,j)-eta(i-1,j))/dx,(eta(i+1,j)-eta(i,j))/dx);
            limhur = minmod((hu(i,j)-hu(i-1,j))/dx,(hu(i+1,j)-hu(i,j))/dx);
            limhvr = minmod((hv(i,j)-hv(i-1,j))/dx,(hv(i+1,j)-hv(i,j))/dx);
    % Data reconstruction
        % left side
            hl = h(i-1,j)+0.5*dx*limhl ;
            etal = eta(i-1,j)+0.5*dx*limetal;
            hul = hu(i-1,j)+0.5*dx*limhul;
            hvl = hv(i-1,j)+0.5*dx*limhvl;
            zl = etal-hl;
        % right side
            hr = h(i,j)-0.5*dx*limhr;
            etar = eta(i,j)-0.5*dx*limetar;
            hur = hu(i,j)-0.5*dx*limhur;
            hvr = hv(i,j)-0.5*dx*limhvr;
            zr = etar-hr;
%--------------------------------------------------------------------------            
    % Single topography value
            zbx(i,j) = max(zl,zr);
    % compute h well-balanced ;
            hl = max(0,etal-zbx(i,j));
            hr = max(0,etar-zbx(i,j));
            hlx(i,j) = hl;
            hrx(i,j) = hr;                   
%--------------------------------------------------------------------------
        % compute vilocity u,v
         % dry bet
            if (hl<1.0e-6)
                ul = 0.0;
                vl = 0.0;
           % wet bed
            else
                ul = hul/hl;
                vl = hvl/hl;
            end
           % dry bed
            if (hr<1.0e-6)
                ur = 0.0;
                vr = 0.0;
            else
                ur = hur/hr;
                vr = hvr/hr;
            end
%--------------------------------------------------------------------------          
    % wave speed
            um = 0.5*(ul+ur)+sqrt(grav*hl)-sqrt(grav*hr);
            hm = (1/grav)*((0.5*(sqrt(grav*hl)+sqrt(grav*hr))+0.25*(ul-ur))^2);
        if (hl>0)&&(hr>0)
            Sl = min(ul-(sqrt(grav*hl)),um-(sqrt(grav*hm)));
            Sr = max(ur+(sqrt(grav*hr)),um+(sqrt(grav*hm)));
            if(hr*(ur-Sr))==(hl*(ul-Sl))
            	Sm = 0;
            else
                Sm = ((Sl*hr*(ur-Sr))-(Sr*hl*(ul-Sl)))/((hr*(ur-Sr))-(hl*(ul-Sl)));
            end              
        elseif (hl==0)&&(hr==0)
            Sl = 0;
            Sm = 0;
            Sr = 0;
        elseif (hl>0)&&(hr==0)
            Sl = ul-sqrt(grav*hl);
            Sr = ul+(2*sqrt(grav*hl));
            Sm = Sr;
        elseif (hl==0)&&(hr>0)
            Sl = ur-(2*sqrt(grav*hr));
            Sr = ur+sqrt(grav*hr);
            Sm = Sl;
        end
%--------------------------------------------------------------------------
    %compute HLLC flux
       % left side
            f1l = hl*ul;
            f2l = hl*ul^2+0.5*grav*hl^2;
            f3l = hl*ul*vl;
       % right side
            f1r = hr*ur;
            f2r = hr*ur^2+0.5*grav*hr^2;
            f3r = hr*ur*vr;
       % middle side
        if(Sr==Sl)
            f1m=0;
            f2m=0;
        else
            f1m = (Sr*f1l-Sl*f1r+Sl*Sr*(hr-hl))/(Sr-Sl);
            f2m = (Sr*f2l-Sl*f2r+Sl*Sr*(hr*ur-hl*ul))/(Sr-Sl);
        end
       % compute F flux
        if (Sl==0)&&(Sr==0)
            f1(i,j) = 0;
            f2(i,j) = 0;
            f3(i,j) = 0;
        elseif (Sl>=0)
            f1(i,j) = f1l;
            f2(i,j) = f2l;
            f3(i,j) = f3l;
        elseif (Sl<=0)&&(Sm>=0)
            f1(i,j) = f1m;
            f2(i,j) = f2m;
            f3(i,j) = f1m*vl;
        elseif (Sm<=0)&&(Sr>=0)
            f1(i,j) = f1m;
            f2(i,j) = f2m;
            f3(i,j) = f1m*vr;
        else %(Sr<=0)
            f1(i,j) = f1r;
            f2(i,j) = f2r;
            f3(i,j) = f3r;
        end 
    end
end
% computation in y direction  
% solution of Riemann ploblem
for i=3:n+2 %start at interface i-1/2
    for j=3:m+3
    % slope limiter
        % left side
            limhl = minmod((h(i,j-1)-h(i,j-2))/dy,(h(i,j)-h(i,j-1))/dy);
            limetal = minmod((eta(i,j-1)-eta(i,j-2))/dy,(eta(i,j)-eta(i,j-1))/dy);
            limhul = minmod((hu(i,j-1)-hu(i,j-2))/dy,(hu(i,j)-hu(i,j-1))/dy);
            limhvl = minmod((hv(i,j-1)-hv(i,j-2))/dy,(hv(i,j)-hv(i,j-1))/dy);
        % right side
            limhr = minmod((h(i,j)-h(i,j-1))/dy,(h(i,j+1)-h(i,j))/dy);
            limetar = minmod((eta(i,j)-eta(i,j-1))/dy,(eta(i,j+1)-eta(i,j))/dy);
            limhur = minmod((hu(i,j)-hu(i,j-1))/dy,(hu(i,j+1)-hu(i,j))/dy);
            limhvr = minmod((hv(i,j)-hv(i,j-1))/dy,(hv(i,j+1)-hv(i,j))/dy);
    % Data reconstruction
        % left side
            hl = h(i,j-1)+0.5*dy*limhl;
            etal = eta(i,j-1)+0.5*dy*limetal;
            hul = hu(i,j-1)+0.5*dy*limhul ;
            hvl = hv(i,j-1)+0.5*dy*limhvl ;
            zl = etal-hl;
        % right side
            hr = h(i,j)-0.5*dy*limhr;
            etar = eta(i,j)-0.5*dy*limetar;
            hur = hu(i,j)-0.5*dy*limhur;
            hvr = hv(i,j)-0.5*dy*limhvr;
            zr = etar-hr;
%--------------------------------------------------------------------------            
    % Single topography value
            zby(i,j)=max(zl,zr);
        % compute h wellbance ;
            hl = max(0,etal-zby(i,j));
            hr = max(0,etar-zby(i,j));
            hly(i,j) = hl;
            hry(i,j) = hr;                       
%--------------------------------------------------------------------------
        % compute vilocity u,v
         % dry bet
            if (hl<1.0e-6)
                 ul = 0.0;
                 vl = 0.0;
           % wet bed
            else
                 ul = hul/hl;
                 vl = hvl/hl;
            end
           % dry bed
            if (hr<1.0e-6)
                 ur = 0.0;
                 vr = 0.0;
            else
                ur = hur/hr;
                vr = hvr/hr;
            end
%--------------------------------------------------------------------------
    % wave speed
            vm=0.5*(vl+vr)+sqrt(grav*hl)-sqrt(grav*hr);
            hm=(1/grav)*((0.5*(sqrt(grav*hl)+sqrt(grav*hr))+0.25*(vl-vr))^2);
        if (hl>0)&&(hr>0)
            Sl = min(vl-(sqrt(grav*hl)),vm-(sqrt(grav*hm)));
            Sr = max(vr+(sqrt(grav*hr)),vm+(sqrt(grav*hm)));
            if (hr*(vr-Sr))==(hl*(vl-Sl))
                Sm = 0;
            else
                Sm = ((Sl*hr*(vr-Sr))-(Sr*hl*(vl-Sl)))/((hr*(vr-Sr))-(hl*(vl-Sl)));
            end
        elseif (hl==0)&&(hr==0)
            Sl = 0;
            Sm = 0;
            Sr = 0;
        elseif (hl>0)&&(hr==0)
            Sl = vl-sqrt(grav*hl);
            Sr = vl+(2*sqrt(grav*hl));
            Sm = Sr;
        elseif (hl==0)&&(hr>0)
            Sl = vr-(2*sqrt(grav*hr));
            Sr = vr+sqrt(grav*hr);
            Sm = Sl;
        end 
       % left side
            g1l = hl*vl;
            g2l = hl*ul*vl;
            g3l = hl*vl^2+0.5*grav*hl^2;
       % right side
            g1r = hr*vr;
            g2r = hr*ur*vr;
            g3r = hr*vr^2+0.5*grav*hr^2;
       % middle side
        if(Sr==Sl)
            g1m = 0;
            g3m = 0;
        else
            g1m = (Sr*g1l-Sl*g1r+Sl*Sr*(hr-hl))/(Sr-Sl);
            g3m = (Sr*g3l-Sl*g3r+Sl*Sr*(hr*vr-hl*vl))/(Sr-Sl);
        end
       % compute G flux
        if (Sl==0)&&(Sr==0)
            g1(i,j) = 0;
            g2(i,j) = 0;
            g3(i,j) = 0;
        elseif (Sl>=0)
            g1(i,j) = g1l;
            g2(i,j) = g2l;
            g3(i,j) = g3l;
        elseif (Sl<=0)&&(Sm>=0)
            g1(i,j) = g1m;
            g2(i,j) = g1m*ul;
            g3(i,j) = g3m;
        elseif (Sm<=0)&&(Sr>=0)
            g1(i,j) = g1m;
            g2(i,j) = g1m*ur;
            g3(i,j) = g3m;
        else %(Sr<=0)
            g1(i,j) = g1r;
            g2(i,j) = g2r;
            g3(i,j) = g3r;
        end  
    end
end
%--------------------------------------------------------------------------
%Updating of solution at dt/2
	for i=3:n+2
        for j=3:m+2
        %source term
            sox = ((hlx(i+1,j)+hrx(i,j))/2)*(((zbx(i+1,j)-zbx(i,j))/dx)); 
            soy = ((hly(i,j+1)+hry(i,j))/2)*(((zby(i,j+1)-zby(i,j))/dy));

            h(i,j) = h(i,j)-(0.5)*(dt/dx)*(f1(i+1,j)-f1(i,j))-(0.5)*(dt/dy)*(g1(i,j+1)-g1(i,j));
            hu(i,j) = hu(i,j)-(0.5)*(dt/dx)*(f2(i+1,j)-f2(i,j))-(0.5)*(dt/dy)*(g2(i,j+1)-g2(i,j))-(0.5)*dt*grav*sox;
            hv(i,j) = hv(i,j)-(0.5)*(dt/dx)*(f3(i+1,j)-f3(i,j))-(0.5)*(dt/dy)*(g3(i,j+1)-g3(i,j))-(0.5)*dt*grav*soy;
        end
	end
%--------------------------------------------------------------------------
    for i=1:n+4
        for j=1:m+4
            if (z(i,j) == max([3-(3/9)*sqrt((x(i)-60)^2+(y(j)-40)^2),...
            3-(3/9)*sqrt((x(i)-60)^2+(y(j)-60)^2),...
            3-(3/9)*sqrt((x(i)-40)^2+(y(j)-50)^2)]))
                if (z(i,j)>=h(i,j))
                    h(i,j) = 0.0;
                end
            end
        end
    end
    eta = h+z;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% computation in x direction  
% solution of Riemann ploblem
for j=3:m+2 %start at interface i-1/2
    for i=3:n+3
    % slope limiter
        % left side
            limhl = minmod((h(i-1,j)-h(i-2,j))/dx,(h(i,j)-h(i-1,j))/dx);
            limetal = minmod((eta(i-1,j)-eta(i-2,j))/dx,(eta(i,j)-eta(i-1,j))/dx);
            limhul = minmod((hu(i-1,j)-hu(i-2,j))/dx,(hu(i,j)-hu(i-1,j))/dx);
            limhvl = minmod((hv(i-1,j)-hv(i-2,j))/dx,(hv(i,j)-hv(i-1,j))/dx);
        % right side
            limhr = minmod((h(i,j)-h(i-1,j))/dx,(h(i+1,j)-h(i,j))/dx);
            limetar = minmod((eta(i,j)-eta(i-1,j))/dx,(eta(i+1,j)-eta(i,j))/dx);
            limhur = minmod((hu(i,j)-hu(i-1,j))/dx,(hu(i+1,j)-hu(i,j))/dx);
            limhvr = minmod((hv(i,j)-hv(i-1,j))/dx,(hv(i+1,j)-hv(i,j))/dx);
    % Data reconstruction
        % left side
            hl = h(i-1,j)+0.5*dx*limhl ;
            etal = eta(i-1,j)+0.5*dx*limetal;
            hul = hu(i-1,j)+0.5*dx*limhul;
            hvl = hv(i-1,j)+0.5*dx*limhvl;
            zl = etal-hl;
        % right side
            hr = h(i,j)-0.5*dx*limhr;
            etar = eta(i,j)-0.5*dx*limetar;
            hur = hu(i,j)-0.5*dx*limhur;
            hvr = hv(i,j)-0.5*dx*limhvr;
            zr = etar-hr;
%--------------------------------------------------------------------------            
    % Single topography value
            zbx(i,j) = max(zl,zr);
    % compute h well-balanced ;
            hl = max(0,etal-zbx(i,j));
            hr = max(0,etar-zbx(i,j));
            hlx(i,j) = hl;
            hrx(i,j) = hr;                   
%--------------------------------------------------------------------------
        % compute vilocity u,v
         % dry bet
            if (hl<1.0e-6)
                ul = 0.0;
                vl = 0.0;
           % wet bed
            else
                ul = hul/hl;
                vl = hvl/hl;
            end
           % dry bed
            if (hr<1.0e-6)
                ur = 0.0;
                vr = 0.0;
            else
                ur = hur/hr;
                vr = hvr/hr;
            end
%--------------------------------------------------------------------------          
    % wave speed
            um = 0.5*(ul+ur)+sqrt(grav*hl)-sqrt(grav*hr);
            hm = (1/grav)*((0.5*(sqrt(grav*hl)+sqrt(grav*hr))+0.25*(ul-ur))^2);
        if (hl>0)&&(hr>0)
            Sl = min(ul-(sqrt(grav*hl)),um-(sqrt(grav*hm)));
            Sr = max(ur+(sqrt(grav*hr)),um+(sqrt(grav*hm)));
            if(hr*(ur-Sr))==(hl*(ul-Sl))
            	Sm = 0;
            else
                Sm = ((Sl*hr*(ur-Sr))-(Sr*hl*(ul-Sl)))/((hr*(ur-Sr))-(hl*(ul-Sl)));
            end              
        elseif (hl==0)&&(hr==0)
            Sl = 0;
            Sm = 0;
            Sr = 0;
        elseif (hl>0)&&(hr==0)
            Sl = ul-sqrt(grav*hl);
            Sr = ul+(2*sqrt(grav*hl));
            Sm = Sr;
        elseif (hl==0)&&(hr>0)
            Sl = ur-(2*sqrt(grav*hr));
            Sr = ur+sqrt(grav*hr);
            Sm = Sl;
        end
%--------------------------------------------------------------------------
    %compute HLLC flux
       % left side
            f1l = hl*ul;
            f2l = hl*ul^2+0.5*grav*hl^2;
            f3l = hl*ul*vl;
       % right side
            f1r = hr*ur;
            f2r = hr*ur^2+0.5*grav*hr^2;
            f3r = hr*ur*vr;
       % middle side
        if(Sr==Sl)
            f1m=0;
            f2m=0;
        else
            f1m = (Sr*f1l-Sl*f1r+Sl*Sr*(hr-hl))/(Sr-Sl);
            f2m = (Sr*f2l-Sl*f2r+Sl*Sr*(hr*ur-hl*ul))/(Sr-Sl);
        end
       % compute F flux
        if (Sl==0)&&(Sr==0)
            f1(i,j) = 0;
            f2(i,j) = 0;
            f3(i,j) = 0;
        elseif (Sl>=0)
            f1(i,j) = f1l;
            f2(i,j) = f2l;
            f3(i,j) = f3l;
        elseif (Sl<=0)&&(Sm>=0)
            f1(i,j) = f1m;
            f2(i,j) = f2m;
            f3(i,j) = f1m*vl;
        elseif (Sm<=0)&&(Sr>=0)
            f1(i,j) = f1m;
            f2(i,j) = f2m;
            f3(i,j) = f1m*vr;
        else %(Sr<=0)
            f1(i,j) = f1r;
            f2(i,j) = f2r;
            f3(i,j) = f3r;
        end 
    end
end
% computation in y direction  
% solution of Riemann ploblem
for i=3:n+2 %start at interface i-1/2
    for j=3:m+3
    % slope limiter
        % left side
            limhl = minmod((h(i,j-1)-h(i,j-2))/dy,(h(i,j)-h(i,j-1))/dy);
            limetal = minmod((eta(i,j-1)-eta(i,j-2))/dy,(eta(i,j)-eta(i,j-1))/dy);
            limhul = minmod((hu(i,j-1)-hu(i,j-2))/dy,(hu(i,j)-hu(i,j-1))/dy);
            limhvl = minmod((hv(i,j-1)-hv(i,j-2))/dy,(hv(i,j)-hv(i,j-1))/dy);
        % right side
            limhr = minmod((h(i,j)-h(i,j-1))/dy,(h(i,j+1)-h(i,j))/dy);
            limetar = minmod((eta(i,j)-eta(i,j-1))/dy,(eta(i,j+1)-eta(i,j))/dy);
            limhur = minmod((hu(i,j)-hu(i,j-1))/dy,(hu(i,j+1)-hu(i,j))/dy);
            limhvr = minmod((hv(i,j)-hv(i,j-1))/dy,(hv(i,j+1)-hv(i,j))/dy);
    % Data reconstruction
        % left side
            hl = h(i,j-1)+0.5*dy*limhl;
            etal = eta(i,j-1)+0.5*dy*limetal;
            hul = hu(i,j-1)+0.5*dy*limhul ;
            hvl = hv(i,j-1)+0.5*dy*limhvl ;
            zl = etal-hl;
        % right side
            hr = h(i,j)-0.5*dy*limhr;
            etar = eta(i,j)-0.5*dy*limetar;
            hur = hu(i,j)-0.5*dy*limhur;
            hvr = hv(i,j)-0.5*dy*limhvr;
            zr = etar-hr;
%--------------------------------------------------------------------------            
    % Single topography value
            zby(i,j)=max(zl,zr);
        % compute h wellbance ;
            hl = max(0,etal-zby(i,j));
            hr = max(0,etar-zby(i,j));
            hly(i,j) = hl;
            hry(i,j) = hr;                       
%--------------------------------------------------------------------------
        % compute vilocity u,v
         % dry bet
            if (hl<1.0e-6)
                 ul = 0.0;
                 vl = 0.0;
           % wet bed
            else
                 ul = hul/hl;
                 vl = hvl/hl;
            end
           % dry bed
            if (hr<1.0e-6)
                 ur = 0.0;
                 vr = 0.0;
            else
                ur = hur/hr;
                vr = hvr/hr;
            end
%--------------------------------------------------------------------------
    % wave speed
            vm=0.5*(vl+vr)+sqrt(grav*hl)-sqrt(grav*hr);
            hm=(1/grav)*((0.5*(sqrt(grav*hl)+sqrt(grav*hr))+0.25*(vl-vr))^2);
        if (hl>0)&&(hr>0)
            Sl = min(vl-(sqrt(grav*hl)),vm-(sqrt(grav*hm)));
            Sr = max(vr+(sqrt(grav*hr)),vm+(sqrt(grav*hm)));
            if (hr*(vr-Sr))==(hl*(vl-Sl))
                Sm = 0;
            else
                Sm = ((Sl*hr*(vr-Sr))-(Sr*hl*(vl-Sl)))/((hr*(vr-Sr))-(hl*(vl-Sl)));
            end
        elseif (hl==0)&&(hr==0)
            Sl = 0;
            Sm = 0;
            Sr = 0;
        elseif (hl>0)&&(hr==0)
            Sl = vl-sqrt(grav*hl);
            Sr = vl+(2*sqrt(grav*hl));
            Sm = Sr;
        elseif (hl==0)&&(hr>0)
            Sl = vr-(2*sqrt(grav*hr));
            Sr = vr+sqrt(grav*hr);
            Sm = Sl;
        end 
       % left side
            g1l = hl*vl;
            g2l = hl*ul*vl;
            g3l = hl*vl^2+0.5*grav*hl^2;
       % right side
            g1r = hr*vr;
            g2r = hr*ur*vr;
            g3r = hr*vr^2+0.5*grav*hr^2;
       % middle side
        if(Sr==Sl)
            g1m = 0;
            g3m = 0;
        else
            g1m = (Sr*g1l-Sl*g1r+Sl*Sr*(hr-hl))/(Sr-Sl);
            g3m = (Sr*g3l-Sl*g3r+Sl*Sr*(hr*vr-hl*vl))/(Sr-Sl);
        end
       % compute G flux
        if (Sl==0)&&(Sr==0)
            g1(i,j) = 0;
            g2(i,j) = 0;
            g3(i,j) = 0;
        elseif (Sl>=0)
            g1(i,j) = g1l;
            g2(i,j) = g2l;
            g3(i,j) = g3l;
        elseif (Sl<=0)&&(Sm>=0)
            g1(i,j) = g1m;
            g2(i,j) = g1m*ul;
            g3(i,j) = g3m;
        elseif (Sm<=0)&&(Sr>=0)
            g1(i,j) = g1m;
            g2(i,j) = g1m*ur;
            g3(i,j) = g3m;
        else %(Sr<=0)
            g1(i,j) = g1r;
            g2(i,j) = g2r;
            g3(i,j) = g3r;
        end  
    end
end
%--------------------------------------------------------------------------
%Updating of solution at n+1
	for i=3:n+2
        for j=3:m+2
        %source term
            sox = ((hlx(i+1,j)+hrx(i,j))/2)*(((zbx(i+1,j)-zbx(i,j))/dx)); 
            soy = ((hly(i,j+1)+hry(i,j))/2)*(((zby(i,j+1)-zby(i,j))/dy));

            h(i,j) = h(i,j)-(0.5)*(dt/dx)*(f1(i+1,j)-f1(i,j))-(0.5)*(dt/dy)*(g1(i,j+1)-g1(i,j));
            hu(i,j) = hu(i,j)-(0.5)*(dt/dx)*(f2(i+1,j)-f2(i,j))-(0.5)*(dt/dy)*(g2(i,j+1)-g2(i,j))-(0.5)*dt*grav*sox;
            hv(i,j) = hv(i,j)-(0.5)*(dt/dx)*(f3(i+1,j)-f3(i,j))-(0.5)*(dt/dy)*(g3(i,j+1)-g3(i,j))-(0.5)*dt*grav*soy;
        end
	end
%--------------------------------------------------------------------------
    for i=1:n+4
        for j=1:m+4
            if (z(i,j) == max([3-(3/9)*sqrt((x(i)-60)^2+(y(j)-40)^2),...
            3-(3/9)*sqrt((x(i)-60)^2+(y(j)-60)^2),...
            3-(3/9)*sqrt((x(i)-40)^2+(y(j)-50)^2)]))
                if (z(i,j)>=h(i,j))
                    h(i,j) = 0.0;
                end
            end
        end
    end
    eta = h + z;
%--------------------------------------------------------------------------
    hh(:,:,tstep) = eta(:,:);
end
    %set(gca,'nextplot','replacechildren')
    %v=VideoWriter('2D_SWE_FVM_HLLC_SGM_1.avi');
    %v.FrameRate = 10;
    %open(v);
for nnn=1:nt
    ss = surf(xx,yy,hh(:,:,nnn)');
    ss.EdgeColor = 'none';
    %ss.FaceColor = 'b';
    hold on
    zz = surf(xx,yy,z(:,:)');
    zz.FaceColor = [0.8500 0.3250 0.0980];
    hold off
    axis([xmin xmax ymin ymax 0 5])
%     xlabel ('x(m)')
%     ylabel ('y(m)')
    zlabel ('h+z(m)')
    title(['Time = ',num2str(time(nnn+1))])
%     view(0,90);
    pause(0.01)
    %frame=getframe(gcf);
    %writeVideo(v,frame);
end
    %close(v)
   
figure
[X,Y] = meshgrid(xx,yy);
contour(X,Y,hh(:,:,nt)')%,'ShowText','on')
axis([xmin xmax ymin ymax 0 20])
title('contour line')
view(0,90)

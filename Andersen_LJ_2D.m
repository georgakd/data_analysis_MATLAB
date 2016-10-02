%%%%%%%%%%%%%%%%%%%%%    MD program for 2D (Velocity Verlet with Andersen thermostat + LJ interactions)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SYSTEM = SAMPLE
% GEOMETRY = FLAT

% KINDS OF FORCES: LJ POTENTIAL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ASSIGNEMENT OF GLOBAL VARIABLES

%global npartx nparty npart Tref dt a rc rr tsim N d

%%%%%%%%%%%%%%%% subroutine init: Initialization of MD program%%%%%%%%%%%%%%

% We assume m=1, k=1, sigma=1, epsilon=1 (REDUCED UNITS)

% length sigma
% energy epsilon
% mass m
% all other variables are expressed with combination of the above units

% VARIABLES
clear all;
npartx=4;
nparty=4;
npart=npartx*nparty; % Total number of particles
Tref=0.3;
dt=0.01; % step 
tsim=20; % simulation time 
N=tsim/dt; % number of steps
a=2; % crystal parameter
rr=0.5; % radius of spherical particle
rc=2.5*5;   % rc >= 2.5*sigma  cutoff
rcsq=rc^2; 
d=2;  % dimensions
dr=0.05; %RDF parameter
nbath=N;


% SET PARAMETERS TO ZERO-equilibrium at t=0 or n=1
sumvx=0;sumvy=0;vx=zeros();vy=zeros();x=zeros();y=zeros();
sumvx2=0;sumvy2=0;ax=zeros();ay=zeros();fx=zeros();fy=zeros();
K=zeros();U=zeros();V=zeros();E=zeros();T=zeros();temp=zeros();

% Build a square lattice - Draw a box and put the atoms inside it
i=0;

for col=1:nparty
        for row=1:npartx
   
        i=i+1;
        x(1,i)=a*(row-1);    % FIRST STEPS - EQUILIBRIUM 
        y(1,i)=a*(col-1);
    
        end
end
    
xmin = - a/2; xmax = (npartx - 1)*a + a/2; % Dimensions of the box
ymin = - a/2; ymax = (nparty - 1)*a + a/2;
Lx = xmax - xmin;
Ly = ymax - ymin;    

%Graphics of equilibrium state
figure(1);
clf;
theta=(0:0.1:2*pi);
n=1;
clf;
hold on;    
for k=1:npart
    xcircle=rr*cos(theta)+x(n,k);
    ycircle=rr*sin(theta)+y(n,k);
    fill(xcircle,ycircle,'k')
    axis ([xmin xmax ymin ymax]);
    axis square;
    rectangle('Position',[xmin,ymin,Lx,Ly])
    title(n);
end
hold off;

%%
%%%%%%%%%%%%Rescaling velocities %%%%%%%%%%%%%
    
        rand('state',6789);
    sumkin=0;n=1;
   for k=1:npart
        vx(n,k)=rand-0.5;vy(n,k)=rand-0.5;
    
        sumvx=sumvx+vx(n,k);sumvy=sumvy+vy(n,k);  %velocity center of mass
        sumvx2=sumvx2+vx(n,k)^2; sumvy2=sumvy2+vy(n,k)^2;  %kinetic energy
   end

        sumvx=sumvx/npart; sumvy=sumvy/npart;   %   <v>
        sumvx2=sumvx2/npart; sumvy2=sumvy2/npart;   %  <v^2>

        fsx = sqrt((d*Tref)/sumvx2);  %scaling factor due to fixed Temp simulation
        fsy = sqrt((d*Tref)/sumvy2);  %scaling factor due to fixed Temp simulation

        
   for  k=1:npart
        vx(n,k)=(vx(n,k)-sumvx)*fsx; %rescaling velocities
        vy(n,k)=(vy(n,k)-sumvy)*fsy;

        sumkin=sumkin+vx(n,k)^2+vy(n,k)^2;

   end
   
   
U(n)=0;K(n)=0.5*sumkin;E(n)=K(n)+U(n);T(n)=Tref;
%%
%%%%%%%%%%%%Bootstrap step at n=2%%%%%%%%%%%%%%%%%
for k=1:npart
    ax(1,k)=0; 
    ay(1,k)=0;
end

n=2;

for k=1:npart
    x(n,k)=x(n-1,k)+dt*vx(n-1,k); % Euler Method
    y(n,k)=y(n-1,k)+dt*vy(n-1,k);
    if  x(n,k)>xmax           %restrain the atoms in the box
        x(n,k)=x(n,k)-Lx;
    end
    if x(n,k)<xmin
       x(n,k)=x(n,k)+Lx;
    end
    
    if  y(n,k)>ymax            %restrain the atoms in the box
        y(n,k)=y(n,k)-Ly;
    end
    if  y(n,k)<ymin
        y(n,k)=y(n,k)+Ly;
    end
end
%%
%Graphics of n=2
figure(n);
clf;
theta=(0:0.5:2*pi);
for n=1:2;
clf;
hold on;    
for k=1:npart
    xcircle=rr*cos(theta)+x(n,k);
    ycircle=rr*sin(theta)+y(n,k);
    fill(xcircle,ycircle,'k')
    axis ([xmin xmax ymin ymax]);
    axis square;
    rectangle('Position',[xmin,ymin,Lx,Ly])
    title(n);
end
hold off;
pause(n);
end

%%
%%%%%%%%%%%%%%%%Force calculation at n=2%%%%%%%%%%%%
sumV=0;mysum=0;
     for i=1:npart-1
         for j=i+1:npart
                dx=x(n,j)-x(n,i);
                dy=y(n,j)-y(n,i);
                xforce=x(n,j);
                yforce=y(n,j);
%%%%%%%%%%%% use of PBC according to minimum image convention
            if dx>0.5*Lx
                xforce=x(n,j)-Lx;
            end
            if dx<-0.5*Lx
                xforce=x(n,j)+Lx;
            end
            if dy>0.5*Ly
                yforce=y(n,j)-Ly;
            end
            if dy<-0.5*Ly;
                yforce=y(n,j)+Ly;
            end
                      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
             r2sq=(xforce-x(n,i))^2+(yforce-y(n,i))^2;
            if r2sq<rcsq
             fx(i,j)=-4*((-6*(xforce-x(n,i)))/(sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2)^8)+(12*(xforce-x(n,i))/(sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2)^14)));
             fy(i,j)=-4*((-6*(yforce-y(n,i)))/(sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2)^8)+(12*(yforce-y(n,i))/(sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2)^14)));
             fx(j,i)=-fx(i,j);
             fy(j,i)=-fy(i,j);
             V(i,j) = 4*(-1/sqrt((xforce - x(n, i))^2 + (yforce - y(n, i))^2)^(6) + 1/sqrt((xforce - x(n, i))^2 + (yforce - y(n, i))^2)^(12));      
             sumV=sumV+V(i,j);
              
                
            elseif r2sq>rcsq
                fx(i,j)=0;fy(i,j)=0;V(i,j)=0;
            end
         end
     end
    U(n)=sumV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%% Calculation of acceleration and update velocity of bootstrap step%%%%%%%%%%%%%%
sumkin=0;tempa=0;
for k=1:npart
    sumx=0;
    for l=1:k-1
        sumx=sumx+fx(k,l);    
    end
    for l=k+1:npart
        sumx=sumx+fx(k,l);    
    end
    sumy=0;
    for l=1:k-1
        sumy=sumy+fy(k,l);    
    end
    for l=k+1:npart
        sumy=sumy+fy(k,l);    
    end
    ax(n,k)=sumx;
    ay(n,k)=sumy;
    vx(n,k)=vx(n-1,k)+ax(n,k)*dt;
    vy(n,k)=vy(n-1,k)+ay(n,k)*dt;
    temp(n,k)=(vx(n,k)^2+vy(n,k)^2)/(d*npart);
    
    sumkin=sumkin+vx(n,k)^2+vy(n,k)^2;
    tempa=tempa+temp(n,k);
end
K(n)=0.5*sumkin;E(n)=K(n)+U(n);T(n)=tempa;
%%      
 %%%%%%%%%%%%%%%%%%%%%% LOOP OVER TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     nu=10; %collision parameter
     Tnew=0.1;sig=sqrt(Tnew);  % because sigma=sqrt(T);   
     counter=0;
    
     
for n=3:N
    
    
    sumkin=0;sumV=0;tempa=0;  %for the calculation of energies
    sumvx=0;sumvy=0;
    for k=1:npart     
            x(n,k)=x(n-1,k)+dt*vx(n-1,k)+0.5*ax(n-1,k)*dt^2; % first part of velocity verlet algorithm 
            y(n,k)=y(n-1,k)+dt*vy(n-1,k)+0.5*ay(n-1,k)*dt^2;
           
        
         if  x(n,k)>xmax           %restrain the atoms in the box
            x(n,k)=x(n,k)-Lx;
         end
        if  x(n,k)<xmin
            x(n,k)=x(n,k)+Lx;
        end
    
    
        if  y(n,k)>ymax           %restrain the atoms in the box
            y(n,k)=y(n,k)-Ly;
         end
        if  y(n,k)<ymin
            y(n,k)=y(n,k)+Ly;
        end
      
    end 
              
 for i=1:npart-1              % update forces
         for j=i+1:npart
                
                dx=x(n,j)-x(n,i);
                dy=y(n,j)-y(n,i);
                
                xforce=x(n,j);
                yforce=y(n,j);
%%%%%%%%%%%% use of PBC according to minimum image convention
            if dx>0.5*Lx
                xforce=x(n,j)-Lx;
            end
            if dx<-0.5*Lx
                xforce=x(n,j)+Lx;
            end
            if dy>0.5*Ly
                yforce=y(n,j)-Ly;
            end
            if dy<-0.5*Ly;
                yforce=y(n,j)+Ly;
            end
                      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                r2sq=(xforce-x(n,i))^2+(yforce-y(n,i))^2; %improving speed
               if r2sq<rcsq
                fx(i,j)=-4*((-6*(xforce-x(n,i)))/(sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2)^8)+(12*(xforce-x(n,i))/(sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2)^14)));
                fy(i,j)=-4*((-6*(yforce-y(n,i)))/(sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2)^8)+(12*(yforce-y(n,i))/(sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2)^14)));
                fx(j,i)=-fx(i,j);
                fy(j,i)=-fy(i,j);
                
                V(i,j) = 4*(-1/sqrt((xforce - x(n, i))^2 + (yforce - y(n, i))^2)^(6) + 1/sqrt((xforce - x(n, i))^2 + (yforce - y(n, i))^2)^(12));
                sumV=sumV+V(i,j);
                 r2(i,j)=sqrt((xforce-x(n,i))^2+(yforce-y(n,i))^2);
       
               elseif r2sq>rcsq
                   fx(i,j)=0;fy(i,j)=0;V(i,j)=0;
               end
               
         end
 end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for k=1:npart       % update accelerations
    sumx=0;
    for l=1:k-1
        sumx=sumx+fx(k,l);    
    end
    for l=k+1:npart
        sumx=sumx+fx(k,l);    
    end
    sumy=0;
    for l=1:k-1
        sumy=sumy+fy(k,l);    
    end
    for l=k+1:npart
        sumy=sumy+fy(k,l);    
    end
    
    ax(n,k)=sumx;
    ay(n,k)=sumy;
    vx(n,k)=vx(n-1,k)+0.5*dt*(ax(n,k)+ax(n-1,k));  %% second part of velocity verlet algorithm 
    vy(n,k)=vy(n-1,k)+0.5*dt*(ay(n,k)+ay(n-1,k));
    temp(n,k)=(vx(n,k)^2+vy(n,k)^2)/(d*npart);
 
    number=rand;  % Andersen thermostat
    if n>nbath
    if number < (nu*dt) % nu=collision rate
        counter=counter+1;
        
        v1=rand;
        v2=rand;
        
        l1=sqrt(-2*log(v1))*cos(2*pi*v2);l2=sqrt(-2*log(v1))*sin(2*pi*v2); %Box-Muller mehod
  % for producing gaussian numbers from 2 different uniform distributions
        vx(n,k)=sig*l1;vy(n,k)=sig*l2;
        temp(n,k)=(vx(n,k)^2+vy(n,k)^2)/(d*npart);
        
   
    end
    end
        sumkin=sumkin+vx(n,k)^2+vy(n,k)^2; %kinetic energy 
        tempa=tempa+temp(n,k);
    
   end 
        
    U(n)=sumV;
    K(n)=0.5*sumkin;
    E(n)=U(n)+K(n);
    T(n)=tempa; 
        
    
end % end of simulation

%%
%%% RDF%%%
h=zeros();g=zeros();r1=zeros();

for p=1:100
    r1(p) = (p - 0.5)*dr; mysum = 0;
  for i=1:npart-1              
         for j=i+1:npart
     if (r2(i, j) < p*dr && r2(i, j) >= (p - 1)*dr)
        mysum = mysum + 1;
     end
         end
  end
        h(p) = mysum;
        g(p) = (Lx*Ly*h(p))/(pi*npart^2*r1(p)*dr);
end

figure(3)
plot(r1(1:p),g(1:p))

%%
%%%%%Plot histograms%%%%%%%%%%%%%%%%
varray=sqrt(vx(N,:).^2+vy(N,:).^2);
Px(N,:)=(sqrt(1./(2*pi.*Tnew))).*exp(-(0.5*vx(N,:).^2)/Tnew); %distribution for components
Py(N,:)=(sqrt(1./(2*pi.*Tnew))).*exp(-(0.5*vy(N,:).^2)/Tnew); %distribution for components
Psquared=(sqrt((2/pi)*Tnew^-3)).*varray.^2.*exp(-(0.5*varray.^2)/Tnew); %distr for velocity vector

figure(4)
subplot(3,1,1)
hist(vx(N,:))
title('vx last time step of all particles')
subplot(3,1,2)
hist(vy(N,:))
title('vy last time step of all particles')
subplot(3,1,3)
hist(varray)
title('v last time step of all particles')

figure(5)
subplot(3,1,1)
scatter(vx(N,:),Px(N,:),'b')
title('vx last time step of all particles')
subplot(3,1,2)
scatter(vy(N,:),Py(N,:),'r')
title('vy last time step of all particles')
subplot(3,1,3)
scatter(varray,Psquared,'k')
title('v last time step of all particles')

%%
%%%%%%%Plot Energies%%%%%%%%%%%%%%%
figure (6)
plot(1:N,K(1:N),'r',1:N,U(1:N),'b',1:N,E(1:N),'k')
xlabel('Time Steps');
ylabel('E(t)');
title('Total Energy over time')
%%
%%%%%%%Plot Temperature%%%%%%%%%%%%%%%
figure (7)
plot(1:N,T(1:N),'r')
xlabel('Time Steps');
ylabel('T(t)');
title('Temperature over time')
%%
%{
%Graphics of system evolution
figure(1);
clf;
theta=(0:0.5:2*pi);
for n=1:N;
clf;
hold on;    
for k=1:npart
    xcircle=rr*cos(theta)+x(n,k);
    ycircle=rr*sin(theta)+y(n,k);
    fill(xcircle,ycircle,'k')
    axis ([xmin xmax ymin ymax]);
    axis square;
    rectangle('Position',[xmin,ymin,Lx,Ly])
    title(n);
end
hold off;
pause(0.001);
end
%}

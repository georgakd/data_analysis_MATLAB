% Simple MD simulation
% System with many identical atoms that interact with each other in a
% spherically symmetric potential. Will they bond to form a solid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Define constants of the problem
r0=2; % distance at which the force is zero (in A)
alpha=1; %amount proportional to the force
N=64; % number of atoms
L=10; % size of the box where the atoms move (in A)
m=1; %attractive term
A=1; %constant of attractive term
n=8; %repulsive term
B=(m*A*r0^(-m+n))/n;
 

rand('seed',1);

% Generate random initial atom positions

fid = fopen('InitialPositions.txt', 'a');


for i=1:N
    x(i)= (L*rand);
    y(i)= (L*rand);
    fprintf(fid,'%f %f\n',x(i),y(i));
end

fclose(fid)
%%
figure (1)
scatter(x,y,'r','filled')
xlabel('x(Angstrom)')
ylabel('y(Angstrom)')
hold on
%%
for k=1:1000 % counter
    for i=2:N % center of gravity is fixed upon atom 1
    Fxi=0; Fyi=0; % initialize cumulative sums
    for j=1:N
        if j~=i

    rij=sqrt((x(j)-x(i))*(x(j)-x(i))+(y(j)-y(i))*(y(j)-y(i))); % Calc atomic seperations
    Fij=interatomicforce(rij);
 Fxi=Fxi+Fij*(x(j)-x(i))/rij;
 Fyi=Fyi+Fij*(y(j)-y(i))/rij;
        end
    end
    Fnorm=interatomicforce(2*r0);
    dx = (alpha*Fxi)/Fnorm; % divide by Fij(2r0) for normalisation
    dy = (alpha*Fyi)/Fnorm;
    if dx<(-L/100)
        dx=-L/100;
    end
    if dx>L/100
        dx=L/100;
    end
    if dy<(-L/100)
        dy=-L/100;
    end
    if dy>L/100
        dy=L/100;
    end
     x(i)=x(i)+dx;
     y(i)=y(i)+dy;
    
    
    end
end
%%

scatter(x,y,'b','filled')
xlabel('x(Angstrom)')
ylabel('y(Angstrom)')
hold off

% Driver script for solving the 1D Burgers equations using linear basis
% functions.
% Both strong form and weak form are implemented.
% Option=1: strong form, 2: weakform. 
clc; clear; clf;
global N J nx K M lmln 
global vmapM vmapP vmapI vmapO


N =50; %Number of elements
xL=-1; xR=1;
X = linspace(xL, xR, N+1);
FinalTime =0.5;

% We use linear basis function.
% We compute the mass matrix on the standard/reference element [0,1]; 
% All elements are mapped to this standard element;
syms t;
phi=[(1-t)/2, (1+t)/2];
phip=[diff(phi(1),1);diff(phi(2),1)];

J =((xR-xL)/N/2);
M=zeros(2,2);
for i=1:2
    for j=1:2
   M(i,j)=int(phi(i)*phi(j)*J,t,-1,1);     
    end
end

K=zeros(2,2);
for i=1:2
    for j=1:2
   K(i,j)=int(phi(i)*phip(j),t,-1,1);     
    end
end

r=[-1;1];
x = ones(2,1)*X(1:N) + 0.5*(r+1)*(X(2:N+1)-X(1:N));%find the x coordinates for all the elements

%Flux lifting term, i.e., the boundary term in weak formulation
lmln = zeros(2,2); 
lmln(1,1) = 1.0; lmln(2,2) = 1.0;

% Surface normals for each element, direction points to the right are positive
nx = zeros(2, N); nx(1, :) = -1.0; nx(2, :) = 1.0;

%Connectivity
vmapM   = zeros(2*N,1);  vmapP   = zeros(2*N,1); 
for i=1:1:N
vmapM(i*2-1,1)=(i-1)*2+1;
vmapM(i*2,1)=(i-1)*2+2;
end
for i=1:1:N
vmapP(i*2-1,1)=(i-1)*2;
vmapP(i*2,1)=(i-1)*2+2+1;
end
vmapP(1,1)=1; vmapP(2*N,1)=N*2;
vmapI = 1; vmapO = 2*N;% Create specific left (inflow) and right (outflow) maps

% Initial condition 
u = sin(pi*x);
figure(1);plot(x,u,'r');hold on;

time = 0;
% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.2; umax = max(max(abs(u)));
dt = CFL* min(xmin/umax);
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 



for tstep=1:Nsteps
    rhs = BurgersRHS(u); 
    u = u+dt*rhs; %Forward Euler to update u^{n+1} from u^n
    time = time+dt;
    %if (mod(tstep,10))==1 %Plotting every 10th time step
       plot(x,u,'b-');%hold off;
       axis([-1 1 -1.5 1.5]);
       title('Solution at different times');
       xlabel('x axis');
       ylabel('u axis');
    %pause
    %end
end

function val = BurgersFlux(u)
    val = 0.5*u^2;
end

function val = LaxFriedrichsFlux(uL,uR)
    C = max(abs(uL),abs(uR));
    val = 0.5*(BurgersFlux(uL)+BurgersFlux(uR)) + 0.5*C*(uL-uR);
end

function rhs = FluxRHS(u)
global N
%     MaxVel = max(max(abs(u)));% Wave speed=constant C in Lax-Friedrichs flux
    
     rhs = zeros(2,N);
     rhs(1,1) = +LaxFriedrichsFlux(u(1,1),u(1,1));
     rhs(2,1) = -LaxFriedrichsFlux(u(2,1),u(1,2));
     
     for el = 2:N-1
        rhs(1,el) = +LaxFriedrichsFlux(u(2,el-1),u(1,el));
        rhs(2,el) = -LaxFriedrichsFlux(u(2,el),u(1,el+1));
     end
    
     rhs(1,N) = +LaxFriedrichsFlux(u(2,N-1),u(1,N));
     rhs(2,N) = -LaxFriedrichsFlux(u(2,N),u(2,N));
end

function [rhs] = BurgersRHS(u)
global N nx K M lmln 
global vmapM vmapP vmapI vmapO

f = arrayfun(@BurgersFlux,u);
Fl = FluxRHS(u);

rhs = M\(K'*f + Fl);

end




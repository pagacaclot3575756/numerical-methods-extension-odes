function [t,y] = rkvectorial(f, ta, y0, m, n)
% rkvectorial - Classical explicit Runge-Kutta method (RK4), vectorized
%
% Inputs:
%   f  - function handle, f(t,y), defining the ODE system
%   ta - [t0, T], time interval
%   y0 - initial condition vector
%   m  - number of equations of the system
%   n  - number of subintervals (mesh size)
%
% Output:
%   t - time mesh
%   y - approximate solution at each mesh point
%
% Used as the starting procedure for the implicit and explicit
% multistep methods, providing the initial values needed before 
% the multistep recursion begins.

t0=ta(1);
T=ta(2);

h=(T-t0)/n;

t=zeros(n+1,1);
y=zeros(m,n+1);

y(:,1)=y0;
t(1)=t0;

for i=1:n
    t(i+1)=h+t(i); 

    z1=y(:,i);
    s1=t(i);
    k1=f(s1,z1);
    
    z2=z1+(h/2)*k1;
    s2=s1+h/2;
    k2=f(s2,z2);
    
    z3=z1+(h/2)*k2;
    s3=s1+h/2;
    k3=f(s3,z3);
    
    z4=z1+h*k3;
    s4=s1+h;
    k4=f(s4,z4);
    y(:,i+1)=y(:,i)+ (h/6)*(k1+2*k2+2*k3+k4);
end
end
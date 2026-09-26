function [t,y]=Nystrom1(f,ta,y0,w,n)  %Explícito 
% Nystrom1 - Explicit Nystrom-type linear multistep method
%
% Inputs:
%   f  - function handle, f(t,y), defining the ODE system
%   ta - [t0, T], time interval
%   y0 - initial condition vector
%   w  - number of steps of the multistep method
%   n  - number of subintervals (mesh size)
%
% Output:
%   t - time mesh
%   y - approximate solution at each mesh point
%
% Uses an explicit RK method (rkvectorial) to compute the first
% n_arr+1 starting values, then advances explicitly using the
% two-step relation y(i+1) = y(i-1) + s*sum(b_l * f(t(i-l), y(i-l))).
t0=ta(1);
T=ta(2);

s=(T-t0)/n;

d=length(y0); %Numero de ecuaciones del sistema
y0=y0(:);

t=linspace(t0,T,n+1);
y=zeros(d,n+1);

b=zeros(w+1,1);

m_ent=(0:w)';
for l=0:w
    m=m_ent; m(l+1)=[]; %Eliminamos la componente l
    if isempty(m)
        integrando=@(v) ones(size(v));
    else
        integrando=@(v) arrayfun(@(x) prod((x+m)./(m-l)), v); 
        %Para que la función se aplique al vector v
    end
b(l+1)=integral(integrando,-1,1); %Porque d=1, u=-1
end

%Necesitamos y_0,...,y_{w+1} (w+2 puntos), con i=w+1 necesitamos y(:,w)
%(w+1 puntos), pero si w=0 i empieza en 1 y necesitamos y(:,0), tenemos que
%arrancar con w+2 puntos
n_arr=max(w+1,2); %n_arr-w>=1 <--> n_arr>=w+1  Al menos 2 para cubrir y(:,i-1) e y(:,i-w)
ts=[t0, t0+ n_arr*s]; 
[~,ys] = rkvectorial(f, ts, y0, d, n_arr);

y(:,1:n_arr +1)=ys; %Ponemos el arranque que hemos calculado antes, tenemos ya 
               %los primeros w+2 valores
%Ahora ya podemos calcular de y_{w+2} en adelante ya que ya podemos
%calcular f_{w+1},f_w,f_{w-1},...,f_0

for i=(n_arr):n
    Suma=0;
    for l=0:w
        Suma=Suma + b(l+1)*f(t(i-l),y(:,i-l));
    end
    y(:,i+1)=y(:,i-1)+s*Suma;
end
end



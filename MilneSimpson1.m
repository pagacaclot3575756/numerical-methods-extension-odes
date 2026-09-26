function [t,y]=MilneSimpson1(f,ta,y0,w,n)  %Implícito 
% MilneSimpson1 - Implicit linear multistep method (predictor-corrector)
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
% n_arr+1 starting values, then applies fixed-point iteration
% (maxit steps) to solve the implicit equation at each step.

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
b(l+1)=integral(integrando,-2,0); %Porque d=0, u=-2
end


%El bucle empieza en i=n_arr +1 y necesita:
% -y(:,i-2)=y(:,n_arr-1) -> n_arr>=2
% -y(:,i-w)=y(:,n_arr -1 -w) -> n_arr>=w+1
%Así n_arr=max(w,2) que nos da n_arr +1 puntos de arranque 
n_arr=max(w+1,3);
ts=[t0, t0+ n_arr*s]; 
[~,ys] = rkvectorial(f, ts, y0, d, n_arr);

y(:,1:n_arr +1)=ys; %Ponemos el arranque que hemos calculado antes, tenemos ya 
               %los primeros n_arr+1 valores

maxit=5; %Según la precisión que queramos, augmentamos o no
for i=(n_arr +1):n+1 
    yp=y(:,i-2); %Predictor ya que no disponemos de y(:,w+2), formula implícita
    for j=1:w
        yp=yp+s*b(j+1)*f(t(i-j),y(:,i-j)); %Calculamos la predicción con los 
    end                                    %que sabemos de 0 a w(explícitamente A-B)
    y_nueva=yp;                 
   
    %Ahora hay que resolver la ecuación implícita:
    %y(:,i)=y(:,i-2)+s*(b(1)*f(t(i),y(:,i)) + resto) y usamos la predicción
    %para ello

    for k=1:maxit
        Suma=0;
        for j=1:w
            Suma=Suma + b(j+1)*f(t(i-j),y(:,i-j));
        end
        y_nueva= y(:,i-2) +s*(b(1)*f(t(i),y_nueva) +Suma);
    end
y(:,i)=y_nueva; 
end
end








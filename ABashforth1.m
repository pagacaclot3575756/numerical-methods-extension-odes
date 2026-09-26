function [t,y,q,consistency]=ABashforth1(f,ta,y0,w,n)    %Explícito 
% ABashforth1 - Explicit Adams-Bashforth method
%
% Inputs:
%   f  - function handle, f(t,y), defining the ODE system
%   ta - [t0, T], time interval
%   y0 - initial condition vector
%   w  - number of steps of the Adams-Bashforth method
%   n  - number of subintervals (mesh size)
%
% Output:
%   t           - time mesh
%   y           - approximate solution at each mesh point
%   q           - order of consistency of the method (computed from the
%                 linear-multistep g_l coefficients)
%   consistency - true if the method satisfies the consistency conditions
%
% Starting values are computed with rkvectorial (ideally a Runge-Kutta
% method of order w+1 to avoid degrading the global order). Order and
% consistency are checked by expanding the method in general linear-
% multistep form and testing the associated g_l coefficients.
t0=ta(1);
T=ta(2);

h=(T-t0)/n;

d=length(y0); %Numero de ecuaciones del sistema
y0=y0(:);

b=zeros(w+1,1);

m_ent=(0:w)';
for l=0:w
    m=m_ent; m(l+1)=[]; %Eliminamos la componente l

integrando=@(v) arrayfun(@(x) prod((x+m)./(m-l)), v); %Para que la función se aplique al vector v
b(l+1)=integral(integrando,0,1); %Porque d=1, u=0
end

%Propiedades MML: orden y consistencia 
%Adams-Bashforth con parametro w tiene p=w+1 pasos. 
%y_{i+1}-y_i =h\sum^w_{l=0} b_l^w f_{i-l} 
%Forma general MML:
% \sum^p_{j=0} a_j y_{i+j} = h\sum^p_{j=0}b_jf_{i+j}
p=w+1;
%Cogemos como indice base k=i-w ---> 
% y_{i+1} = y_{k+w+1}=y_{k+p} con coef a_p=1
% y_i = y_{k+w}=y_{k+p-1} con coef a_{p-1}=1
a_mml=zeros(p+1,1);
a_mml(w+1)=-1; %Matlab empieza en a(1)
%a_w=a_{p-1}=-1 --> a(w+1)=-1
a_mml(w+2)=1;

b_mml=zeros(p+1,1);
b_mml(1:w+1)=flipud(b);
b_mml(w+2)=0; %El sumatorio es hasta w, b_{w+1}=0

mmax=2*p; %Orden máximo 2p-1 al ser explicito
%Calculamos los g_l hasta 2p para ver el primero no nulo
g=zeros(mmax,1);
j=(0:p)';

%g_0
g(1)=sum(a_mml);

%g_l l=1,...,q+1
for l=1:mmax
    g(l+1)=sum(a_mml.*(j.^l)-l*b_mml.*(j.^(l-1)));
end

tol=1e-10;
primer_no_nulo=find(abs(g)>tol,1,"first");

if isempty(primer_no_nulo)
    q=mmax;
else
    q=primer_no_nulo-2;
end

if abs(g(1))<=tol && abs(g(2))<=tol
    consistency=true;
else
    consistency=false;
end

%Como habíamos mencionado, hay que calcular los primos w valores con otro
%método monopaso y así inicializar todo. Usaremos RK-4 vectorial, aunque lo 
% ideal sería un método Runge-Kutta de nivel w+1 para no degradar el orden global.

ts=[t0, t0+ w*h]; %Solo hay que calcular los primeros w+1 valores iniciales de y y f
[~,ys] = rkvectorial(f, ts, y0, d, w);

t=linspace(t0,T,n+1);
y=zeros(d,n+1);
y(:,1:w+1)=ys; %Ponemos el arranque que hemos calculado antes, tenemos ya 
               %los primeros w+1 valores
%Ahora ya podemos calcular de y_{w+1} en adelante ya que ya podemos
%calcular f_w,f_{w-1},...,f_0

for i=w+1:n
    Suma=0;
    for l=0:w
        Suma=Suma + b(l+1)*f(t(i-l),y(:,i-l));
    end
    y(:,i+1)=y(:,i)+h*Suma;
end
end


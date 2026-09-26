function [t,y,q,consistency]=AMoulton1(f,ta,y0,w,n,mu)  
% AMoulton1 - Implicit Adams-Moulton method (predictor-corrector, P(EC)^mu)
%
% Inputs:
%   f  - function handle, f(t,y), defining the ODE system
%   ta - [t0, T], time interval
%   y0 - initial condition vector
%   w  - number of steps of the Adams-Moulton method (w=0 -> implicit Euler)
%   n  - number of subintervals (mesh size)
%   mu - number of corrector iterations (P(EC)^mu scheme)
%
% Output:
%   t           - time mesh
%   y           - approximate solution at each mesh point
%   q           - order of consistency of the method (computed from the
%                 linear-multistep g_l coefficients)
%   consistency - true if the method satisfies the consistency conditions
%
% For w=0, reduces to implicit Euler, solved by fixed-point iteration.
% For w>=1, uses an explicit Adams-Bashforth predictor and an implicit
% Adams-Moulton corrector, iterated mu times (P(EC)^mu). Starting values
% are computed with rkvectorial. Order and consistency are checked by
% expanding the method in general linear-multistep form and testing the
% associated g_l coefficients.
t0=ta(1);
T=ta(2);
h=(T-t0)/n;

d=length(y0); %Numero de ecuaciones del sistema
y0=y0(:);

t=linspace(t0,T,n+1);
y=zeros(d,n+1);

if w==0
    b0=1; %Integral de 1 en [-1,0]
    y(:,1)=y0;  
    %Propiedades MML: orden y consistencia 
    %Adams-Moulton con parametro w tiene p=w+1 pasos. 

    %CASO w=0: Euler implícito
    %y_i - y_{i-1}=hf_i 
    p=1;

    a_mml=[-1;1];
    b_mml=[0;1];

    mmax=2*p+1;
    %Calculamos los g_l hasta 2p+1 para ver el primero no nulo
    %(en matlab 2p+2).
    g=zeros(mmax+1,1);
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

    for i=2:n+1 
        yp=y(:,i-1) + h*f(t(i-1),y(:,i-1)); 
%Usamos A-B orden 1 para el predictor, es decir, Euler Explícito
        y_nueva=yp;
        for k=1:mu
            y_nueva= y(:,i-1) + h*(b0*f(t(i),y_nueva));
        end
        y(:,i)=y_nueva; 
    end

else
    %Coeficientes del corrector A-M
    bC=zeros(w+1,1);

    m_ent=(0:w)';
    for l=0:w
        m=m_ent; m(l+1)=[]; %Eliminamos la componente l

        integrando=@(v) arrayfun(@(x) prod((x+m)./(m-l)), v); 
        bC(l+1)=integral(integrando,-1,0); %Porque d=0, u=-1
    end

    %Coeficientes del predictor Adams-Bashforth
    bP=zeros(w+1,1);

    m_ent=(0:w)';
    for l=0:w
        m=m_ent; m(l+1)=[];

        integrando=@(v) arrayfun(@(x) prod((x+m)./(m-l)), v);
        bP(l+1)=integral(integrando,0,1);
    end
    %Propiedades MML: orden y consistencia 
    %Adams-Moulton CASO w>=1
    % y_i-y_{i-1}=h\sum^w_{l=0} b_l^w f_{i-l}
    %
    % Forma general MML:
    % \sum^p_{j=0} a_j y_{k+j} = h\sum^p_{j=0} b_j f_{k+j}
    %
    % Tomamos como indice base k=i-w.
    % Entonces y_i=y_{k+w} y y_{i-1}=y_{k+w-1}.
    % Por tanto, en esta escritura general usamos p=w.
    p=w;
    % y_i = y_{k+w}=y_{k+p} con coef a_{p}=1
    % y_{i-1} = y_{k+w-1}=y_{k+p-1} con coef a_{p-1}=-1
    a_mml=zeros(p+1,1);
    a_mml(w)=-1; %Matlab empieza en a(1)
    a_mml(w+1)=1;
    
    % Coeficientes b_j
    % En el corrector:
    % f_i     tiene coeficiente bC(1)
    % f_{i-1} tiene coeficiente bC(2)
    % ...
    % f_{i-w} tiene coeficiente bC(w+1)
    %
    % Como k=i-w, tenemos:
    % f_{i-w}=f_k,
    % ...
    % f_i=f_{k+w}.
    % Por eso hay que invertir bC(2:end).
    b_mml=zeros(p+1,1);
    b_mml(1:w)=flipud(bC(2:end));
    b_mml(w+1)=bC(1); 

    mmax=2*p+1;
    %Orden máximo 2p al ser implicito
    %Calculamos los g_l hasta 2p+1 para ver el primero no nulo
    %(en matlab 2p+2)
    g=zeros(mmax+1,1);
    j=(0:p)';
    %g_0
    g(1)=sum(a_mml);
    %g_l l=1,...,q
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

    %Inicializamos las iteraciones con Runge-Kutta 4 vectorial.
    ts=[t0, t0+ w*h]; 
    %Solo hay que calcular los primeros w+1 valores iniciales de y y f
    [~,ys] = rkvectorial(f, ts, y0, d, w);

    y(:,1:w+1)=ys; %Ponemos el arranque que hemos calculado antes
    
    %Guardamos tambien los f(y_i,t_i)
    F=zeros(d,n+1);
    for i=1:w+1
        F(:,i)=f(t(i),y(:,i));
    end
    %Bucle predictor-corrector
    for i=w+2:n+1 
        %Predictor Adams-Bashforth:
        %y_i^{[0]}=y_{i-1} + h sum_{l=0}^w bP_l f_{i-1-l}
        yp=y(:,i-1);
        %Predictor ya que no disponemos de y(:,w+2), formula implicita
        for l=0:w
            yp=yp+h*bP(l+1)*F(:,i-1-l); 
        end
        %Calculamos la prediccion con A-B explicito
        y_iter=yp;  

        %Parte conocida del corrector Adams-Moulton
        SumaC=zeros(d,1);

        for l=1:w
            SumaC=SumaC + bC(l+1)*F(:,i-l);
        end
        
        %Ahora hay que resolver la ecuacion implicita:
        %y(:,i)=y(:,i-1)+h*(bC(1)*f(t(i),y(:,i)) + SumaC) 
        %Usamos la prediccion yp como aproximación inicial para ello
        parte_conocida=y(:,i-1) + h*SumaC;

        %Corrector iterado: P(EC)^mu
        for r=1:mu
            y_iter=parte_conocida + h*bC(1)*f(t(i),y_iter);
        end
        
        %Evaluación final E
        y(:,i)=y_iter;
        F(:,i)=f(t(i),y(:,i));
    end
end
end






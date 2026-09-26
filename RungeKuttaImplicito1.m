function [t,y,s,consistency,info] = RungeKuttaImplicito1(f,ta,y0,B,a,n,tol,maxit)
% RungeKuttaImplicito1 - Implicit Runge-Kutta method with Newton-based stage solve
%
% Inputs:
%   f      - function handle, f(t,y), defining the ODE system
%   ta     - [t0, T], time interval
%   y0     - initial condition vector
%   B      - Butcher matrix (b_lj coefficients)
%   a      - Butcher weight vector
%   n      - number of time steps
%   tol    - tolerance for Newton's method (on step or residual norm)
%   maxit  - maximum number of Newton iterations per step
%
% Output:
%   t           - time mesh
%   y           - approximate solution at each mesh point
%   s           - number of stages of the method
%   consistency - true if the method satisfies the order-1 condition
%   info        - struct with stepsize, order info, Newton convergence
%                 history and stage values K
%
% At each step, the stage system K = f(t+ch, y+h*B*K) is solved with
% Newton's method, initialized from an explicit Taylor predictor built
% via generarPhik. The final update is y(i+1) = y(i) + h*sum(a_l*K_l).

% Convencion usada: f(t,Y)
% El metodo viene dado por:
%
%   k_l = f(t_i + c_l h, y_i + h sum_{j=1}^s b_{l,j} k_j),
%   y_{i+1} = y_i + h sum_{l=1}^s a_l k_l.
%
% En cada paso resolvemos el sistema de etapas mediante Newton.

% Datos del intervalo
t0 = ta(1);
T  = ta(2);
h  = (T-t0)/n;

% Coeficientes de Butcher
a = a(:);
s = length(a);

if size(B,1) ~= s || size(B,2) ~= s
    error('La matriz B debe ser de tamano s x s, con s=length(a).');
end

% Guardamos una copia exacta para comprobar el orden
Borden = B;
aorden = a;

% Para el cálculo numérico usamos double
B = double(B);
a = double(a);

% Nodos c_l = sum_j b_lj
c = sum(B,2);

% Dimension del sistema
y0 = y0(:);
d  = length(y0);

% Mallado y matriz de soluciones
t = linspace(t0,T,n+1);
y = zeros(d,n+1);
y(:,1) = y0;

% Comprobacion de orden y consistencia a partir de los coeficientes
infoOrden = comprobarOrdenRKImp(Borden,aorden);
consistency = infoOrden.consistency;

% Para inicializar Newton usamos un predictor de Taylor explicito.
% Si el orden detectado es 0, usamos al menos grado 1.
qPred = max(1, infoOrden.ordenDetectado);

% Evitamos generar demasiadas derivadas: para la prediccion basta con pocas.
qPred = min(qPred,4);

% generarPhik se usa para:
%   phi{k}: predictor de Taylor explicito para los valores de etapa
%   dphi{1}: jacobiana J_y f para Newton
[phi,dphi] = generarPhik(f,d,qPred);
Jf = dphi{1};

% Informacion adicional
info.h = h;
info.B = B;
info.a = a;
info.c = c;
info.s = s;
info.tol = tol;
info.maxit = maxit;
info.orden = infoOrden;
info.qPred = qPred;

info.iterNewton = zeros(1,n);
info.resNewton  = zeros(1,n);
info.convNewton = false(1,n);
info.K = zeros(d,s,n);

% Bucle principal
for i = 1:n

    % Aproximacion inicial de las etapas.
    % Primero aproximamos el valor de la solucion en t_i+c_l h
    % mediante Taylor explicito desde (t_i,y_i), y despues evaluamos f.
    K = zeros(d*s,1);

    for ell = 1:s

        Yell = y(:,i);

        for r = 1:qPred
            phir = phi{r}(t(i),y(:,i));
            Yell = Yell + ((c(ell)*h)^r/factorial(r))*phir(:);
        end

        Kell = f(t(i)+c(ell)*h,Yell);
        K((ell-1)*d+1:ell*d) = Kell(:);

    end

    % Newton para resolver el sistema de etapas F_i(K)=0
    for r = 1:maxit

        F = calcularFRK(t(i),y(:,i),h,K,f,B,c,s,d);
        J = calcularJRK(t(i),y(:,i),h,K,Jf,B,c,s,d);

        % Correccion de Newton: J*delta=-F
        delta = -(J\F);

        K = K + delta;

        % Residuo actualizado
        F = calcularFRK(t(i),y(:,i),h,K,f,B,c,s,d);

        errDelta = norm(delta(:),inf);
        errF = norm(F(:),inf);

        if errDelta < tol || errF < tol
            info.convNewton(i) = true;
            break;
        end

    end

    % Guardamos informacion de Newton
    info.iterNewton(i) = r;
    info.resNewton(i)  = norm(F,inf);

    % Actualizacion final:
    % y_{i+1}=y_i+h sum_l a_l k_l
    incr = zeros(d,1);

    for ell = 1:s
        Kell = K((ell-1)*d+1:ell*d);
        info.K(:,ell,i) = Kell;
        incr = incr + a(ell)*Kell;
    end

    y(:,i+1) = y(:,i) + h*incr;

end

end
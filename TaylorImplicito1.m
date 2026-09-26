function [t,y,q,consistency,info]=TaylorImplicito1(f,ta,y0,grado,n,tol,maxit)
% TaylorImplicito1 - Implicit Taylor method for ODE systems y'=f(t,y)
%
% Inputs:
%   f      - function handle, f(t,Y), defining the ODE system
%   ta     - [t0, T], integration interval
%   y0     - initial condition (column vector, length d)
%   grado  - order q of the implicit Taylor method
%   n      - number of steps (mesh has n+1 points)
%   tol    - Newton tolerance (on step correction and residual, inf-norm)
%   maxit  - maximum Newton iterations allowed per step
%
% Output:
%   t           - 1 x (n+1) time mesh
%   y           - d x (n+1) matrix of computed solution values
%   q           - order used (equal to grado)
%   consistency - logical, consistency check from comprobarOrdenTaylorImp
%   info        - struct with method/order info and per-step Newton
%                 diagnostics (iterNewton, resNewton, convNewton)
%
% Method: y_{i+1} = y_i + h*psi_q(y_{i+1},y_i,t_i,h), solved at each
% step by Newton's method on G_i(Y) = Y - y_i - h*psi_q(Y,y_i,t_i,h) = 0.

% Convencion usada: f(t,Y)
% El metodo se escribe como:
%y_{i+1}=y_i+h*psi_q(y_{i+1},y_i,t_i,h)
% donde:
% psi_q(Z,Y,t,h)=sum_{k=1}^q (-1)^(k+1) h^(k-1)/k! phi_k(Z,t+h)
% Por tanto, en cada paso resolvemos mediante Newton:
% G_i(Y)=Y-y_i-h*psi_q(Y,y_i,t_i,h)=0
q = grado;
t0 = ta(1);
T  = ta(2);
h  = (T-t0)/n;

% Dimensión del sistema
y0 = y0(:);
d  = length(y0);

% Mallado y matriz de soluciones
t = linspace(t0,T,n+1);
y = zeros(d,n+1);
y(:,1) = y0;

% Propiedades: orden y consistencia
% generarPhik construye phi_1,...,phi_{q+1}, pero psi_q y dpsi_q
% solo usan los terminos k=1,...,q. El termino q+1 se necesita para
% comprobar si falla o no la condicion de orden siguiente.
[phi,~,psi,dpsi,phiSym,dphiSym,ySym,tSym] = generarPhik(f,d,q);

infoOrden = comprobarOrdenTaylorImp(phiSym,ySym,tSym,q);
consistency = infoOrden.consistency;

% Información adicional
info.grado = q;
info.h = h;
info.tol = tol;
info.maxit = maxit;
info.orden = infoOrden;
info.iterNewton = zeros(1,n);
info.resNewton = zeros(1,n);
info.convNewton = false(1,n);
info.phiSym = phiSym;
info.dphiSym = dphiSym;
info.ySym = ySym;
info.tSym = tSym;

I = eye(d);

% Bucle principal del método
for i=1:n

    % Predicción inicial para Newton.
    % Usamos Taylor explícito de grado q a partir de y_i:
    %
    % Y^{[0]} = y_i + sum_{k=1}^q h^k/k! phi_k(y_i,t_i)
    %
    % Para q=1 esto es Euler explícito.

    Y = y(:,i);

    for k=1:q
        phik = phi{k}(t(i),y(:,i));
        Y = Y + (h^k/factorial(k))*phik(:);
    end

    Y = Y(:);

    % Newton para resolver G_i(Y)=0
    %maxit=número máximo de iteraciones de Newton 
    % permitidas en cada paso
    for r=1:maxit

        % G_i(Y)=Y-y_i-h*psi_q(Y,y_i,t_i,h)
        % Como psi_q no depende realmente de y_i, usamos psi(t_i,Y,h).
        G = Y - y(:,i) - h*psi(t(i),Y,h);

        % JG_i(Y)=I-h*J_Y psi_q(Y,y_i,t_i,h)
        J = I - h*dpsi(t(i),Y,h);

        % Corrección de Newton: J*delta=-G
        %En sistemas, la derivada escalar \(G_i'\) se sustituye
        % por la matriz jacobiana \(J_G\). Por ello, en lugar de 
        % dividir como en el caso escalar, en cada iteración de 
        % Newton se resuelve un sistema lineal para obtener la 
        % corrección \(\delta^{[r]}\), y posteriormente se 
        % actualiza \(Y^{[r+1]}=Y^{[r]}+\delta^{[r]}\).
        delta = -(J\G);

        Y = Y + delta;
        
        % Residuo actualizado, ya evaluado en el nuevo Y.
        G = Y - y(:,i) - h*psi(t(i),Y,h);

        if norm(delta,inf) < tol || norm(G,inf) < tol
            info.convNewton(i) = true;
            break;
        end

    end

    % Guardamos información de Newton
    info.iterNewton(i) = r;
    info.resNewton(i)  = norm(G,inf);

    % Nueva aproximación
    y(:,i+1) = Y;
end
end
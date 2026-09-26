function J = calcularJRK(ti,yi,h,K,Jf,B,c,s,d)
% calcularJRK - Jacobian of the stage system for an implicit RK method
%
% Inputs:
%   ti,yi - current time and solution value
%   h     - step size
%   K     - current iterate of the stacked stage slopes (d*s x 1)
%   Jf    - function handle for the Jacobian of f w.r.t. y
%   B,c   - Butcher matrix and nodes
%   s,d   - number of stages and system dimension
%
% Output:
%   J - block Jacobian (d*s x d*s) of the residual calcularFRK, with
%       block (ell,m) = delta_{ell,m}*I - h*B(ell,m)*Jf(t+c_ell*h,Y_ell)
%
% Used to build the Newton correction J\F at each iteration.

% Calcula la matriz jacobiana del sistema de etapas RK implicito.
%
% Si el residuo de la etapa ell es:
%
%   F_ell(K) = K_ell - f(t_i+c_ell h, Y_ell),
%
% donde
%
%   Y_ell = y_i + h*sum_j b_{ell,j} K_j,
%
% entonces el bloque (ell,m) de la jacobiana es:
%
%   dF_ell/dK_m = delta_{ell,m} I
%                 - h*b_{ell,m}*J_y f(t_i+c_ell h,Y_ell).

% Aseguramos formato columna.
yi = yi(:);
K  = K(:);

% J es una matriz por bloques de tamano (d*s) x (d*s).
% Cada bloque corresponde a la derivada de F_ell respecto de K_m.
J = zeros(d*s,d*s);

% Identidad de dimension d, usada cuando ell = m.
I = eye(d);

for ell = 1:s

    % Filas del bloque correspondiente a F_ell.
    % Como F = [F_1; F_2; ...; F_s] y cada F_ell tiene dimension d:
    %   F_ell ocupa las filas (ell-1)*d+1 : ell*d.
    filasEll = (ell-1)*d + 1 : ell*d;

    % Construimos el punto interno:
    %   Y_ell = y_i + h*sum_j b_{ell,j} K_j.
    Yell = yi;

    for j = 1:s

        % Posiciones de la etapa K_j dentro del vector grande K.
        filasJ = (j-1)*d + 1 : j*d;
        Kj = K(filasJ);

        Yell = Yell + h*B(ell,j)*Kj;
    end

    % Jacobiana de f respecto de y evaluada en la etapa ell:
    %   J_y f(t_i+c_ell h,Y_ell).
    JfEll = Jf(ti + c(ell)*h, Yell);
    JfEll = double(reshape(JfEll,d,d));

    for m = 1:s

        % Columnas del bloque correspondiente a la variable K_m.
        % Como K = [K_1; K_2; ...; K_s], K_m ocupa estas columnas.
        colsM = (m-1)*d + 1 : m*d;

        % Parte comun del bloque:
        %   -h*b_{ell,m}*J_y f(t_i+c_ell h,Y_ell).
        bloque = -h*B(ell,m)*JfEll;

        % Si ell = m, aparece la derivada de K_ell respecto de K_m,
        % que aporta la identidad:
        %   delta_{ell,m} I.
        if ell == m
            bloque = I + bloque;
        end

        % Colocamos el bloque (ell,m) dentro de la jacobiana total.
        J(filasEll,colsM) = bloque;

    end
end

end
function F = calcularFRK(ti,yi,h,K,f,B,c,s,d)
% calcularFRK - Residual of the stage system for an implicit RK method
%
% Inputs:
%   ti,yi - current time and solution value
%   h     - step size
%   K     - current iterate of the stacked stage slopes (d*s x 1)
%   f     - function handle, f(t,y), defining the ODE system
%   B,c   - Butcher matrix and nodes
%   s,d   - number of stages and system dimension
%
% Output:
%   F - stacked residual F_ell(K) = K_ell - f(t+c_ell*h, Y_ell), with
%       Y_ell = y_i + h*sum_j B(ell,j)*K_j
%
% Used as the nonlinear system Newton's method drives to zero to find
% the internal stage slopes K_1,...,K_s.

% Calcula el residuo del sistema de etapas de un metodo RK implicito.
%
% Para cada etapa ell:
%   K_ell = f(t_i+c_ell h, Y_ell),
%   Y_ell = y_i + h*sum_{j=1}^s b_{ell,j} K_j.
%
% Por tanto, Newton resuelve F(K)=0 con:
%
%   F_ell(K) = K_ell - f(t_i+c_ell h, Y_ell).
%
% La solucion del sistema F(K)=0 proporciona las pendientes internas
% K_1,...,K_s del metodo. K es el vector de estas pendientes, y d la
% dimensión del sistema

% Aseguramos formato columna.
yi = yi(:);
K  = K(:);
%Inicializamos residuo (tamaño ds al ser s etapas de dimensión d cada una
F = zeros(d*s,1);

for ell = 1:s
   % Posiciones del vector grande K que corresponden a la etapa ell.
    % Si cada etapa tiene d componentes:
    %   K_ell = K((ell-1)*d+1 : ell*d).
    filasEll = (ell-1)*d + 1 : ell*d;
    Kell = K(filasEll);

    % Construimos:
    %   Y_ell = y_i + h*sum_j b_{ell,j} K_j.
    Yell = yi;

    for j = 1:s
        % Posiciones correspondientes a la etapa j dentro de K.
        filasJ = (j-1)*d + 1 : j*d;
        Kj = K(filasJ);
        Yell = Yell + h*B(ell,j)*Kj;
    end

    % Residuo de la etapa ell. Evaluacion de f en la etapa ell:
    %   f(t_i+c_ell h, Y_ell).
    fEll = f(ti + c(ell)*h, Yell);
    % Residuo:
    %   F_ell(K) = K_ell - f(t_i+c_ell h, Y_ell).
    F(filasEll) = Kell - fEll(:);
end
end

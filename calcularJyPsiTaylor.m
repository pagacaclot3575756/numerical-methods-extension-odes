function J = calcularJyPsiTaylor(t,Z,h,dphi,grado)
% calcularJyPsiTaylor - Jacobian of the implicit Taylor slope function
%
% Inputs:
%   t     - current time (evaluated at t+h inside)
%   Z     - point at which the Jacobians dphi_k are evaluated
%   h     - step size
%   dphi  - cell array of function handles dphi{k}(t,Y), the Jacobians
%           of the total time-derivatives phi_k (as built by generarPhik)
%   grado - order q of the Taylor method (number of terms to sum)
%
% Output:
%   J - d x d matrix, J_y psi_q(Z,t,h) = sum_{k=1}^q (-1)^(k+1) *
%       h^(k-1)/k! * dphi_k(t+h,Z)
%
% Computes the Jacobian, with respect to Z, of calcularPsiTaylor, for
% use in a Newton iteration that solves the implicit Taylor method.

Z = Z(:);
%Dimensión del sistema
d = length(Z);
%Inicializamos J_y \psi_q como matriz nula
J = zeros(d,d);

for k=1:grado
    coef = (-1)^(k+1)*h^(k-1)/factorial(k);
    Jk = dphi{k}(t+h,Z);
    J = J + coef*reshape(Jk,d,d);
end
end
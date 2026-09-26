function val = calcularPsiTaylor(t,Z,h,phi,grado)
% calcularPsiTaylor - Slope function of the implicit Taylor method
%
% Inputs:
%   t     - current time (evaluated at t+h inside)
%   Z     - point at which the total derivatives phi_k are evaluated
%   h     - step size
%   phi   - cell array of function handles phi{k}(t,Y), the total
%           time-derivatives of the ODE (as built by generarPhik)
%   grado - order q of the Taylor method (number of terms to sum)
%
% Output:
%   val - column vector, psi_q(Z,t,h) = sum_{k=1}^q (-1)^(k+1) *
%         h^(k-1)/k! * phi_k(t+h,Z)
%
% Evaluates the implicit slope function of a degree-q Taylor method at
% (t+h,Z). In this formulation psi_q does not depend on the base point
% Y, only on the implicit unknown Z.

% Calcula psi_q(Z,Y,t,h). En esta formulacion psi_q no depende de Y.
Z = Z(:);
val = zeros(length(Z),1); 
%Inicializamos psi_q como vector columna de ceros

for k=1:grado
    coef = (-1)^(k+1)*h^(k-1)/factorial(k);
    phik = phi{k}(t+h,Z);
    val = val + coef*phik(:);
end
end

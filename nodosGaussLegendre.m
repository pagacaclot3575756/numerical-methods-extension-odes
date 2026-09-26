function c = nodosGaussLegendre(s)
% nodosGaussLegendre - Gauss-Legendre nodes on the interval [0,1]
%
% Input:
%   s - number of nodes (positive integer)
%
% Output:
%   c - s x 1 symbolic vector of nodes in [0,1], sorted ascending
%
% For s=1,2,3 returns exact closed-form symbolic expressions. For s>=4
% uses the Golub-Welsch algorithm (eigenvalues of the Jacobi matrix for
% Legendre polynomials on [-1,1]), then maps to [0,1] and returns
% high-precision symbolic approximations (vpa, 32 digits).

% Calcula los nodos de Gauss-Legendre en el intervalo [0,1].
%
% Para s=1,2,3 se devuelven expresiones simbolicas exactas.
% Para s>=4 se usa el algoritmo de Golub-Welsch y se devuelven nodos
% simbolicos aproximados con alta precision.

if s < 1 || s ~= floor(s)
    error('El numero de etapas s debe ser un entero positivo.');
end

switch s
    case 1
        c = sym(1)/2;

    case 2
        c = [sym(1)/2 - sqrt(sym(3))/6;
             sym(1)/2 + sqrt(sym(3))/6];

    case 3
        c = [sym(1)/2 - sqrt(sym(15))/10;
             sym(1)/2;
             sym(1)/2 + sqrt(sym(15))/10];

    otherwise
        % Para s general, los nodos de Gauss-Legendre se obtienen como
        % los valores propios de la matriz tridiagonal de Jacobi asociada
        % a los polinomios de Legendre, y luego se trasladan de [-1,1] a [0,1].
        k = (1:s-1)';
        beta = k ./ sqrt(4*k.^2 - 1);
        Jac = diag(beta,1) + diag(beta,-1);

        x = sort(eig(Jac));

        % Cambio de variable de [-1,1] a [0,1].
        c = (x + 1)/2;
        c = sym(vpa(c,32));
end
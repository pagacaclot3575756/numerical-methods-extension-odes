function [B,a,c] = construirRKColocacion(c)
% construirRKColocacion - Build collocation Runge-Kutta coefficients
% from a set of nodes
%
% Input:
%   c - length-s vector of collocation nodes c_1,...,c_s in [0,1]
%
% Output:
%   B - s x s coefficient matrix, b_{lj} = int_0^{c_l} L_j(x) dx
%   a - s x 1 weight vector, a_j = int_0^1 L_j(x) dx
%   c - the same node vector, returned as a symbolic column
%
% L_j is the Lagrange basis polynomial for node c_j. This construction
% reduces to the Gauss-Legendre method when c are the Gauss-Legendre
% nodes on [0,1] (see nodosGaussLegendre).

% Construye los coeficientes RK.
% Dados c_1,...,c_s, se considera la base de Lagrange L_j(x) y se define:
% a_j = int_0^1 L_j(x) dx,
% b_lj = int_0^{c_l} L_j(x) dx.
% Esta construccion incluye, como caso particular, los metodos de Gauss
% cuando los nodos c son los nodos gaussianos en [0,1].

%Cogemos el vector c y lo volvemos un columna simbolico
c = sym(c(:));
s = length(c);

%Inicializamos
x = sym('x','real');
B = sym(zeros(s,s));
a = sym(zeros(s,1));

%Cosntruimos el polinomio de Lagrange asociado a cada c_j
for j = 1:s

    Lj = sym(1);

    for m = 1:s
        if m ~= j
            Lj = Lj * (x - c(m))/(c(j)-c(m));
        end
    end
    %Limpiamos la expresión para la integración posterior
    Lj = simplify(expand(Lj));
    
    %a_j = \int_0^1 L_j(x) dx
    a(j) = simplify(int(Lj,x,0,1));

    for ell = 1:s
        %b_{\ell j} = \int_0^{c_\ell} L_j(x) dx
        B(ell,j) = simplify(int(Lj,x,0,c(ell)));
    end
end
end
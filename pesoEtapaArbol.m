function F = pesoEtapaArbol(tau,B)
% pesoEtapaArbol - Elementary stage weight vector associated to a tree
%
% Inputs:
%   tau - rooted tree (as built by crearArbol)
%   B   - Butcher matrix (symbolic or numeric)
%
% Output:
%   F - column vector (size s) with the elementary weight F_i(tau)
%
% For a single node, F(tau) = 1. For tau = [tau_1,...,tau_m],
% F_i(tau) = prod_k (B*F(tau_k))_i, computed recursively over children.

% Calcula el vector de pesos elementales de las etapas asociado a un arbol.
%tau:arbol creado por crear arbol
%B:matriz de coeficientes b_i,j de RK implícito
%F:vector columna simbolico de tamaño s
%
% Si tau es el arbol de un solo nodo:
%       F_i(tau) = 1.
%
% Si tau tiene hijos tau_1,...,tau_m:
%       F_i(tau) = prod_k (B F(tau_k))_i.

%Numero de etapas?
s = size(B,1);

%Caso base: arbol de un solo nodo
if isempty(tau.hijos)
    F = sym(ones(s,1));
%Caso general: arbol con hijos
else
    %Inicializamos
    F = sym(ones(s,1));

    for k = 1:length(tau.hijos)
        %Para cada hijo del arbol se calcula recursivamente
        %su vector de pesos F(tau_k)
        Fhijo = pesoEtapaArbol(tau.hijos{k},B);
        % Si tau=[tau_1,...,tau_m], cada hijo aporta el factor B*F(tau_k)
        %a la etapa i mediante la combinacion lineal de etapas dada por 
        % la matriz B:
        %
        %(B*F(tau_k))_i = sum_j b_ij F_j(tau_k).
        %
        % Como los hijos de la raiz aparecen como factores independientes,
        % sus contribuciones se multiplican componente a componente:
        %
        %F_i(tau) = prod_k (B*F(tau_k))_i.
        % y, al colgar todos de la misma raiz, dichos factores se multiplican
        % componente a componente:
        %
        %F_i(tau)=prod_k (B*F(tau_k))_i
        F = simplify(F .* (B*Fhijo));
        % Este vector F(tau) se usa despues en la condicion de orden
        %
        %a^T F(tau) = 1/gamma(tau).
    end
end
end
function tau = crearArbol(hijos)
% crearArbol - Build a rooted tree node from its subtrees
%
% Inputs:
%   hijos - cell array of child trees (empty for a single-node tree)
%
% Output:
%   tau - struct with fields hijos, orden, gamma (Butcher's tree
%         factorial) and clave (canonical string key, invariant
%         under the order children are given in)
%
% gamma is defined recursively: gamma(single node) = 1,
% gamma([tau_1,...,tau_m]) = |tau| * prod_j gamma(tau_j).

% Crea un arbol enraizado a partir de sus hijos.
%El número \(\gamma(\tau)\), denominado factorial del árbol 
% de Butcher, se define recursivamente por
%\gamma(\bullet)=1,
%\gamma([\tau_1,\dots,\tau_m])=|\tau|\prod_{j=1}^m\gamma(\tau_j).
%
%Este número aparece en el término exacto asociado al árbol 
% \(\tau\) dentro de las condiciones de orden de Runge-Kutta, 
% que toman la forma
%
% a^T\Phi(\tau)=\frac{1}{\gamma(\tau)}.

%Caso base: si hijos vacio, solo una raiz
if isempty(hijos)
    tau.hijos = {};
    tau.orden = 1;
    tau.gamma = sym(1);
    tau.clave = 'o';
    return
end

% Ordenamos los hijos por clave para que el arbol no dependa del orden.
claves = cell(1,length(hijos));
for j = 1:length(hijos)
    claves{j} = hijos{j}.clave;
end
[claves,idx] = sort(claves);
hijos = hijos(idx);

%Inicializamos orden y gamma
orden = 1;
gamma = sym(1);

for j = 1:length(hijos)
    orden = orden + hijos{j}.orden;
    gamma = gamma*hijos{j}.gamma;
end

%Guardamos el arbol completo
tau.hijos = hijos;
tau.orden = orden;
tau.gamma = sym(orden)*gamma;
tau.clave = ['[', strjoin(claves,','), ']'];
end
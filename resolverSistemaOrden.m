function [valores,ok,resMax] = resolverSistemaOrden(exprs,vars,s)
% resolverSistemaOrden - Numerically solve a system of RK order conditions
%
% Inputs:
%   exprs - symbolic column vector of expressions that must equal zero
%           (the order conditions)
%   vars  - symbolic column vector of unknowns (a_i's and B_ij's) to
%           solve for
%   s     - number of stages, used to build the initial guess
%
% Output:
%   valores - symbolic vector of solution values for vars, in the same
%             order as vars (empty sym if no acceptable solution found)
%   ok      - logical, true if a solution with small enough residual
%             (< 1e-20) was found
%   resMax  - maximum absolute residual of the accepted solution
%             (Inf if none was found)
%
% Since the system generally has many solutions, tries up to 25
% deterministically perturbed initial points around a base guess
% (equally spaced nodes, uniform weights) and calls vpasolve at each,
% via extraerValoresSolucion to parse the result, accepting the first
% solution whose residual (evaluated at 50-digit precision) is below
% the tolerance.

% Resuelve numericamente el sistema de condiciones de orden.
%
%expr:vector de expresiones simbolicas que tienen que ser 0
%vars:vecyor de variables sombolicas respecto a las que resolver
%s:numero de nieveles
%
%ok:true si se ha encontrado una solucion aceptable
%valores:valores de las variables vars encontados
%resMax:residuo maximo de la solucion aceptada despues de 
%    sustituirla en la ecuacion.
%
% Como el sistema puede tener muchas soluciones, se prueban varios puntos
% iniciales. Se acepta una solucion si el residuo maximo es pequeno.

%Inicializamos
ok = false;
resMax = inf;
valores = sym([]);

%max|exprs|<tol
tol = 1e-20;
%Se probaran 25 puntos iniciales distintos
numIntentos = 25;

% Punto inicial base.
%Tomamos nodos repartidos en [0,1]
cBase = (1:s)'/(s+1);
aBase = ones(s,1)/s;
%construye una matriz inicial B usando los valores de cBase
BBase = repmat(cBase/s,1,s);

%Se junta todo en un unico vector
xBase = [aBase; BBase(:)];

for intento = 1:numIntentos

    % Perturbacion determinista para no depender de numeros aleatorios.
    perturbacion = 0.1*sin((1:length(xBase))'*(intento+1));
    x0 = xBase + perturbacion;

    try
        sol = vpasolve(exprs == 0, vars, x0);
    catch
        continue
    end

    if isempty(sol)
        continue
    end
    %Lo volvemso un vector de valores ordenado
    vals = extraerValoresSolucion(sol,vars);

    if isempty(vals)
        continue
    end

    % Comprobamos residuo poniendo los valores encontrados en 
    %las ecuaciones exprs con 50 dígitos de precisión
    res = vpa(subs(exprs,vars,vals),50);

    try
        resActual = double(max(abs(res)));
    catch
        continue
    end

    if resActual < tol
        valores = vals;
        resMax = resActual;
        ok = true;
        return
    end
end

end
function [B,a,c,infoOrden] = calcularCoefRKOrden(s,p)
% calcularCoefRKOrden - Numerically solve the RK order conditions for a
% given number of stages and target order
%
% Inputs:
%   s - number of stages
%   p - desired order (must satisfy p < 2*s; use mode 'gauss' in
%       coefRKImplicito to reach the maximal order 2*s)
%
% Output:
%   B         - s x s coefficient matrix (symbolic)
%   a         - s x 1 weight vector (symbolic)
%   c         - s x 1 node vector, c = B*ones(s,1) (symbolic)
%   infoOrden - struct with fields modo, ordenSolicitado,
%               ordenGarantizado, etapas, familia, numeroCondiciones,
%               numeroVariables, residuoMaximo, descripcion, observacion
%
% Builds the order conditions a^T*F(tau) = 1/gamma(tau) for every rooted
% tree tau up to order p (via generarArbolesHastaOrden and
% pesoEtapaArbol), then finds one particular numeric solution with
% resolverSistemaOrden. The solution is not unique, and for high p the
% underlying nonlinear system may be hard to solve; errors out if no
% solution with small enough residual is found.

%Construye coeficientes RK implicitos para alcanzar
%al menos un orden p, sin introducir B ni a manualmente.
%
% En este modo se plantean automaticamente las condiciones algebraicas
% de orden de Runge-Kutta hasta orden p, y se busca una solucion particular.
%
% Importante:
%   - La eleccion de coeficientes no es unica.
%   - Para alcanzar el orden maximo 2s debe usarse el modo 'gauss'.
%   - Para ordenes altos, el sistema no lineal puede ser dificil de resolver.

if p < 1 || p ~= floor(p)
    error('El orden p debe ser un entero positivo.');
end

if s < 1 || s ~= floor(s)
    error('El numero de etapas s debe ser un entero positivo.');
end

if p > 2*s
    error('Con s etapas no se puede alcanzar orden mayor que 2s.');
end

if p == 2*s
    error('Para alcanzar el orden maximo 2s usa directamente el modo ''gauss''.');
end

% Variables simbolicas:
%   a_i     pesos del metodo
%   B_ij    coeficientes internos
aSym = sym('a',[s 1],'real');
BSym = sym('b',[s s],'real');

% En un metodo RK se toma c = B*1.
cSym = simplify(BSym*sym(ones(s,1)));

% Generamos todos los arboles enraizados hasta orden p.
arboles = generarArbolesHastaOrden(p);

% Construimos las condiciones de orden.
exprs = sym.empty(0,1);

for q = 1:p
    lista = arboles{q};

    for r = 1:length(lista)
        tau = lista{r};

        % Peso elemental del metodo RK asociado al arbol tau.
        Ftau = pesoEtapaArbol(tau,BSym);

        % Condicion de orden:
        %       a^T F(tau) = 1/gamma(tau)
        exprs(end+1,1) = simplify(aSym.'*Ftau - sym(1)/tau.gamma);
    end
end

vars = [aSym(:); BSym(:)];

[valores,ok,resMax] = resolverSistemaOrden(exprs,vars,s);

if ~ok
    error(['No se ha encontrado una solucion numerica para las condiciones ', ...
           'de orden solicitadas. Prueba con otro numero de etapas, ', ...
           'con el modo ''gauss'', o introduce los coeficientes manualmente.']);
end

% Sustituimos la solucion encontrada.
a = simplify(subs(aSym,vars,valores));
B = simplify(subs(BSym,vars,valores));
c = simplify(subs(cSym,vars,valores));

infoOrden.modo = 'orden';
infoOrden.ordenSolicitado = p;
infoOrden.ordenGarantizado = p;
infoOrden.etapas = s;
infoOrden.familia = 'condiciones de orden por arboles enraizados';
infoOrden.numeroCondiciones = length(exprs);
infoOrden.numeroVariables = length(vars);
infoOrden.residuoMaximo = resMax;
infoOrden.descripcion = ['Coeficientes obtenidos resolviendo automaticamente ', ...
    'las condiciones algebraicas de orden de Runge-Kutta hasta orden ', ...
    num2str(p), '.'];
infoOrden.observacion = ['La solucion no es unica; el programa devuelve una ', ...
    'solucion particular encontrada numericamente.'];
end
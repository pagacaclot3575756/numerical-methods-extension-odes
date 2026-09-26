function info = comprobarOrdenRKImp(B,a,maxOrden,tol)
% comprobarOrdenRKImp - Order and consistency check for an implicit RK method
%
% Inputs:
%   B        - Butcher matrix (b_ij coefficients)
%   a        - Butcher weight vector
%   maxOrden - maximum order to check (optional, default 2s+1)
%   tol      - numerical tolerance for order conditions (optional)
%
% Output:
%   info - struct with consistency flag, detected order, which order
%          conditions hold/fail, their residuals and a text summary
%
% Verifies the rooted-tree order conditions a^T*F(tau) = 1/gamma(tau)
% for every tree up to maxOrden, using generarArbolesHastaOrden to
% enumerate the trees and pesoEtapaArbol to evaluate each F(tau).

% Comprueba consistencia y orden de un metodo Runge-Kutta implicito.
%
% Se usan las condiciones de orden mediante arboles enraizados:
%
%   a^T F(tau) = 1/gamma(tau),
%
% donde F(tau) es el peso elemental de etapas asociado al arbol tau.
%
% Entradas:
%   B        matriz de coeficientes b_{ij}
%   a        vector de pesos a_i
%   maxOrden orden maximo que se quiere comprobar, opcional
%   tol      tolerancia numerica, opcional
%
% Salida:
%   info estructura con los campos principales:
%       consistency
%       ordenDetectado
%       cumpleOrden
%       residuosOrden
%       descripcionOrden

if nargin < 4 || isempty(tol)
    tol = 1e-10;
end

% Vector columna de pesos.
a = a(:);
s = length(a);

if ~isequal(size(B),[s,s])
    error('La matriz B debe tener dimension s x s, con s=length(a).');
end

if nargin < 3 || isempty(maxOrden)
    % Para un metodo RK de s etapas el orden no puede superar 2s.
    % Comprobamos una condicion mas para poder distinguir orden exacto
    % frente a orden al menos 2s en los casos habituales.
    maxOrden = 2*s + 1;
end

if maxOrden < 1 || maxOrden ~= floor(maxOrden)
    error('maxOrden debe ser un entero positivo.');
end

% Trabajamos en simbolico para aprovechar simplificacion exacta cuando sea
% posible, pero aceptamos tambien coeficientes numericos aproximados.
Bsym = sym(B);
asym = sym(a);

arboles = generarArbolesHastaOrden(maxOrden);

cumpleOrden = false(1,maxOrden);
residuosOrden = inf(1,maxOrden);
numeroCondicionesOrden = zeros(1,maxOrden);
condiciones = cell(maxOrden,1);

for q = 1:maxOrden

    lista = arboles{q};
    numeroCondicionesOrden(q) = length(lista);

    cumpleEsteOrden = true;
    residuoMaxQ = 0;
    condiciones{q} = struct('clave',{},'gamma',{},'residuo',{},'cumple',{});

    for r = 1:length(lista)

        tau = lista{r};
        Ftau = pesoEtapaArbol(tau,Bsym);

        % Condicion de orden asociada a tau.
        expr = simplify(asym.'*Ftau - sym(1)/tau.gamma);

        [cumpleCondicion,residuo] = esCeroOrdenRK(expr,tol);

        condiciones{q}(r).clave = tau.clave;
        condiciones{q}(r).gamma = tau.gamma;
        condiciones{q}(r).residuo = residuo;
        condiciones{q}(r).cumple = cumpleCondicion;

        residuoMaxQ = max(residuoMaxQ,residuo);

        if ~cumpleCondicion
            cumpleEsteOrden = false;
        end
    end

    cumpleOrden(q) = cumpleEsteOrden;
    residuosOrden(q) = residuoMaxQ;

end

% El orden detectado es el mayor q para el que se cumplen todas las
% condiciones hasta q.
ordenDetectado = 0;

for q = 1:maxOrden
    if cumpleOrden(q)
        ordenDetectado = q;
    else
        break
    end
end

consistency = ordenDetectado >= 1;

if ordenDetectado == 0
    tipoOrden = 'no consistente';
    descripcionOrden = 'El metodo no es consistente: falla alguna condicion de orden 1.';
elseif ordenDetectado == maxOrden
    tipoOrden = 'al menos';
    descripcionOrden = ['Orden al menos ',num2str(ordenDetectado), ...
        ' dentro de las condiciones comprobadas.'];
else
    tipoOrden = 'exactamente';
    descripcionOrden = ['Orden exactamente ',num2str(ordenDetectado), ...
        ', porque falla alguna condicion de orden ',num2str(ordenDetectado+1),'.'];
end

info = struct();
info.consistency = consistency;
info.ordenDetectado = ordenDetectado;
info.tipoOrden = tipoOrden;
info.descripcionOrden = descripcionOrden;
info.maxOrdenComprobado = maxOrden;
info.tol = tol;
info.cumpleOrden = cumpleOrden;
info.residuosOrden = residuosOrden;
info.numeroCondicionesOrden = numeroCondicionesOrden;
info.condiciones = condiciones;
info.B = B;
info.a = a;
info.c = simplify(Bsym*sym(ones(s,1)));

% Campos utiles para mensajes posteriores.
info.ordenAlMenos1 = consistency;
info.ordenAlMenosDetectado = ordenDetectado;
end



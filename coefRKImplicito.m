 function [B,a,c,infoCoef] = coefRKImplicito(s,modo,p,varargin)
% coefRKImplicito - Build or register the coefficients of an implicit
% Runge-Kutta method
%
% Inputs:
%   s        - number of stages
%   modo     - how the coefficients are chosen:
%              'manual' - coefficients supplied directly via
%                         'B',B,'a',a name-value pairs (see varargin)
%              'orden'  - coefficients computed by imposing the
%                         algebraic order conditions (rooted trees) up
%                         to order p
%              'gauss'  - Gauss-Legendre collocation method with s
%                         stages (theoretical order 2s)
%   p        - desired order, only used (and required) in 'orden' mode;
%              omit or pass [] otherwise
%   varargin - optional name-value pairs, only used in 'manual' mode:
%              'B',B - s x s coefficient matrix b_{lj}
%              'a',a - length-s weight vector a_l
%              (nodes are defined as c_l = sum_j b_{lj})
%
% Output:
%   B        - s x s coefficient matrix
%   a        - s x 1 weight vector
%   c        - s x 1 node vector, c = sum(B,2)
%   infoCoef - struct with fields modo, descripcion, familia,
%              ordenSolicitado, ordenTeorico (NaN for manual mode,
%              since order/consistency are not checked here; use
%              comprobarOrdenRKImp separately)

% Construye o registra los coeficientes de un metodo Runge-Kutta implicito.
% MODOS IMPLEMENTADOS
% 1) manual:
%       [B,a,c,infoCoef] = coefRKImplicito(s,'manual','B',B,'a',a)
%    El usuario introduce los coeficientes y el programa solo comprueba
%    consistencia y orden.
% 2) orden:
%       [B,a,c,infoCoef] = coefRKImplicito(s,'orden',p)
%    Se plantean automaticamente las condiciones algebraicas de orden
%    de Runge-Kutta mediante arboles enraizados hasta orden p,
%    y se busca una solucion particular.
% 3) gauss:
%       [B,a,c,infoCoef] = coefRKImplicito(s,'gauss')
%    Se construye el metodo de Gauss-Legendre de s etapas, de orden teorico 2s.
%
%Argumentos de entrada:
% s: número de etapas o niveles del método.
% p: orden deseado. Solo se usa en modo 'orden'.
% modo: forma de elegir los coeficientes. Puede ser:
%      'manual'-> se introducen B y a directamente.
%      'orden'-> se calculan coeficientes imponiendo condiciones de orden.
%      'gauss'-> se usan los coeficientes de Gauss-Legendre,
%                obteniendo orden teórico 2s.
% varargin: argumentos opcionales. Se pasan por parejas nombre-valor.
% Por ejemplo:
% 'B',B matriz de coeficientes b_{lj}, usada en modo manual.
% 'a',a vector de pesos a_l, usado en modo manual.
%Usamos la convención siguiente
%   c_l = sum_j b_lj.

modo = normalizarModo(modo);

% Permite llamar:
%   coefRKImplicito(s,'manual','B',B,'a',a)
% sin tener que poner p=[].
if nargin < 3
    p = [];
    argumentosOpcionales = {};
else
    if strcmp(modo,'manual') && (ischar(p) || isstring(p))
        argumentosOpcionales = [{p}, varargin];
        p = [];
    else
        argumentosOpcionales = varargin;
    end
end

if s < 1 || s ~= floor(s)
    error('El numero de etapas s debe ser un entero positivo.');
end

switch modo

    case 'manual'

        opts = leerOpciones(argumentosOpcionales{:});

        if ~isempty(p)
            error('En modo manual no debes pasar orden p. Usa solo ''B'',B,''a'',a.');
        end

        if ~isfield(opts,'B') || ~isfield(opts,'a')
            error('En modo manual debes pasar ''B'',B,''a'',a.');
        end

        B = opts.B;
        a = opts.a(:);
        %Comprobamos que los tamaños de B y a sean compatibles con 
        %el numero de etapas s
        comprobarDimensiones(B,a,s);
        %Calculamos el vector de nodos sumando cada fila de B
        c = sum(B,2);

        descripcion = 'Coeficientes introducidos manualmente.';
        %No sabemos automáticamente el orden, 
        %en este caso hay que usar luego
        %comprobarOrdenRKImp
        ordenTeorico = NaN;
        ordenSolicitado = [];
        familia = 'manual';

    case 'orden'

        if isempty(p)
            error('En modo orden debes pasar el orden deseado p como tercer argumento.');
        end

        if ~isempty(argumentosOpcionales)
            error('En modo orden no deben pasarse coeficientes manuales ni opciones extra.');
        end

        [B,a,c,infoOrden] = calcularCoefRKOrden(s,p);

        descripcion = infoOrden.descripcion;
        ordenTeorico = infoOrden.ordenGarantizado;
        ordenSolicitado = p;
        familia = infoOrden.familia;

    case 'gauss'

        if ~isempty(p)
            error('En modo gauss no debes pasar p. El orden teorico es automaticamente 2s.');
        end

        if ~isempty(argumentosOpcionales)
            error('En modo gauss no deben pasarse coeficientes manuales ni opciones extra.');
        end

        c = nodosGaussLegendre(s);
        [B,a,c] = construirRKColocacion(c);

        descripcion = 'Metodo de colocacion de Gauss-Legendre.';
        ordenTeorico = 2*s;
        ordenSolicitado = 2*s;
        familia = 'Gauss-Legendre';


    otherwise
        error('Modo no reconocido. Usa ''manual'', ''orden'' o ''gauss''.');
end
%Guardamos todo para la salida
infoCoef = struct();
infoCoef.modo = modo;
infoCoef.descripcion = descripcion;
infoCoef.familia = familia;
infoCoef.ordenSolicitado = ordenSolicitado;
infoCoef.ordenTeorico = ordenTeorico;
end
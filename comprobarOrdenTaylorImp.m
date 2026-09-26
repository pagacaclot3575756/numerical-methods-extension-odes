function infoOrden = comprobarOrdenTaylorImp(phiSym,ySym,tSym,grado)
% comprobarOrdenTaylorImp - Symbolic order/consistency check for the
% implicit Taylor method
%
% Inputs:
%   phiSym - cell array of symbolic total derivatives phi_1,...,phi_{q+1}
%            (as produced by generarPhik)
%   ySym   - symbolic column vector of state variables
%   tSym   - symbolic time variable
%   grado  - order q of the method being checked
%
% Output:
%   infoOrden - struct with fields:
%     consistency          - logical, varphi_q(0) == phi_1 (i.e. f)
%     condOrden             - 1 x q logical, order condition satisfied at
%                              each step s=1,...,q
%     diferenciasOrden      - cell array, symbolic difference for each s
%     ordenAlMenosQ         - logical, all condOrden(1:q) true
%     diferenciaOrdenSiguiente - symbolic difference for the order-(q+1)
%                              condition
%     fallaOrdenQmas1       - logical, true if order q+1 is NOT reached
%     grado                 - q
%     ordenDetectado        - detected order for this problem
%     tipoOrden             - 'exactamente' | 'al menos' | 'ninguno'
%     descripcionOrden      - human-readable description string
%
% Checks symbolically, via Taylor expansion in h of
% varphi_q(h) = psi_q(y(t+h),y(t),t,h), whether the order conditions of
% a degree-q implicit Taylor method hold for the given ODE.

% Comprueba simbolicamente las condiciones de orden y consistencia
% del metodo de Taylor implicito.
%
%Las condiciones de orden son:
% varphi_q(h) = psi_q(y(t+h),y(t),t,h)
%
% 1/s! * y^(s)(t) = 1/(s-1)! * phi_q^(s-1)(0),
% s=1,...,q.
%
% Como y^(s)(t)=phi_s(y(t),t), esto equivale a que el coeficiente
% de h^(s-1) en varphi_q(h) sea 1/s! * phi_s(y,t).
q=grado;
hSym = sym('h','real');

% Construimos el desarrollo de y(t+h) hasta orden q+1.
% y(t+h) = y + h phi_1(y,t) + h^2/2! phi_2(y,t) + ...
% para comprobar simbolicamente las condiciones.
Yh = ySym;

%Vemos hasta q+1 para ver si falla o no la condición asociada.
for r=1:q+1
    Yh = Yh + hSym^r/factorial(r)*phiSym{r};
end

% Construimos varphi_q(h)=psi_q(y(t+h),y(t),t,h):
%
% psi_q = sum_{k=1}^q (-1)^(k+1) h^(k-1)/k! phi_k(y(t+h),t+h)

varphiq = sym(zeros(length(ySym),1));

varsAntiguas = [tSym; ySym(:)];
varsNuevas   = [tSym + hSym; Yh(:)];

for k=1:q

    phi_k_tmash = subs(phiSym{k},varsAntiguas,varsNuevas);
    %Hacemos taylor en phi_k(y(t+h),t+h) respecto h alrededor de 0
    phi_k_tmash = taylor(phi_k_tmash,hSym,'ExpansionPoint',0,'Order',q+1);

    factor_k = (-1)^(k+1)*hSym^(k-1)/factorial(k);
    varphiq = varphiq + factor_k*phi_k_tmash;

end
%Reordenamos y desarrollamos la suma completa en potencias de h.
varphiq = simplify(taylor(varphiq,hSym,'ExpansionPoint',0,'Order',q+1));

% Comprobacion de consistencia:
% varphi_q(0)=phi_1(y,t)=f(y,t)

varphi0 = simplify(subs(varphiq,hSym,0));%Evaluas en h=0 y simplificas

infoOrden.consistency = esCeroSimbolico(varphi0 - phiSym{1});

% Comprobacion de las condiciones de orden s=1,...,q

infoOrden.condOrden = false(q,1);
infoOrden.diferenciasOrden = cell(q,1);

for s=1:q

    % Coeficiente de h^(s-1) en phi_q(h)
    coef_s = simplify(subs(diff(varphiq,hSym,s-1),hSym,0)/factorial(s-1));

    % Coeficiente teorico esperado: 1/s! phi_s(y,t)
    esperado_s = simplify(phiSym{s}/factorial(s));

    diferencia = simplify(coef_s - esperado_s);

    infoOrden.diferenciasOrden{s} = diferencia;
    infoOrden.condOrden(s) = esCeroSimbolico(diferencia);

end
%Comprobamos si todas las condiciones de 1 a q son verdad
infoOrden.ordenAlMenosQ = all(infoOrden.condOrden);

% Comprobacion del termino siguiente.
% Para que el metodo tuviera orden al menos q+1, tendria que cumplirse que
% el coeficiente de h^q en phi_q(h) fuera 1/(q+1)! phi_{q+1}(y,t).

coef_q = simplify(subs(diff(varphiq,hSym,q),hSym,0)/factorial(q));
esperado_qmas1 = simplify(phiSym{q+1}/factorial(q+1));

infoOrden.diferenciaOrdenSiguiente = simplify(coef_q - esperado_qmas1);

% Si esta diferencia no es cero, el metodo no tiene orden al menos q+1 para el
% problema concreto. En general, por la teoria de main, el metodo de Taylor
% implicito de grado q es de orden exactamente q.

infoOrden.fallaOrdenQmas1 = ~esCeroSimbolico(infoOrden.diferenciaOrdenSiguiente);

infoOrden.grado = q;

if infoOrden.ordenAlMenosQ

    if infoOrden.fallaOrdenQmas1
        % Se cumplen las condiciones hasta q,
        % pero falla la condición de orden q+1.
        infoOrden.ordenDetectado = q;
        infoOrden.tipoOrden = 'exactamente';
        infoOrden.descripcionOrden = ['Orden exactamente ', num2str(q), ...
            ' para este problema.'];

    else
        % Se cumplen las condiciones hasta q
        % y también se cumple la condición siguiente.
        infoOrden.ordenDetectado = q + 1;
        infoOrden.tipoOrden = 'al menos';
        infoOrden.descripcionOrden = ['Orden al menos ', num2str(q+1), ...
            ' para este problema.'];
    end

else

    % Si falla antes de llegar a q, buscamos hasta dónde se cumplen
    % las condiciones consecutivamente.
    qaux = 0;

    for s = 1:q
        if infoOrden.condOrden(s)
            qaux = s;
        else
            break;
        end
    end

    infoOrden.ordenDetectado = qaux;

    if qaux == 0
        infoOrden.tipoOrden = 'ninguno';
        infoOrden.descripcionOrden = ...
            'No se cumple ni la primera condición de orden.';
    else
        % Si falla la condición siguiente, el orden detectado es exacto.
        infoOrden.tipoOrden = 'exactamente';
        infoOrden.descripcionOrden = ['Orden exactamente ', num2str(qaux), ...
            ' para este problema.'];
    end

end
end
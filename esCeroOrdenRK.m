function [tf,residuo] = esCeroOrdenRK(expr,tol)
% esCeroOrdenRK - Test whether a symbolic/numeric expression is zero
%
% Inputs:
%   expr - symbolic or numeric expression to test
%   tol  - tolerance used for the numeric fallback check
%
% Output:
%   tf      - true if expr is (numerically or symbolically) zero
%   residuo - residual used to decide, 0 for an exact symbolic zero
%
% Tries an exact symbolic zero test first (isAlways); if that fails,
% falls back to a high-precision numeric evaluation (vpa).

% Decide si una expresion simbolica/numerica es cero con tolerancia.

expr = simplify(expr);

% Primero intentamos una comprobacion simbolica exacta.
try
    tfExacto = all(isAlways(expr(:)==0,'Unknown','false'));
    if tfExacto
        tf = true;
        residuo = 0;
        return
    end
catch
    % Si la version de MATLAB no acepta la opcion Unknown, probamos abajo.
end

try
    tfExacto = all(isAlways(expr(:)==0));
    if tfExacto
        tf = true;
        residuo = 0;
        return
    end
catch
end

% Si no se puede demostrar simbolicamente, evaluamos numericamente.
try
    val = double(vpa(expr,50));
    residuo = max(abs(val(:)));

    if isempty(residuo)
        residuo = 0;
    end

    tf = isfinite(residuo) && residuo < tol;
catch
    residuo = inf;
    tf = false;
end

end
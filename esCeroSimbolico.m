function tf = esCeroSimbolico(expr)
% esCeroSimbolico - Check whether a symbolic expression is identically zero
%
% Input:
%   expr - symbolic scalar, vector, or matrix expression
%
% Output:
%   tf - logical scalar, true if every entry of expr simplifies to zero
%        (checked with isAlways, falling back to a version without the
%        'Unknown','false' option if that syntax is unsupported)

% Comprueba si una expresion simbolica vectorial/matricial es cero.

expr = simplify(expr);

try
    tf = all(isAlways(expr(:)==0,'Unknown','false'));
catch
    tf = all(isAlways(expr(:)==0));
end
tf = logical(tf);
end
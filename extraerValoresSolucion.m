function vals = extraerValoresSolucion(sol,vars)
% extraerValoresSolucion - Extract a value vector from a vpasolve result
%
% Inputs:
%   sol  - output of vpasolve: either a symbolic array or a struct with
%          one field per solved variable
%   vars - symbolic vector of variables the solution should provide
%          values for, in the desired output order
%
% Output:
%   vals - length(vars) x 1 symbolic vector of values, in the order of
%          vars. Returns an empty sym (sym.empty(0,1)) if sol is not a
%          struct/sym, has the wrong size, or is missing any variable
%          named in vars.

% Extrae los valores de una solucion sol devuelta por vpasolve.
%Inicializamos
vals = sym.empty(0,1);
%Caso en que vpasolve devuelve vector simbolico
if isa(sol,'sym')
    if numel(sol) == numel(vars)
        vals = sol(:);
    end
    return
end

%Si no es ni simbolica ni estructura, se descarta
if ~isstruct(sol)
    return
end

vals = sym(zeros(length(vars),1));

for k = 1:length(vars)

    nombre = char(vars(k));
    %Comprobamos si el nombre existe, y en ese caso 
    % se guarda su valor
    if isfield(sol,nombre)
        vals(k) = sol.(nombre);
    %Si falta alguna variable, la solucion se rechaza
    else
        vals = sym.empty(0,1);
        return
    end
end
end
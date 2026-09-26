function combinaciones = combinarHijos(lista,restante,idxMin)
% combinarHijos - Enumerate multisets of trees with a given total order
%
% Inputs:
%   lista    - available trees to use as children
%   restante - remaining order that the multiset must sum to
%   idxMin   - minimum index in lista allowed (enforces non-decreasing
%              indices, avoiding repeated multisets in different order)
%
% Output:
%   combinaciones - cell array, each entry a cell array of child trees
%                   whose orders sum to restante (base case restante==0
%                   returns {{}}, the empty multiset)


% Genera combinaciones con repeticion de arboles cuyo orden total
% sea igual a restante.
%
%lista:arboles disponibles para usar como hijos.
%restante:orden total que falta por completar.
%idxMin: indica desde que posicion de la lista se puede seguir.
%
% idxMin impone que los indices sean no decrecientes, para evitar
% repetir el mismo multiconjunto en distinto orden.

%Inicializamos la celda de combinaciones, donde cada
%elemento será una celda con varios arboles
combinaciones = {};

if restante == 0
    combinaciones = {{}};
    return
end

for idx = idxMin:length(lista)
    %Vemos el orden del arbol que estamos considerando añadir.
    ordenActual = lista{idx}.orden;
    %Solo se añade si su orden no supera el que falta.
    if ordenActual <= restante
        %El orden restante nuevo es restante-ordenActual
        %Seguimos lalmando recursivamente hasta que no entren más
        subcombinaciones = combinarHijos(lista,restante-ordenActual,idx);

        for j = 1:length(subcombinaciones)
            %Se añade el arbol actual delante de cada subcombinación
            %De esta forma se hacen todas las combinaciones completas.
            combinaciones{end+1} = [{lista{idx}}, subcombinaciones{j}];
        end
    end
end
end
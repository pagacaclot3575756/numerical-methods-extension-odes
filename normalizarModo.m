function modo = normalizarModo(modo)
% normalizarModo - Normalize a mode string to a canonical form
%
% Input:
%   modo - mode name as char or string (e.g. 'Gauss', 'gauss-legendre')
%
% Output:
%   modo - lowercased char, with leading/trailing whitespace trimmed
%          and internal spaces, underscores and hyphens removed, so
%          that simple spelling/formatting variants are all accepted
%          by coefRKImplicito's switch statement.

% Limpia el texto del modo para aceptar variantes sencillas.

modo = lower(strtrim(char(modo)));
modo = strrep(modo,' ','');
modo = strrep(modo,'_','');
modo = strrep(modo,'-','');

end
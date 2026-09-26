function opts = leerOpciones(varargin)
% leerOpciones - Parse name-value pair arguments into a struct
%
% Input:
%   varargin - arguments given as name-value pairs (even count required)
%
% Output:
%   opts - struct with one field per option. Names 'b'/'a' (any case)
%          are normalized to fields 'B'/'a'; any other name is passed
%          through matlab.lang.makeValidName to become a valid field
%          name.
%
% Errors if the number of arguments is odd or a name is not text.

% Lee opciones dadas por parejas nombre-valor.
%Creamos estructura vacía para guardar las opciones.
opts = struct();

if mod(length(varargin),2) ~= 0
    error('Las opciones deben pasarse por parejas nombre-valor.');
end

for k = 1:2:length(varargin)

    nombre = varargin{k};

    if ~(ischar(nombre) || isstring(nombre))
        error('El nombre de cada opcion debe ser texto.');
    end

    %Convertimos a caracter la string nombre
    nombre = char(nombre);

    switch lower(nombre) %Lo ponemos en minusculas
        case 'b'
            nombreCampo = 'B';
        case 'a'
            nombreCampo = 'a';
        otherwise
            nombreCampo = matlab.lang.makeValidName(nombre);
    end

    opts.(nombreCampo) = varargin{k+1};
end
end
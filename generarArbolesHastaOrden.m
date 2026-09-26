function arboles = generarArbolesHastaOrden(p)
% generarArbolesHastaOrden - Enumerate rooted trees up to a given order
%
% Inputs:
%   p - maximum tree order to generate
%
% Output:
%   arboles - cell array, arboles{q} holds all rooted trees of order q
%
% Builds each new order by attaching a root to every valid multiset of
% smaller trees (via combinarHijos) and wrapping it with crearArbol,
% discarding duplicates by their canonical key.

% Genera todos los arboles enraizados no ordenados hasta orden p.
%
% Un arbol se representa mediante una estructura con campos:
%   hijos  : celda con sus subarboles
%   orden  : numero de nodos
%   gamma  : factorial/arbol de Butcher
%   clave  : cadena que identifica el arbol de forma unica

%Creamos una celda de p posiciones con p el orden deseado
% para guardar todos los arboeles que creemos.
arboles = cell(p,1);

% Arbol de un solo nodo.
arboles{1} = {crearArbol({})};

for n = 2:p
    % Lista de todos los arboles de orden menor que n.
    %La idea es crear un arbol de orden n creando una raiz nueva
    % y colgando subarboles suya suma de ordenes de n-1.
    lista = {};
    for q = 1:n-1
        %Creamos lista con todos los arboles conocidos de orden
        %menor que n
        lista = [lista, arboles{q}];
    end

    % Ordenamos por clave para generar multiconjuntos de forma canonica.
    claves = cell(1,length(lista));
    for j = 1:length(lista)
        claves{j} = lista{j}.clave;
    end
    [~,idx] = sort(claves);
    lista = lista(idx);

    % Un arbol de orden n se obtiene poniendo una raiz sobre un
    % multiconjunto de subarboles con orden total n-1.
    combinaciones = combinarHijos(lista,n-1,1);
    %Evitamos duplicados con el diccionario mapa de claves de los arboles 
    %que ya han aparecido
    mapa = containers.Map('KeyType','char','ValueType','logical');
    nuevos = {};
    %Para cada combinacion de hijos anterior, creamos un arbol nuevo
    for j = 1:length(combinaciones)
        tau = crearArbol(combinaciones{j});
           %Ese arbol ya existe? Lo mira en mapa, y 
           % si no existia se guarda
        if ~isKey(mapa,tau.clave)
            mapa(tau.clave) = true;
            nuevos{end+1} = tau; 
        end
    end
    arboles{n} = nuevos;
end

end
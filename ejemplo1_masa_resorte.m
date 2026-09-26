% ejemplo1_masa_resorte - Hamiltonian mass-spring example
%
% Compares Adams-Bashforth, Adams-Moulton, implicit Taylor and implicit
% Runge-Kutta on the harmonic oscillator system
%
%       y1' = y2,
%       y2' = -(k/m)*y1,
%
% against the exact solution, in both state components and phase plane.
% Also builds error/order tables over a mesh-refinement sweep
% (n = 50, 100, ..., 1600) for both the solution and the energy H(y).
%
% Outputs (written to figDir = <script folder>/Imagenes):
%   - PNG figures: state components, phase plane, error/order tables
%   - CSV tables of error and observed order vs. n
%   - LaTeX .tex table files and a .tex figure/table fragment for the TFG
%
% No inputs; run as a script. Requires ABashforth1, AMoulton1,
% TaylorImplicito1, coefRKImplicito, RungeKuttaImplicito1 on the path.

% Ejemplo 1: sistema hamiltoniano masa-resorte.
% Aplica Adams-Bashforth, Adams-Moulton, Taylor implicito y
% Runge-Kutta implicito al sistema
%
%       y1' = y2,
%       y2' = -(k/m)y1,
%
% y compara cada aproximacion con la solucion exacta.

clear; close all; clc;

% Fondo blanco en todas las figuras y ejes.
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesColor','w');

% Colores oscuros para que ejes, numeracion, titulos y etiquetas
% se lean bien al insertar las figuras en LaTeX.
set(groot,'defaultAxesXColor','k');
set(groot,'defaultAxesYColor','k');
set(groot,'defaultTextColor','k');

%% 1. Cargar carpetas con las funciones
baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir)
    baseDir = pwd;
end

% Se incluyen dos posibles nombres porque, al descomprimir, algunos sistemas
% cambian la codificacion de las carpetas con acentos.
carpetasFunciones = { ...
    'MML', ...
    'RK implicito', 'RK implícito', 'RK impl#U00edcito', ...
    'Taylor implicito', 'Taylor implícito', 'Taylor impl#U00edcito'};

for j = 1:numel(carpetasFunciones)
    carpeta = fullfile(baseDir, carpetasFunciones{j});
    if isfolder(carpeta)
        addpath(carpeta);
    end
end

% Carpeta donde se guardaran las figuras.
figDir = fullfile(baseDir,'Imagenes');
if ~isfolder(figDir)
    mkdir(figDir);
end

%% 2. Datos del modelo masa-resorte
m = 1;          % masa
k = 1;          % constante elastica
omega = sqrt(k/m);

% Condicion inicial: posicion inicial x0 y velocidad inicial v0.
x0 = 1;
v0 = 0;
y0 = [x0; v0];

% Intervalo temporal y numero de pasos.
t0 = 0;
T  = 20*pi;
ta = [t0,T];
n  = 1000;

% Sistema escrito como y'=f(t,y), con y=(y1,y2)^T=(x,x')^T.
f = @(t,Y) [Y(2); -(k/m)*Y(1)];

% Solucion exacta del oscilador armonico.
solExacta = @(t) [x0*cos(omega*t) + (v0/omega)*sin(omega*t); ...
                 -x0*omega*sin(omega*t) + v0*cos(omega*t)];

% Energia exacta del sistema.
H = @(Y) 0.5*m*(Y(2,:).^2) + 0.5*k*(Y(1,:).^2);
H0 = H(y0);

%% 3. Parametros de los metodos
% Adams-Bashforth con w=3: metodo de 4 pasos.
wAB = 3;

% Adams-Moulton con w=3 y mu correcciones: predictor-corrector.
wAM = 3;
muAM = 4;

% Taylor implicito de grado 4.
gradoTaylor = 4;

% Runge-Kutta implicito de Gauss-Legendre con s=2 etapas, orden teorico 4.
sRK = 2;

% Parametros de Newton para los metodos implicitos.
tol = 1e-12;
maxit = 30;

%% 4. Ejecucion de los metodos
fprintf('Ejecutando Adams-Bashforth...\n');
[tAB,yAB,qAB,consAB] = ABashforth1(f,ta,y0,wAB,n);

fprintf('Ejecutando Adams-Moulton...\n');
[tAM,yAM,qAM,consAM] = AMoulton1(f,ta,y0,wAM,n,muAM);

fprintf('Ejecutando Taylor implicito...\n');
[tTI,yTI,qTI,consTI,infoTI] = TaylorImplicito1(f,ta,y0,gradoTaylor,n,tol,maxit);

fprintf('Construyendo coeficientes RK implicito de Gauss...\n');
[B,a,c,infoCoefRK] = coefRKImplicito(sRK,'gauss');

fprintf('Ejecutando Runge-Kutta implicito...\n');
[tRK,yRK,sRKsalida,consRK,infoRK] = RungeKuttaImplicito1(f,ta,y0,B,a,n,tol,maxit);

%% 5. Guardar resultados en una estructura comun
metodos = struct([]);
metodos(1).nombre = 'Adams-Bashforth';
metodos(1).nombreCorto = 'ABashforth';
metodos(1).t = tAB(:).';
metodos(1).y = yAB;
metodos(1).orden = qAB;
metodos(1).consistente = consAB;

metodos(2).nombre = 'Adams-Moulton';
metodos(2).nombreCorto = 'AMoulton';
metodos(2).t = tAM(:).';
metodos(2).y = yAM;
metodos(2).orden = qAM;
metodos(2).consistente = consAM;

metodos(3).nombre = 'Taylor implicito';
metodos(3).nombreCorto = 'TaylorImplicito';
metodos(3).t = tTI(:).';
metodos(3).y = yTI;
metodos(3).orden = qTI;
metodos(3).consistente = consTI;

metodos(4).nombre = 'Runge-Kutta implicito';
metodos(4).nombreCorto = 'RKImplicito';
metodos(4).t = tRK(:).';
metodos(4).y = yRK;
metodos(4).orden = infoCoefRK.ordenTeorico;
metodos(4).consistente = consRK;

%% 6. Graficas comparativas: componentes y plano de fases
% Se representan todas las aproximaciones en una misma figura.
% La solucion exacta se dibuja con linea continua y las aproximaciones
% numericas con linea discontinua y marcadores.
estilos = {'--o','--s','--d','--^'};

% Limites comunes para que todos los paneles de cada mosaico se
% representen exactamente con la misma escala. Para calcularlos se tienen
% en cuenta tanto la solucion exacta como todas las aproximaciones numericas.
tLimComun = [t0,T];
valoresComponentes = [];

for comp = 1:2
    for j = 1:numel(metodos)
        tAux = metodos(j).t;
        yAux = metodos(j).y;
        yExAux = solExacta(tAux);
        valoresComponentes = [valoresComponentes, yExAux(comp,:), yAux(comp,:)]; %#ok<AGROW>
    end
end
yLimComponentesComun = limitesConMargen(valoresComponentes,0.08);

valoresFase1 = [];
valoresFase2 = [];
for j = 1:numel(metodos)
    tAux = metodos(j).t;
    yAux = metodos(j).y;
    yExAux = solExacta(tAux);

    valoresFase1 = [valoresFase1, yExAux(1,:), yAux(1,:)]; %#ok<AGROW>
    valoresFase2 = [valoresFase2, yExAux(2,:), yAux(2,:)]; %#ok<AGROW>
end

limFaseComun = limitesConMargen([valoresFase1, valoresFase2],0.08);

% 6.1 Componentes y_1 e y_2.
% En lugar de representar todos los metodos en los mismos ejes,
% se generan dos mosaicos 2x2. En cada panel se compara la solucion
% exacta con un unico metodo, y el titulo indica el metodo empleado.
for comp = 1:2
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 650]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    for j = 1:numel(metodos)
        nexttile;
        hold on; grid on; box on;

        t = metodos(j).t;
        y = metodos(j).y;
        yEx = solExacta(t);
        idxMarca = unique(round(linspace(1,length(t),25)));

        plot(t,yEx(comp,:),'k-','LineWidth',1.8, ...
            'DisplayName','Solucion exacta');
        plot(t,y(comp,:),estilos{j},'LineWidth',1.1,'MarkerSize',4, ...
            'MarkerIndices',idxMarca,'DisplayName',metodos(j).nombre);

        xlabel('$t$','Interpreter','latex','Color','k');
        ylabel(sprintf('$y_%d(t)$',comp),'Interpreter','latex','Color','k');
        title(metodos(j).nombre,'Interpreter','latex','Color','k');

        legend('Location','best','TextColor','k','Color','w','EdgeColor','k');
        formatearEjes();
        xlim(tLimComun);
        ylim(yLimComponentesComun);
    end

    if comp == 1
        exportgraphics(fig,fullfile(figDir,'masa_resorte_mosaico_y1_metodos.png'),'Resolution',300);
    else
        exportgraphics(fig,fullfile(figDir,'masa_resorte_mosaico_y2_metodos.png'),'Resolution',300);
    end
    close(fig);
end

% 6.2 Plano de fases comparativo en una unica figura.
% Se vuelve a definir una malla de referencia para la solucion exacta.
tRef = metodos(1).t;
yExRef = solExacta(tRef);
fig = figure('Visible','off','Color','w');
hold on; grid on; box on;
set(gca,'Color','w');
plot(yExRef(1,:),yExRef(2,:),'k-','LineWidth',1.8,'DisplayName','Solucion exacta');
for j = 1:numel(metodos)
    y = metodos(j).y;
    idxMarca = unique(round(linspace(1,size(y,2),25)));
    plot(y(1,:),y(2,:),estilos{j},'LineWidth',1.1,'MarkerSize',4, ...
        'MarkerIndices',idxMarca,'DisplayName',metodos(j).nombre);
end
xlabel('$y_1$','Interpreter','latex','Color','k');
ylabel('$y_2$','Interpreter','latex','Color','k');
title('Plano de fases comparativo','Interpreter','latex','Color','k');
legend('Location','eastoutside','TextColor','k','Color','w','EdgeColor','k');
formatearEjes();
axis equal;
xlim(limFaseComun);
ylim(limFaseComun);
exportgraphics(fig,fullfile(figDir,'masa_resorte_fase_comparada.png'),'Resolution',300);
close(fig);

%% 7. Tablas de error y orden observado
% Para estudiar la dependencia del error con el tamano de paso, se repite
% el experimento con n = 50,100,...,1600. En cada caso se calcula:
%   - error frente a la solucion exacta: max_i ||y_i-y(t_i)||_infty,
%   - error relativo de energia: max_i |H(y_i)-H(y_0)|/|H(y_0)|.
% A partir de errores consecutivos se estima el orden observado mediante
%   p_n = log(E_{n_anterior}/E_n)/log(h_{n_anterior}/h_n).
% Como h=T/n, esta expresion equivale a
%   p_n = log(E_{n_anterior}/E_n)/log(n/n_anterior).

nGrid = [50, 100:100:1600];
numN = numel(nGrid);
numMetodos = numel(metodos);

errSolucion = zeros(numN,numMetodos);
errEnergia = zeros(numN,numMetodos);
ordenSolucion = NaN(numN,numMetodos);
ordenEnergia = NaN(numN,numMetodos);

fprintf('\nCalculando tablas de error para n=50:50:1600...\n');
for r = 1:numN
    nActual = nGrid(r);
    fprintf('  n = %d\n', nActual);

    [tABn,yABn] = ABashforth1(f,ta,y0,wAB,nActual);
    [tAMn,yAMn] = AMoulton1(f,ta,y0,wAM,nActual,muAM);
    [tTIn,yTIn] = TaylorImplicito1(f,ta,y0,gradoTaylor,nActual,tol,maxit);
    [tRKn,yRKn] = RungeKuttaImplicito1(f,ta,y0,B,a,nActual,tol,maxit);

    tLista = {tABn(:).', tAMn(:).', tTIn(:).', tRKn(:).'};
    yLista = {yABn, yAMn, yTIn, yRKn};

    for j = 1:numMetodos
        tAux = tLista{j};
        yAux = yLista{j};
        yExAux = solExacta(tAux);

        errSolucion(r,j) = max(vecnorm(yAux-yExAux,inf,1));
        errEnergia(r,j) = max(abs(H(yAux)-H0)/abs(H0));
    end

    if r > 1
        cocienteMallas = log(nGrid(r)/nGrid(r-1));
        for j = 1:numMetodos
            if errSolucion(r,j) > 0 && errSolucion(r-1,j) > 0
                ordenSolucion(r,j) = log(errSolucion(r-1,j)/errSolucion(r,j))/cocienteMallas;
            end
            if errEnergia(r,j) > 0 && errEnergia(r-1,j) > 0
                ordenEnergia(r,j) = log(errEnergia(r-1,j)/errEnergia(r,j))/cocienteMallas;
            end
        end
    end
end

% Guardar tablas en formato CSV, por si se quieren revisar fuera de LaTeX.
nombresVariables = {'Adams_Bashforth','Adams_Moulton','Taylor_implicito','RK_implicito'};

T_errSol = array2table(errSolucion,'VariableNames',nombresVariables);
T_errSol.n = nGrid(:);
T_errSol = movevars(T_errSol,'n','Before',1);
writetable(T_errSol,fullfile(figDir,'tabla_masa_resorte_error_solucion.csv'));

T_ordSol = array2table(ordenSolucion,'VariableNames',nombresVariables);
T_ordSol.n = nGrid(:);
T_ordSol = movevars(T_ordSol,'n','Before',1);
writetable(T_ordSol,fullfile(figDir,'tabla_masa_resorte_orden_solucion.csv'));

T_errEnergia = array2table(errEnergia,'VariableNames',nombresVariables);
T_errEnergia.n = nGrid(:);
T_errEnergia = movevars(T_errEnergia,'n','Before',1);
writetable(T_errEnergia,fullfile(figDir,'tabla_masa_resorte_error_energia.csv'));

T_ordEnergia = array2table(ordenEnergia,'VariableNames',nombresVariables);
T_ordEnergia.n = nGrid(:);
T_ordEnergia = movevars(T_ordEnergia,'n','Before',1);
writetable(T_ordEnergia,fullfile(figDir,'tabla_masa_resorte_orden_energia.csv'));

% Guardar las mismas tablas directamente en formato LaTeX.
nombresCabecera = {'AB','AM','Taylor imp.','RK imp.'};
escribirTablaLatex(fullfile(figDir,'tabla_masa_resorte_error_solucion.tex'), ...
    'Errores maximos frente a la solucion exacta para el sistema masa-resorte.', ...
    'tab:masa_resorte_error_solucion', nGrid, errSolucion, nombresCabecera, 'sci');

escribirTablaLatex(fullfile(figDir,'tabla_masa_resorte_orden_solucion.tex'), ...
    'Orden observado a partir de los errores frente a la solucion exacta para el sistema masa-resorte.', ...
    'tab:masa_resorte_orden_solucion', nGrid, ordenSolucion, nombresCabecera, 'dec');

escribirTablaLatex(fullfile(figDir,'tabla_masa_resorte_error_energia.tex'), ...
    'Errores relativos maximos de energia para el sistema masa-resorte.', ...
    'tab:masa_resorte_error_energia', nGrid, errEnergia, nombresCabecera, 'sci');

escribirTablaLatex(fullfile(figDir,'tabla_masa_resorte_orden_energia.tex'), ...
    'Orden observado a partir de los errores relativos de energia para el sistema masa-resorte.', ...
    'tab:masa_resorte_orden_energia', nGrid, ordenEnergia, nombresCabecera, 'dec');

%% 8. Tabla resumen por pantalla
fprintf('\nResumen del experimento:\n');
fprintf('Metodo\t\t\tOrden/parametro\tConsistente\tError final\tError energia final\n');
for j = 1:numel(metodos)
    t = metodos(j).t;
    y = metodos(j).y;
    yEx = solExacta(t);
    errFinal = norm(y(:,end)-yEx(:,end),inf);
    errHFinal = abs(H(y(:,end))-H0)/abs(H0);

    fprintf('%-22s\t%-8g\t\t%d\t\t%.3e\t%.3e\n', ...
        metodos(j).nombre, metodos(j).orden, metodos(j).consistente, errFinal, errHFinal);
end

%% 9. Generar fragmento LaTeX para insertar las figuras y tablas en el TFG
texFile = fullfile(figDir,'fragmento_masa_resorte.tex');
fid = fopen(texFile,'w');

fprintf(fid,'%% Fragmento generado automaticamente por ejemplo1_masa_resorte.m\n');
fprintf(fid,'%% Para las tablas se requiere \\usepackage{longtable}.\n\n');

fprintf(fid,'\\begin{figure}[!htbp]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.95\\textwidth,height=0.39\\textheight,keepaspectratio]{Imagenes/masa_resorte_mosaico_y1_metodos.png}\n');
fprintf(fid,'    \\vspace{0.5em}\n');
fprintf(fid,'    \\includegraphics[width=0.95\\textwidth,height=0.39\\textheight,keepaspectratio]{Imagenes/masa_resorte_mosaico_y2_metodos.png}\n');
fprintf(fid,'    \\caption{Comparacion de las componentes del sistema masa-resorte para los metodos descritos en la Tabla \\ref{tab:parametros_masa_resorte}. En cada mosaico, cada panel compara la solucion exacta con uno de los metodos numericos. La solucion exacta se representa con linea continua y la aproximacion numerica con linea discontinua y marcadores.}\n');
fprintf(fid,'    \\label{fig:masa_resorte_componentes}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\begin{figure}[!htbp]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.72\\textwidth]{Imagenes/masa_resorte_fase_comparada.png}\n');
fprintf(fid,'    \\caption{Plano de fases del sistema masa-resorte para los metodos descritos en la Tabla \\ref{tab:parametros_masa_resorte}. La solucion exacta se representa con linea continua y las aproximaciones numericas con lineas discontinuas y marcadores.}\n');
fprintf(fid,'    \\label{fig:masa_resorte_fase_comparada}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\input{Imagenes/tabla_masa_resorte_error_solucion.tex}\n\n');
fprintf(fid,'\\input{Imagenes/tabla_masa_resorte_orden_solucion.tex}\n\n');
fprintf(fid,'\\input{Imagenes/tabla_masa_resorte_error_energia.tex}\n\n');
fprintf(fid,'\\input{Imagenes/tabla_masa_resorte_orden_energia.tex}\n');

fclose(fid);

fprintf('\nFiguras guardadas en:\n%s\n', figDir);
fprintf('Fragmento LaTeX guardado en:\n%s\n', texFile);


function escribirTablaLatex(nombreArchivo, captionTabla, labelTabla, nGrid, valores, nombresCabecera, tipoFormato)
% escribirTablaLatex - Write a results matrix as a LaTeX table
%
% Inputs:
%   nombreArchivo   - output .tex file path
%   captionTabla    - table caption text
%   labelTabla      - LaTeX \label key
%   nGrid           - vector of mesh sizes n (first column)
%   valores         - numMesh x numMetodos matrix of values to tabulate
%   nombresCabecera - cell array of column header names (one per method)
%   tipoFormato     - 'sci' for %.3e formatting, anything else for %.3f
%
% No output. NaN/Inf entries are written as '--'.

%ESCRIBIRTABLALATEX Escribe una matriz de resultados como tabla LaTeX.
% La primera columna es n y las columnas restantes corresponden a los metodos.

fid = fopen(nombreArchivo,'w');
if fid == -1
    error('No se pudo crear el archivo %s.', nombreArchivo);
end

numMetodos = numel(nombresCabecera);
alineacion = ['|c|' repmat('c|',1,numMetodos)];

fprintf(fid,'\\begin{table}[!htbp]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\small\n');
fprintf(fid,'\\setlength{\\tabcolsep}{4pt}\n');
fprintf(fid,'\\begin{tabular}{%s}\n', alineacion);

fprintf(fid,'\\hline\n');
fprintf(fid,'$n$');
for j = 1:numMetodos
    fprintf(fid,' & %s', nombresCabecera{j});
end
fprintf(fid,'\\\\\n');
fprintf(fid,'\\hline\n');

for r = 1:numel(nGrid)
    fprintf(fid,'%d', nGrid(r));
    for j = 1:numMetodos
        valor = valores(r,j);

        if isnan(valor) || isinf(valor)
            cadena = '--';
        elseif strcmp(tipoFormato,'sci')
            cadena = sprintf('%.3e', valor);
        else
            cadena = sprintf('%.3f', valor);
        end

        fprintf(fid,' & %s', cadena);
    end
    fprintf(fid,'\\\\\n');
    fprintf(fid,'\\hline\n');
end

fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\caption{%s}\n', captionTabla);
fprintf(fid,'\\label{%s}\n', labelTabla);
fprintf(fid,'\\end{table}\n');

fclose(fid);
end

function formatearEjes()
% formatearEjes - Darken axes, ticks, labels and grid for figure export
%
% No inputs/outputs. Operates on the current axes (gca) so that exported
% figures read well when embedded in LaTeX with a white background.

%FORMATEAREJES Oscurece ejes, marcas, etiquetas y rejilla para exportar figuras.
ax = gca;
set(ax,'Color','w', ...
    'XColor','k','YColor','k', ...
    'GridColor',[0.25 0.25 0.25], ...
    'MinorGridColor',[0.45 0.45 0.45], ...
    'GridAlpha',0.35, ...
    'MinorGridAlpha',0.25, ...
    'LineWidth',0.8, ...
    'FontSize',10);
ax.Title.Color = 'k';
ax.XLabel.Color = 'k';
ax.YLabel.Color = 'k';
end



function lim = limitesConMargen(valores,margenRelativo)
% limitesConMargen - Common axis limits with a relative margin
%
% Inputs:
%   valores        - vector of values the limits must cover
%   margenRelativo - fractional margin added on each side (e.g. 0.08)
%
% Output:
%   lim - [min,max] limits, expanded by margenRelativo*range(valores).
%         If the data is constant, an absolute margin is used instead
%         to avoid degenerate limits; returns [-1,1] if valores is empty.

%LIMITESCONMARGEN Devuelve limites comunes con un margen relativo.
% Si los datos son constantes, se anade un margen absoluto para evitar
% limites degenerados.
valores = valores(:);
valores = valores(isfinite(valores));

if isempty(valores)
    lim = [-1,1];
    return;
end

minValor = min(valores);
maxValor = max(valores);
rango = maxValor - minValor;

if rango == 0
    margen = max(1,abs(minValor))*margenRelativo;
else
    margen = margenRelativo*rango;
end

lim = [minValor - margen, maxValor + margen];
end

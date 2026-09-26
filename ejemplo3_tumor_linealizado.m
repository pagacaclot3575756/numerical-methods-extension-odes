%% EJEMPLO 3: TUMOR LINEALIZADO PROLIFERANTES-QUIESCENTES
% ejemplo3_tumor_linealizado - Linearized proliferating/quiescent tumor model
%
% Applies Adams-Bashforth, Adams-Moulton, implicit Taylor and implicit
% Runge-Kutta to the two-population linear system
%
%       P'(t) = (beta - muP - r0)*P(t) + r1*Q(t),
%       Q'(t) = r0*P(t) - (r1 + muQ)*Q(t),
%
% where P(t) is the proliferating-cell population and Q(t) the
% quiescent-cell population, and compares each method against the
% reference solution obtained via the matrix exponential. Also tracks
% the proliferating fraction G = P/(P+Q) and builds error/order tables
% over a mesh-refinement sweep.
%
% Outputs (written to figDir = <script folder>/Imagenes):
%   - PNG figures: population components, proliferating-fraction curves,
%     error/order behavior
%   - CSV and LaTeX .tex tables of error and observed order vs. n
%   - A LaTeX .tex fragment with \includegraphics/\input blocks for the TFG
%
% No inputs; run as a script. Requires ABashforth1, AMoulton1,
% TaylorImplicito1, coefRKImplicito, RungeKuttaImplicito1 on the path.

% Modelo lineal de dos poblaciones celulares:
%
%       P'(t) = (beta - muP - r0) P(t) + r1 Q(t),
%       Q'(t) = r0 P(t) - (r1 + muQ) Q(t),
%
% donde P(t) representa celulas proliferantes y Q(t) celulas quiescentes.
% Se compara cada metodo numerico con la solucion de referencia obtenida
% mediante la exponencial matricial.

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

% Se incluyen varios posibles nombres porque, al descomprimir, algunos
% sistemas cambian la codificacion de las carpetas con acentos.
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

% Carpeta donde se guardaran las figuras y tablas.
figDir = fullfile(baseDir,'Imagenes');
if ~isfolder(figDir)
    mkdir(figDir);
end

%% 2. Datos del modelo tumoral linealizado

% beta  : tasa de division de las celulas proliferantes.
% muP   : tasa de mortalidad de las celulas proliferantes.
% muQ   : tasa de mortalidad de las celulas quiescentes.
% r0    : tasa de paso de proliferantes a quiescentes.
% r1    : tasa de paso de quiescentes a proliferantes.

beta = 0.16;
muP  = 0.03;
muQ  = 0.04;
r0   = 0.08;
r1   = 0.03;

% Matriz del sistema y' = A y, con y=(P,Q)^T.
A = [ beta - muP - r0,       r1;
                 r0, -(r1 + muQ) ];

f = @(t,Y) A*Y;

% Condicion inicial: poblaciones relativas iniciales.
P0 = 1.0;
Q0 = 0.2;
y0 = [P0; Q0];

% Intervalo temporal y numero de pasos.
t0 = 0;
T  = 60;
ta = [t0,T];
n  = 1000;

%% 3. Solucion de referencia

% Como el sistema es lineal, la solucion exacta se obtiene con expm(A(t-t0)).
% Se define como funcion para usar la misma referencia en graficas y tablas.
solExacta = @(t) solucionReferenciaLineal(A,y0,t,t0);

%% 4. Parametros de los metodos

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

%% 5. Ejecucion de los metodos

fprintf('Ejecutando Adams-Bashforth...\n');
[tAB,yAB,qAB,consAB] = ABashforth1(f,ta,y0,wAB,n);

fprintf('Ejecutando Adams-Moulton...\n');
[tAM,yAM,qAM,consAM] = AMoulton1(f,ta,y0,wAM,n,muAM);

fprintf('Ejecutando Taylor implicito...\n');
[tTI,yTI,qTI,consTI,infoTI] = TaylorImplicito1(f,ta,y0,gradoTaylor,n,tol,maxit); %#ok<ASGLU>

fprintf('Construyendo coeficientes RK implicito de Gauss...\n');
[B,a,cRKcoef,infoCoefRK] = coefRKImplicito(sRK,'gauss'); %#ok<ASGLU>

fprintf('Ejecutando Runge-Kutta implicito...\n');
[tRK,yRK,sRKsalida,consRK,infoRK] = RungeKuttaImplicito1(f,ta,y0,B,a,n,tol,maxit); %#ok<ASGLU>

%% 6. Guardar resultados en una estructura comun

metodos = struct([]);

metodos(1).nombre = 'Adams-Bashforth';
metodos(1).nombreCorto = 'ABashforth';
metodos(1).t = tAB(:).';
metodos(1).y = orientarSolucion(yAB,y0,tAB);
metodos(1).orden = qAB;
metodos(1).consistente = consAB;

metodos(2).nombre = 'Adams-Moulton';
metodos(2).nombreCorto = 'AMoulton';
metodos(2).t = tAM(:).';
metodos(2).y = orientarSolucion(yAM,y0,tAM);
metodos(2).orden = qAM;
metodos(2).consistente = consAM;

metodos(3).nombre = 'Taylor implicito';
metodos(3).nombreCorto = 'TaylorImplicito';
metodos(3).t = tTI(:).';
metodos(3).y = orientarSolucion(yTI,y0,tTI);
metodos(3).orden = qTI;
metodos(3).consistente = consTI;

metodos(4).nombre = 'Runge-Kutta implicito';
metodos(4).nombreCorto = 'RKImplicito';
metodos(4).t = tRK(:).';
metodos(4).y = orientarSolucion(yRK,y0,tRK);
metodos(4).orden = infoCoefRK.ordenTeorico;
metodos(4).consistente = consRK;

%% 7. Graficas comparativas con escala comun
% Como en el ejemplo 1, se evita representar todos los metodos sobre los
% mismos ejes. En su lugar, se generan mosaicos 2x2. En cada panel se
% compara la solucion de referencia con un unico metodo numerico, y todos
% los paneles de una misma magnitud usan exactamente los mismos limites.

estilos = {'--o','--s','--d','--^'};
tLimComun = [t0,T];

valoresPoblaciones = [];
valoresTotales = [];
valoresFracciones = [];
valoresErrores = [];

for j = 1:numel(metodos)
    t = metodos(j).t;
    Y = metodos(j).y;
    Yex = solExacta(t);

    Nnum = sum(Y,1);
    Nex  = sum(Yex,1);

    Gnum = fraccionProliferante(Y);
    Gex  = fraccionProliferante(Yex);

    err = vecnorm(Y-Yex,inf,1);

    valoresPoblaciones = [valoresPoblaciones, Yex(1,:), Yex(2,:), Y(1,:), Y(2,:)]; %#ok<AGROW>
    valoresTotales = [valoresTotales, Nex, Nnum]; %#ok<AGROW>
    valoresFracciones = [valoresFracciones, Gex, Gnum]; %#ok<AGROW>
    valoresErrores = [valoresErrores, max(err,eps)]; %#ok<AGROW>
end

yLimPoblacionesComun = limitesConMargen(valoresPoblaciones,0.08);
yLimTotalComun = limitesConMargen(valoresTotales,0.08);
yLimFraccionComun = limitesConMargen(valoresFracciones,0.08);
yLimErrorComun = limitesLogConMargen(valoresErrores,0.08);

% 7.1 Poblaciones P y Q.
for comp = 1:2
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 650]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    for j = 1:numel(metodos)
        nexttile;
        hold on; grid on; box on;

        t = metodos(j).t;
        Y = metodos(j).y;
        Yex = solExacta(t);
        idxMarca = unique(round(linspace(1,length(t),25)));

        plot(t,Yex(comp,:),'k-','LineWidth',1.8, ...
            'DisplayName','Solucion de referencia');
        plot(t,Y(comp,:),estilos{j},'LineWidth',1.1,'MarkerSize',4, ...
            'MarkerIndices',idxMarca,'DisplayName',metodos(j).nombre);

        xlabel('$t$','Interpreter','latex','Color','k');
        if comp == 1
            ylabel('$P(t)$','Interpreter','latex','Color','k');
            title([metodos(j).nombre, ': proliferantes'],'Interpreter','latex','Color','k');
        else
            ylabel('$Q(t)$','Interpreter','latex','Color','k');
            title([metodos(j).nombre, ': quiescentes'],'Interpreter','latex','Color','k');
        end

        legend('Location','best','TextColor','k','Color','w','EdgeColor','k');
        formatearEjes();
        xlim(tLimComun);
        ylim(yLimPoblacionesComun);
    end

    if comp == 1
        exportgraphics(fig,fullfile(figDir,'tumor_mosaico_P_metodos.png'),'Resolution',300);
    else
        exportgraphics(fig,fullfile(figDir,'tumor_mosaico_Q_metodos.png'),'Resolution',300);
    end
    close(fig);
end

% 7.2 Poblacion tumoral total.
fig = figure('Visible','off','Color','w','Position',[100 100 1100 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:numel(metodos)
    nexttile;
    hold on; grid on; box on;

    t = metodos(j).t;
    Y = metodos(j).y;
    Yex = solExacta(t);
    idxMarca = unique(round(linspace(1,length(t),25)));

    Nnum = sum(Y,1);
    Nex  = sum(Yex,1);

    plot(t,Nex,'k-','LineWidth',1.8, ...
        'DisplayName','Solucion de referencia');
    plot(t,Nnum,estilos{j},'LineWidth',1.1,'MarkerSize',4, ...
        'MarkerIndices',idxMarca,'DisplayName',metodos(j).nombre);

    xlabel('$t$','Interpreter','latex','Color','k');
    ylabel('$N(t)=P(t)+Q(t)$','Interpreter','latex','Color','k');
    title(metodos(j).nombre,'Interpreter','latex','Color','k');
    legend('Location','best','TextColor','k','Color','w','EdgeColor','k');
    formatearEjes();
    xlim(tLimComun);
    ylim(yLimTotalComun);
end

exportgraphics(fig,fullfile(figDir,'tumor_mosaico_total_metodos.png'),'Resolution',300);
close(fig);

% 7.3 Fraccion de celulas proliferantes.
fig = figure('Visible','off','Color','w','Position',[100 100 1100 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:numel(metodos)
    nexttile;
    hold on; grid on; box on;

    t = metodos(j).t;
    Y = metodos(j).y;
    Yex = solExacta(t);
    idxMarca = unique(round(linspace(1,length(t),25)));

    Gnum = fraccionProliferante(Y);
    Gex  = fraccionProliferante(Yex);

    plot(t,Gex,'k-','LineWidth',1.8, ...
        'DisplayName','Solucion de referencia');
    plot(t,Gnum,estilos{j},'LineWidth',1.1,'MarkerSize',4, ...
        'MarkerIndices',idxMarca,'DisplayName',metodos(j).nombre);

    xlabel('$t$','Interpreter','latex','Color','k');
    ylabel('$P(t)/(P(t)+Q(t))$','Interpreter','latex','Color','k');
    title(metodos(j).nombre,'Interpreter','latex','Color','k');
    legend('Location','best','TextColor','k','Color','w','EdgeColor','k');
    formatearEjes();
    xlim(tLimComun);
    ylim(yLimFraccionComun);
end

exportgraphics(fig,fullfile(figDir,'tumor_mosaico_fraccion_metodos.png'),'Resolution',300);
close(fig);

% 7.4 Error infinito frente a la solucion de referencia.
fig = figure('Visible','off','Color','w','Position',[100 100 1100 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:numel(metodos)
    nexttile;
    hold on; grid on; box on;

    t = metodos(j).t;
    Y = metodos(j).y;
    Yex = solExacta(t);
    err = max(vecnorm(Y-Yex,inf,1),eps);

    semilogy(t,err,'LineWidth',1.2,'DisplayName',metodos(j).nombre);

    xlabel('$t$','Interpreter','latex','Color','k');
    ylabel('$\|y_i-y(t_i)\|_{\infty}$','Interpreter','latex','Color','k');
    title(metodos(j).nombre,'Interpreter','latex','Color','k');
    legend('Location','best','TextColor','k','Color','w','EdgeColor','k');
    formatearEjes();
    xlim(tLimComun);
    ylim(yLimErrorComun);
end

exportgraphics(fig,fullfile(figDir,'tumor_mosaico_error_global_metodos.png'),'Resolution',300);
close(fig);

%% 8. Tablas de error y orden observado
% Para estudiar la dependencia del error con el tamano de paso, se repite
% el experimento con n = 50,100,...,1600. En cada caso se calcula:
%   - error frente a la solucion exacta obtenida mediante la exponencial matricial:
%     max_i ||y_i-y(t_i)||_infty,
%   - error en la poblacion total: max_i |N_i-N(t_i)|,
%   - error en la fraccion proliferante: max_i |G_i-G(t_i)|.
% A partir de errores consecutivos se estima el orden observado mediante
%   p_n = log(E_{n_anterior}/E_n)/log(n/n_anterior),
% ya que h=T/n.

nGrid = [50, 100:100:1600];
numN = numel(nGrid);
numMetodos = numel(metodos);

errSolucion = zeros(numN,numMetodos);
errTotal = zeros(numN,numMetodos);
errFraccion = zeros(numN,numMetodos);
ordenSolucion = NaN(numN,numMetodos);
ordenTotal = NaN(numN,numMetodos);
ordenFraccion = NaN(numN,numMetodos);

fprintf('\nCalculando tablas de error para n=50,100,...,1600...\n');
for r = 1:numN
    nActual = nGrid(r);
    fprintf('  n = %d\n', nActual);

    [tABn,yABn] = ABashforth1(f,ta,y0,wAB,nActual);
    [tAMn,yAMn] = AMoulton1(f,ta,y0,wAM,nActual,muAM);
    [tTIn,yTIn] = TaylorImplicito1(f,ta,y0,gradoTaylor,nActual,tol,maxit);
    [tRKn,yRKn] = RungeKuttaImplicito1(f,ta,y0,B,a,nActual,tol,maxit);

    tLista = {tABn(:).', tAMn(:).', tTIn(:).', tRKn(:).'};
    yLista = { ...
        orientarSolucion(yABn,y0,tABn), ...
        orientarSolucion(yAMn,y0,tAMn), ...
        orientarSolucion(yTIn,y0,tTIn), ...
        orientarSolucion(yRKn,y0,tRKn)};

    for j = 1:numMetodos
        tAux = tLista{j};
        yAux = yLista{j};
        yExAux = solExacta(tAux);

        Nnum = sum(yAux,1);
        Nex  = sum(yExAux,1);

        Gnum = fraccionProliferante(yAux);
        Gex  = fraccionProliferante(yExAux);

        errSolucion(r,j) = max(vecnorm(yAux-yExAux,inf,1));
        errTotal(r,j) = max(abs(Nnum-Nex));
        errFraccion(r,j) = max(abs(Gnum-Gex));
    end

    if r > 1
        cocienteMallas = log(nGrid(r)/nGrid(r-1));
        for j = 1:numMetodos
            if errSolucion(r,j) > 0 && errSolucion(r-1,j) > 0
                ordenSolucion(r,j) = log(errSolucion(r-1,j)/errSolucion(r,j))/cocienteMallas;
            end
            if errTotal(r,j) > 0 && errTotal(r-1,j) > 0
                ordenTotal(r,j) = log(errTotal(r-1,j)/errTotal(r,j))/cocienteMallas;
            end
            if errFraccion(r,j) > 0 && errFraccion(r-1,j) > 0
                ordenFraccion(r,j) = log(errFraccion(r-1,j)/errFraccion(r,j))/cocienteMallas;
            end
        end
    end
end

% Guardar tablas en formato CSV, por si se quieren revisar fuera de LaTeX.
nombresVariables = {'Adams_Bashforth','Adams_Moulton','Taylor_implicito','RK_implicito'};

T_errSol = array2table(errSolucion,'VariableNames',nombresVariables);
T_errSol.n = nGrid(:);
T_errSol = movevars(T_errSol,'n','Before',1);
writetable(T_errSol,fullfile(figDir,'tabla_tumor_error_solucion.csv'));

T_ordSol = array2table(ordenSolucion,'VariableNames',nombresVariables);
T_ordSol.n = nGrid(:);
T_ordSol = movevars(T_ordSol,'n','Before',1);
writetable(T_ordSol,fullfile(figDir,'tabla_tumor_orden_solucion.csv'));

T_errTotal = array2table(errTotal,'VariableNames',nombresVariables);
T_errTotal.n = nGrid(:);
T_errTotal = movevars(T_errTotal,'n','Before',1);
writetable(T_errTotal,fullfile(figDir,'tabla_tumor_error_total.csv'));

T_ordTotal = array2table(ordenTotal,'VariableNames',nombresVariables);
T_ordTotal.n = nGrid(:);
T_ordTotal = movevars(T_ordTotal,'n','Before',1);
writetable(T_ordTotal,fullfile(figDir,'tabla_tumor_orden_total.csv'));

T_errFraccion = array2table(errFraccion,'VariableNames',nombresVariables);
T_errFraccion.n = nGrid(:);
T_errFraccion = movevars(T_errFraccion,'n','Before',1);
writetable(T_errFraccion,fullfile(figDir,'tabla_tumor_error_fraccion.csv'));

T_ordFraccion = array2table(ordenFraccion,'VariableNames',nombresVariables);
T_ordFraccion.n = nGrid(:);
T_ordFraccion = movevars(T_ordFraccion,'n','Before',1);
writetable(T_ordFraccion,fullfile(figDir,'tabla_tumor_orden_fraccion.csv'));

% Guardar las mismas tablas directamente en formato LaTeX.
nombresCabecera = {'AB','AM','Taylor imp.','RK imp.'};
escribirTablaLatex(fullfile(figDir,'tabla_tumor_error_solucion.tex'), ...
    'Errores maximos frente a la solucion de referencia para el modelo tumoral linealizado.', ...
    'tab:tumor_error_solucion', nGrid, errSolucion, nombresCabecera, 'sci');

escribirTablaLatex(fullfile(figDir,'tabla_tumor_orden_solucion.tex'), ...
    'Orden observado a partir de los errores frente a la solucion de referencia para el modelo tumoral linealizado.', ...
    'tab:tumor_orden_solucion', nGrid, ordenSolucion, nombresCabecera, 'dec');

escribirTablaLatex(fullfile(figDir,'tabla_tumor_error_total.tex'), ...
    'Errores maximos en la poblacion tumoral total para el modelo tumoral linealizado.', ...
    'tab:tumor_error_total', nGrid, errTotal, nombresCabecera, 'sci');

escribirTablaLatex(fullfile(figDir,'tabla_tumor_orden_total.tex'), ...
    'Orden observado a partir de los errores en la poblacion tumoral total para el modelo tumoral linealizado.', ...
    'tab:tumor_orden_total', nGrid, ordenTotal, nombresCabecera, 'dec');

escribirTablaLatex(fullfile(figDir,'tabla_tumor_error_fraccion.tex'), ...
    'Errores maximos en la fraccion de celulas proliferantes para el modelo tumoral linealizado.', ...
    'tab:tumor_error_fraccion', nGrid, errFraccion, nombresCabecera, 'sci');

escribirTablaLatex(fullfile(figDir,'tabla_tumor_orden_fraccion.tex'), ...
    'Orden observado a partir de los errores en la fraccion de celulas proliferantes para el modelo tumoral linealizado.', ...
    'tab:tumor_orden_fraccion', nGrid, ordenFraccion, nombresCabecera, 'dec');

%% 9. Tabla resumen por pantalla

fprintf('\nResumen del experimento:\n');
fprintf('Metodo\t\t\tOrden/parametro\tConsistente\tError final\tError total final\tError fraccion final\n');

for j = 1:numel(metodos)
    t = metodos(j).t;
    Y = metodos(j).y;
    Yex = solExacta(t);

    errFinal = norm(Y(:,end)-Yex(:,end),inf);

    Nnum = sum(Y,1);
    Nex  = sum(Yex,1);
    Gnum = fraccionProliferante(Y);
    Gex  = fraccionProliferante(Yex);

    errTotalFinal = abs(Nnum(end)-Nex(end));
    errFraccionFinal = abs(Gnum(end)-Gex(end));

    fprintf('%-22s\t%-8g\t\t%d\t\t%.3e\t%.3e\t%.3e\n', ...
        metodos(j).nombre, metodos(j).orden, metodos(j).consistente, ...
        errFinal, errTotalFinal, errFraccionFinal);
end

%% 10. Generar fragmento LaTeX para insertar las figuras y tablas en el TFG

texFile = fullfile(figDir,'fragmento_tumor_linealizado.tex');
fid = fopen(texFile,'w');

fprintf(fid,'%% Fragmento generado automaticamente por ejemplo3_tumor_linealizado.m\n');
fprintf(fid,'%% Para las tablas se requiere \\usepackage{longtable}.\n\n');

fprintf(fid,'\\begin{figure}[!htbp]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.95\\textwidth,height=0.39\\textheight,keepaspectratio]{Imagenes/tumor_mosaico_P_metodos.png}\n');
fprintf(fid,'    \\vspace{0.5em}\n');
fprintf(fid,'    \\includegraphics[width=0.95\\textwidth,height=0.39\\textheight,keepaspectratio]{Imagenes/tumor_mosaico_Q_metodos.png}\n');
fprintf(fid,'    \\caption{Comparacion de las poblaciones proliferantes y quiescentes del modelo tumoral linealizado. En cada mosaico, cada panel compara la solucion de referencia con uno de los metodos numericos. La solucion de referencia se representa con linea continua y la aproximacion numerica con linea discontinua y marcadores.}\n');
fprintf(fid,'    \\label{fig:tumor_poblaciones}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\begin{figure}[!htbp]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.95\\textwidth,height=0.40\\textheight,keepaspectratio]{Imagenes/tumor_mosaico_total_metodos.png}\n');
fprintf(fid,'    \\caption{Evolucion de la poblacion tumoral total en el modelo linealizado. Cada panel compara la solucion de referencia con uno de los metodos numericos, manteniendo la misma escala en todos los paneles.}\n');
fprintf(fid,'    \\label{fig:tumor_total}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\begin{figure}[!htbp]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.95\\textwidth,height=0.40\\textheight,keepaspectratio]{Imagenes/tumor_mosaico_fraccion_metodos.png}\n');
fprintf(fid,'    \\caption{Evolucion de la fraccion de celulas proliferantes en el modelo tumoral linealizado. Cada panel compara la solucion de referencia con uno de los metodos numericos, manteniendo la misma escala en todos los paneles.}\n');
fprintf(fid,'    \\label{fig:tumor_fraccion}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\begin{figure}[!htbp]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.95\\textwidth,height=0.40\\textheight,keepaspectratio]{Imagenes/tumor_mosaico_error_global_metodos.png}\n');
fprintf(fid,'    \\caption{Error global de los metodos numericos frente a la solucion de referencia del modelo tumoral linealizado. Cada panel representa un unico metodo con la misma escala logaritmica vertical.}\n');
fprintf(fid,'    \\label{fig:tumor_error_global}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\input{Imagenes/tabla_tumor_error_solucion.tex}\n\n');
fprintf(fid,'\\input{Imagenes/tabla_tumor_orden_solucion.tex}\n\n');
fprintf(fid,'\\input{Imagenes/tabla_tumor_error_total.tex}\n\n');
fprintf(fid,'\\input{Imagenes/tabla_tumor_orden_total.tex}\n\n');
fprintf(fid,'\\input{Imagenes/tabla_tumor_error_fraccion.tex}\n\n');
fprintf(fid,'\\input{Imagenes/tabla_tumor_orden_fraccion.tex}\n');

fclose(fid);

fprintf('\nFiguras y tablas guardadas en:\n%s\n', figDir);
fprintf('Fragmento LaTeX guardado en:\n%s\n', texFile);

%% FUNCIONES AUXILIARES

function Y = solucionReferenciaLineal(A,y0,t,t0)
% solucionReferenciaLineal - Exact solution of the linear system y'=Ay
%
% Inputs:
%   A  - d x d system matrix
%   y0 - initial condition (column vector, length d)
%   t  - row/column vector of times at which to evaluate the solution
%   t0 - initial time
%
% Output:
%   Y - d x length(t) matrix, Y(:,i) = expm(A*(t(i)-t0))*y0

%SOLUCIONREFERENCIALINEAL Calcula la solucion exacta de y'=Ay.
% La solucion se evalua en los puntos de la malla t mediante expm(A(t-t0)).

    t = t(:).';
    m = length(y0);
    Y = zeros(m,length(t));

    for i = 1:length(t)
        Y(:,i) = expm(A*(t(i)-t0))*y0;
    end
end

function Y = orientarSolucion(Y,y0,t)
% orientarSolucion - Ensure a solution matrix is oriented variables x time
%
% Inputs:
%   Y  - solution matrix, either m x Nt or Nt x m
%   y0 - initial condition, used only to get the system dimension m
%   t  - time vector, used only to get Nt
%
% Output:
%   Y - the same matrix, transposed to m x Nt if needed. Errors if
%       neither orientation matches size(Y).

%ORIENTARSOLUCION Asegura el formato variables x tiempos.

    m = length(y0);
    Nt = length(t);

    if isequal(size(Y),[m,Nt])
        return

    elseif isequal(size(Y),[Nt,m])
        Y = Y.';

    else
        error(['La solucion tiene tamano inesperado. ', ...
               'Se esperaba %d x %d o %d x %d, pero se obtuvo %d x %d.'], ...
               m,Nt,Nt,m,size(Y,1),size(Y,2));
    end
end

function G = fraccionProliferante(Y)
% fraccionProliferante - Proliferating-cell fraction P/(P+Q)
%
% Input:
%   Y - 2 x Nt matrix with Y(1,:)=P(t) and Y(2,:)=Q(t)
%
% Output:
%   G - 1 x Nt row vector, G = P./(P+Q). Denominators with
%       |P+Q| < eps are replaced by eps to avoid division by zero.

%FRACCIONPROLIFERANTE Calcula P/(P+Q) evitando divisiones numericamente nulas.

    denominador = sum(Y,1);
    mascara = abs(denominador) < eps;

    if any(mascara)
        denominador(mascara) = eps;
    end

    G = Y(1,:)./denominador;
end

function escribirTablaLatex(nombreArchivo, captionTabla, labelTabla, nGrid, valores, nombresCabecera, tipoFormato)
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
% limitesLogConMargen - Common positive axis limits for semilogy plots
%
% Inputs:
%   valores        - vector of (expected positive) values to cover
%   margenRelativo - fractional margin added on each side, in log10 scale
%
% Output:
%   lim - [min,max] limits (in linear scale, i.e. 10.^[...]), expanded
%         by margenRelativo*range(log10(valores)). Non-positive and
%         non-finite entries are discarded before computing the range;
%         returns [1e-16,1] if nothing usable remains.

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

function lim = limitesLogConMargen(valores,margenRelativo)
%LIMITESLOGCONMARGEN Devuelve limites positivos comunes para graficas semilogy.
valores = valores(:);
valores = valores(isfinite(valores) & valores > 0);

if isempty(valores)
    lim = [1e-16,1];
    return;
end

logValores = log10(valores);
minValor = min(logValores);
maxValor = max(logValores);
rango = maxValor - minValor;

if rango == 0
    margen = max(1,abs(minValor))*margenRelativo;
else
    margen = margenRelativo*rango;
end

lim = 10.^[minValor - margen, maxValor + margen];
end

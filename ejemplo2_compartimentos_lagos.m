%% EJEMPLO 2: GRANDES LAGOS
% ejemplo2_compartimentos_lagos - Great Lakes linear compartment model
%
% Applies Adams-Bashforth, Adams-Moulton, implicit Taylor and implicit
% Runge-Kutta to the linear compartmental system
%
%       c'(t) = A*c(t),
%
% where c_i(t) is the relative pollutant concentration in lake i
% (Superior, Michigan, Huron, Erie, Ontario), and compares each method
% against the exact solution c(t) = expm(A*(t-t0))*c0.
%
% Outputs (written to figDir = <script folder>/figuras_grandes_lagos):
%   - PNG figures: per-lake concentrations, total relative pollution,
%     and global error vs. the reference solution, per method and combined
%   - A LaTeX .tex fragment with \includegraphics blocks for the TFG
%
% No inputs; run as a script. Requires ABashforth1, AMoulton1,
% TaylorImplicito1, coefRKImplicito, RungeKuttaImplicito1 on the path.

% Modelo compartimental lineal:
%
%       c'(t) = A c(t)
%
% donde c_i(t) representa la concentración relativa de contaminante
% en cada lago.

clear; close all; clc;

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

% Carpeta donde se guardaran las figuras.
figDir = fullfile(baseDir,'figuras_grandes_lagos');
if ~isfolder(figDir)
    mkdir(figDir);
end

%% Parámetros del modelo

% Volúmenes de los lagos
VS = 2900;   % Superior
VM = 1180;   % Michigan
VH = 850;    % Huron
VE = 116;    % Erie
VO = 393;    % Ontario

% Flujos entre lagos
qS = 15;     % Superior -> Huron
qM = 38;     % Michigan -> Huron
qH = 68;     % Huron -> Erie
qE = 85;     % Erie -> Ontario
qO = 99;     % Ontario -> salida

V = [VS; VM; VH; VE; VO];

% Matriz del sistema c' = A c
A = [ -qS/VS,      0,       0,       0,       0;
          0,  -qM/VM,       0,       0,       0;
      qS/VH,   qM/VH,  -qH/VH,       0,       0;
          0,       0,   qH/VE,  -qE/VE,       0;
          0,       0,       0,   qE/VO,  -qO/VO ];

f = @(t,c) A*c;

%% Datos de mallado

t0 = 0;
T  = 100;
n  = 1000;
ta=[t0,T];

t = linspace(t0,T,n+1);
h = (T-t0)/n;

% Concentración inicial relativa
c0 = ones(5,1);

%% Solución de referencia

C_ref = zeros(5,n+1);

for i = 1:n+1
    C_ref(:,i) = expm(A*(t(i)-t0))*c0;
end

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
[tAB,cAB,qAB,consAB] = ABashforth1(f,ta,c0,wAB,n);

fprintf('Ejecutando Adams-Moulton...\n');
[tAM,cAM,qAM,consAM] = AMoulton1(f,ta,c0,wAM,n,muAM);

fprintf('Ejecutando Taylor implicito...\n');
[tTI,cTI,qTI,consTI,infoTI] = TaylorImplicito1(f,ta,c0,gradoTaylor,n,tol,maxit);

fprintf('Construyendo coeficientes RK implicito de Gauss...\n');
[B,a,cRKcoef,infoCoefRK] = coefRKImplicito(sRK,'gauss');

fprintf('Ejecutando Runge-Kutta implicito...\n');
[tRK,cRK,sRKsalida,consRK,infoRK] = RungeKuttaImplicito1(f,ta,c0,B,a,n,tol,maxit);

%% 5. Solucion de referencia y estructura comun

% Mallado de referencia. Se toma el de Adams-Bashforth, que coincide con
% el resto porque todos los metodos usan el mismo intervalo y el mismo n.
tRef = tAB(:).';
cRef = solucionReferenciaLineal(A,c0,tRef,t0);

metodos = struct([]);

metodos(1).nombre = 'Adams-Bashforth';
metodos(1).nombreCorto = 'ABashforth';
metodos(1).t = tAB(:).';
metodos(1).c = orientarSolucion(cAB,c0,tAB);
metodos(1).orden = qAB;
metodos(1).consistente = consAB;

metodos(2).nombre = 'Adams-Moulton';
metodos(2).nombreCorto = 'AMoulton';
metodos(2).t = tAM(:).';
metodos(2).c = orientarSolucion(cAM,c0,tAM);
metodos(2).orden = qAM;
metodos(2).consistente = consAM;

metodos(3).nombre = 'Taylor implicito';
metodos(3).nombreCorto = 'TaylorImplicito';
metodos(3).t = tTI(:).';
metodos(3).c = orientarSolucion(cTI,c0,tTI);
metodos(3).orden = qTI;
metodos(3).consistente = consTI;

metodos(4).nombre = 'Runge-Kutta implicito';
metodos(4).nombreCorto = 'RKImplicito';
metodos(4).t = tRK(:).';
metodos(4).c = orientarSolucion(cRK,c0,tRK);
metodos(4).orden = infoCoefRK.ordenTeorico;
metodos(4).consistente = consRK;

%% 6. Graficas individuales: concentraciones y contaminacion total

M0 = V.'*c0;

for j = 1:numel(metodos)

    t = metodos(j).t;
    C = metodos(j).c;
    Cex = solucionReferenciaLineal(A,c0,t,t0);

    Mnum = (V.'*C)/M0;
    Mex  = (V.'*Cex)/M0;

    % Concentraciones frente al tiempo.
    fig = figure('Visible','off');
    tiledlayout(5,1,'TileSpacing','compact','Padding','compact');

    nombresLagos = {'Superior','Michigan','Huron','Erie','Ontario'};

    for r = 1:5
        nexttile;
        plot(t,Cex(r,:),'LineWidth',1.2); hold on;
        plot(t,C(r,:),'--','LineWidth',1.0);
        grid on;
        ylabel(['$c_',num2str(r),'(t)$'],'Interpreter','latex');
        title(nombresLagos{r},'Interpreter','latex');

        if r == 1
            legend('Referencia','Numerica','Location','best');
        end

        if r < 5
            set(gca,'XTickLabel',[]);
        else
            xlabel('$t$','Interpreter','latex');
        end
    end

    sgtitle([metodos(j).nombre, ': concentraciones'],'Interpreter','latex');

    nombreFig = fullfile(figDir, ['grandes_lagos_',metodos(j).nombreCorto,'_concentraciones.png']);
    exportgraphics(fig,nombreFig,'Resolution',300);
    close(fig);

    % Contaminacion total relativa.
    fig = figure('Visible','off');
    plot(t,Mex,'LineWidth',1.3); hold on;
    plot(t,Mnum,'--','LineWidth',1.1);
    grid on;
    xlabel('$t$','Interpreter','latex');
    ylabel('Contaminacion total relativa','Interpreter','latex');
    title([metodos(j).nombre, ': contaminacion total relativa'],'Interpreter','latex');
    legend('Referencia','Numerica','Location','best');

    nombreFig = fullfile(figDir, ['grandes_lagos_',metodos(j).nombreCorto,'_total.png']);
    exportgraphics(fig,nombreFig,'Resolution',300);
    close(fig);
end

%% 7. Graficas comparativas globales

% 7.1 Mosaico de concentraciones para todos los metodos.
fig = figure('Visible','off');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:numel(metodos)
    t = metodos(j).t;
    C = metodos(j).c;

    nexttile;
    plot(t,C.','LineWidth',1.0);
    grid on;
    xlabel('$t$','Interpreter','latex');
    ylabel('$c_i(t)$','Interpreter','latex');
    title(metodos(j).nombre,'Interpreter','latex');
    ylim([0,1.05]);
end

legend({'Superior','Michigan','Huron','Erie','Ontario'},'Location','best');
exportgraphics(fig,fullfile(figDir,'grandes_lagos_mosaico_concentraciones.png'),'Resolution',300);
close(fig);

% 7.2 Contaminacion total relativa.
fig = figure('Visible','off');
hold on; grid on;

Mref = (V.'*cRef)/M0;
plot(tRef,Mref,'LineWidth',1.3);

for j = 1:numel(metodos)
    t = metodos(j).t;
    C = metodos(j).c;
    Mnum = (V.'*C)/M0;
    plot(t,Mnum,'--','LineWidth',1.0);
end

xlabel('$t$','Interpreter','latex');
ylabel('Contaminacion total relativa','Interpreter','latex');
title('Contaminacion total relativa en el sistema','Interpreter','latex');
legend(['Referencia',{metodos.nombre}],'Location','best');
exportgraphics(fig,fullfile(figDir,'grandes_lagos_contaminacion_total.png'),'Resolution',300);
close(fig);

% 7.3 Error infinito respecto de la solucion de referencia.
fig = figure('Visible','off');
hold on; grid on;

for j = 1:numel(metodos)
    t = metodos(j).t;
    C = metodos(j).c;
    Cex = solucionReferenciaLineal(A,c0,t,t0);
    err = vecnorm(C-Cex,inf,1);
    semilogy(t,max(err,eps),'LineWidth',1.1);
end

xlabel('$t$','Interpreter','latex');
ylabel('$\|c_i-c(t_i)\|_{\infty}$','Interpreter','latex');
title('Error frente a la solucion de referencia','Interpreter','latex');
legend({metodos.nombre},'Location','best');
exportgraphics(fig,fullfile(figDir,'grandes_lagos_error_global.png'),'Resolution',300);
close(fig);

%% 8. Tabla resumen por pantalla

fprintf('\nResumen del experimento:\n');
fprintf('Metodo\t\t\tOrden/parametro\tConsistente\tError final\tError total final\n');

for j = 1:numel(metodos)
    t = metodos(j).t;
    C = metodos(j).c;
    Cex = solucionReferenciaLineal(A,c0,t,t0);

    errFinal = norm(C(:,end)-Cex(:,end),inf);

    Mnum = (V.'*C)/M0;
    Mex  = (V.'*Cex)/M0;
    errTotalFinal = abs(Mnum(end)-Mex(end));

    fprintf('%-22s\t%-8g\t\t%d\t\t%.3e\t%.3e\n', ...
        metodos(j).nombre, metodos(j).orden, metodos(j).consistente, errFinal, errTotalFinal);
end

%% 9. Generar fragmento LaTeX para insertar las figuras en el TFG

texFile = fullfile(figDir,'fragmento_grandes_lagos.tex');
fid = fopen(texFile,'w');

fprintf(fid,'%% Fragmento generado automaticamente por ejemplo2_compartimentos_lagos.m\n');
fprintf(fid,'%% Requiere \\usepackage{float}. Si usas [H].\n\n');

fprintf(fid,'\\begin{figure}[H]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.90\\textwidth]{figuras_grandes_lagos/grandes_lagos_mosaico_concentraciones.png}\n');
fprintf(fid,'    \\caption{Evolucion de las concentraciones relativas de contaminante en el modelo de los Grandes Lagos.}\n');
fprintf(fid,'    \\label{fig:grandes_lagos_concentraciones}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\begin{figure}[H]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.80\\textwidth]{figuras_grandes_lagos/grandes_lagos_contaminacion_total.png}\n');
fprintf(fid,'    \\caption{Evolucion de la contaminacion total relativa en el sistema de los Grandes Lagos.}\n');
fprintf(fid,'    \\label{fig:grandes_lagos_total}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\begin{figure}[H]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\includegraphics[width=0.80\\textwidth]{figuras_grandes_lagos/grandes_lagos_error_global.png}\n');
fprintf(fid,'    \\caption{Error global de los metodos numericos frente a la solucion de referencia del modelo de los Grandes Lagos.}\n');
fprintf(fid,'    \\label{fig:grandes_lagos_error_global}\n');
fprintf(fid,'\\end{figure}\n');

fclose(fid);

fprintf('\nFiguras guardadas en:\n%s\n', figDir);
fprintf('Fragmento LaTeX guardado en:\n%s\n', texFile);

%% FUNCIONES AUXILIARES DE GRAFICAS Y FORMATO

function C = solucionReferenciaLineal(A,c0,t,t0)
% solucionReferenciaLineal - Exact solution of the linear system c'=Ac
%
% Inputs:
%   A  - d x d system matrix
%   c0 - initial condition (column vector, length d)
%   t  - row/column vector of times at which to evaluate the solution
%   t0 - initial time
%
% Output:
%   C - d x length(t) matrix, C(:,i) = expm(A*(t(i)-t0))*c0

% Calcula la solucion exacta del sistema lineal c'=Ac en los puntos t.

    t = t(:).';
    m = length(c0);
    C = zeros(m,length(t));

    for i = 1:length(t)
        C(:,i) = expm(A*(t(i)-t0))*c0;
    end
end

function Y = orientarSolucion(Y,y0,t)
% orientarSolucion - Ensure a solution matrix is oriented variables x time
%
% Inputs:
%   Y  - solution matrix, either d x Nt or Nt x d
%   y0 - initial condition, used only to get the system dimension d
%   t  - time vector, used only to get Nt
%
% Output:
%   Y - the same matrix, transposed to d x Nt if needed. Errors if
%       neither orientation matches size(Y).

% Asegura el formato variables x tiempos.

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

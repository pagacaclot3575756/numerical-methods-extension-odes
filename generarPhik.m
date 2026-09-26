function [phi,dphi,psi,dpsi,phiSym,dphiSym,ySym,tSym]=generarPhik(f,d,grado)
% generarPhik - Symbolic Taylor-derivative and Taylor-method-slope generator
%
% Inputs:
%   f     - function handle, f(t,Y), defining the ODE system
%   d     - system dimension (number of equations)
%   grado - order q up to which the total derivatives are generated
%
% Output:
%   phi     - cell array of numeric function handles; phi{k}(t,Y) is the
%             (k-1)-th total derivative of y'=f(t,y), with phi{1}=f.
%             Generated for k=1,...,q+1 (the extra term is needed to
%             check the order of a degree-q method)
%   dphi    - cell array of numeric function handles; dphi{k}(t,Y) is
%             the Jacobian of phi{k} with respect to Y (dphi{1} is J_y f)
%   psi     - function handle psi(t,Z,h); implicit-Taylor-method slope,
%             psi_q(Z,t,h) = sum_{k=1}^q (-1)^(k+1) h^(k-1)/k! * phi_k(Z,t+h)
%   dpsi    - function handle dpsi(t,Z,h); Jacobian of psi w.r.t. Z
%   phiSym  - cell array with the symbolic expressions of phi{k}
%   dphiSym - cell array with the symbolic expressions of dphi{k}
%   ySym    - symbolic column vector of state variables (length d)
%   tSym    - symbolic time variable
%
% Uses the recurrence phi_1 = f, phi_{k+1} = J_y(phi_k)*f + d(phi_k)/dt
% to build the total-time-derivatives of the ODE needed for a Taylor
% predictor/method, then converts the symbolic expressions to numeric
% function handles via matlabFunction. psi/dpsi are built from phi/dphi
% via the helpers calcularPsiTaylor and calcularJyPsiTaylor.

%d=dimensión del sistema o núm de ecuaciones
%dphi=las matrices jacobianas respecto de y de phi_k
%dpsi=matriz jacobiana de psi_q respecto de y
%phiSym,dphiSym=versiones simbólicas de las funciones

%Seguimos la recurrencia:
% phi_1(y,t)=f(y,t)
% phi_{k+1}(y,t)=J_y(phi_k)(y,t)*f(y,t)+partial_t(phi_k)(y,t)
% Convencion MATLAB usada: f(t,Y)

%Y también calculamos:
% psi_q(Z,Y,t,h)=sum_{k=1}^q (-1)^(k+1) h^(k-1)/k! phi_k(Z,t+h)
%
% Para comprobar el orden de un metodo de grado q se necesita tambien
% phi_{q+1}. Por eso aqui se generan phi_1,...,phi_{q+1}, pero psi_q
% y su jacobiana se construyen solo con los terminos k=1,...,q.
q = grado;
numPhi = q + 1;
% Variables simbolicas
tSym = sym('t','real');
ySym = sym('Y',[d,1],'real'); %vector columna simbólico

% Evaluamos la function handle en variables simbolicas
fSym = f(tSym,ySym);
fSym = fSym(:); %Para que sea columna

%Inicializamos las celdas de los resultados
phiSym  = cell(1,numPhi);
dphiSym = cell(1,numPhi);

phi  = cell(1,numPhi);
dphi = cell(1,numPhi);

% phi_1(Y,t)=f(Y,t)
phiSym{1} = simplify(fSym);

% Construccion recursiva de phi_2,...,phi_grado
for k=1:numPhi-1

    % Jacobiana de phi_k respecto de Y
    Jy_phi = jacobian(phiSym{k},ySym);

    % Derivada parcial de phi_k respecto de t
    dt_phi = diff(phiSym{k},tSym);

    % phi_{k+1}=J_y(phi_k)*f + partial_t(phi_k)
    phiSym{k+1} = simplify(Jy_phi*fSym + dt_phi);

end

% Jacobianas respecto de Y de las phi_k y conversion a funciones numericas
for k=1:numPhi

    dphiSym{k} = simplify(jacobian(phiSym{k},ySym));
    %Las pasamos de funciones simbólicas a numéricas
    phi{k}  = matlabFunction(phiSym{k}, 'Vars', {tSym,ySym});
    dphi{k} = matlabFunction(dphiSym{k},'Vars', {tSym,ySym});

end

% Funcion pendiente del metodo Taylor implicito (Y no se usa)
psi = @(t,Z,h) calcularPsiTaylor(t,Z,h,phi,q);

% Jacobiana de psi respecto de la variable implicita Z
dpsi = @(t,Z,h) calcularJyPsiTaylor(t,Z,h,dphi,q);
end
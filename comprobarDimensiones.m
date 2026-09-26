function comprobarDimensiones(B,a,s)
% comprobarDimensiones - Check that B and a match the number of stages s
%
% Inputs:
%   B - candidate coefficient matrix
%   a - candidate weight vector
%   s - expected number of stages
%
% No output. Errors if size(B) ~= [s,s] or numel(a) ~= s.

% Comprueba que B y a tengan dimensiones compatibles con s.
if ~isequal(size(B),[s,s])
    error('La matriz B debe tener dimension s x s.');
end

if numel(a) ~= s
    error('El vector a debe tener s componentes.');
end
end
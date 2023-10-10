clearvars
close all

f = @(x) 8*x + 4;

nodes = (0:0.5:1.5)';
elem = [1, 2; 2, 3; 3,4];

numNodes = size(nodes,1);
numElem = size(elem,1);

kc = 0.57;    % for item (a)
T1 = 20.0;    % for item (b)
T4 = 25.3;    % for item (b)
Tinf = 22.0;  % for item (c)

M = [1 -1 0 0;
    -1 2 -1 0;
    0 -1 2 -1;
    0 0 -1 1];

h = 0.5;

% Item (a) 
K = kc*M/h;
sum(sum(K));
K(2,3);

clc
fprintf('Problema 1\n')
fprintf('(a) Sum of all coefficients of the stiffness matrix K: %f\n', sum(K(:)))
fprintf('    Hint 1. K(2,3) = %f\n', K(2,3))

% Item (b) 
kc = 0.5;
K = kc*M/h;
F = zeros(numNodes,1);
Q = zeros(numNodes,1);

for e = 1: numElem
    rows = [elem(e,1), elem(e,2)];
    x = (nodes(e) + nodes(e+1))/2;
    Fe = f(x)/4;
    F(rows) = F(rows) + Fe;
end

fixedNods = [1,4];
freeNods = setdiff(1:numNodes, fixedNods);

% B.C.
% Natural
Q(freeNods) = 0; % Redundant

% Essential
u = zeros(numNodes,1);
u(fixedNods) = [T1; T4];

Fm = F(freeNods) + Q(freeNods) - K(freeNods, fixedNods)* u(fixedNods);
Km = K(freeNods, freeNods);

um = Km\Fm;

u(freeNods) = um;

% Post-process
Q = K*u - F;

fprintf('(b) Q(4) = %.4f\n', Q(4))
fprintf('    Hint 2. Q(1) = %.4f\n', Q(1))

% Item (c)
A = [Km, [0;0]; 0, -1, T4-Tinf];
b = [Fm; F(4) - T4];
v = A\b;
fprintf('(c) beta = %.4f\n', v(3))

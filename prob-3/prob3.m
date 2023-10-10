clearvars;
close all;

f = 150.0;                   % (W/m^3) Heat generation at region A
indNodHint = 81;             % We give the temperature at node indNodHint

meshFile = 'mesh8x8Quad';    % (.m) file with the mesh's nodes and elements 
kc = 2.0;                    % Thermal conductivity
tempTop = 15.0;              % Temperature at the top boundary
tempBot = -10.0;             % Temperature at the bottom boundary
TinfLeft = 100.0;            % Bulk temperature at left edge
betaLeft = 30.0;             % Convection coefficient at left edge 
TinfRight = 55.0;            % Bulk temperature at right edge
betaRight = 45.0;            % Convection coefficient at right edge

%Point for interpolated temperature 
interpPointP = [0.31250, 0.31250];

%Point for interpolated temperature (1D interpolation & nodes on diagonal)
interpPointR = [0.81250,0.81250];

showPlt = 'Y';               % Show plots: yes

eval(meshFile);
numNod= size(nodes,1);
numElem= size(elem,1);

indNodTop = find(nodes(:,2) > 0.99);
indNodBot = find(nodes(:,2) < 0.01);
indNodLeft = find(nodes(:,1) < 0.01);
indNodRight = find(nodes(:,1) > 0.99);
indNodDiag = find( abs(nodes(:,1)-nodes(:,2)) < 0.01);

% Define Coefficients vector of the model equation
% In this case we use the Poisson coefficients defined in the problem above
a11 = kc;
a12=0.0;
a21=a12;
a22=a11;
a00=0.0;
coeff=[a11,a12,a21,a22,a00,0.0];

K=zeros(numNod);
F=zeros(numNod,1);
Q=zeros(numNod,1);

elemHeatGen = [];

for e=1:numElem
    %
    % x and y coordinates of the nodes of the element
    %
    x=nodes(elem(e,:),1); y=nodes(elem(e,:),2);

    if (min(y) < min(x) - 0.01) % Is element e (strictly) below the diagonal?
        elemHeatGen=[elemHeatGen,e];
        coeff(6) = f;
    else
        coeff(6) = 0.0;
    end  

    [Ke,Fe]=bilinearQuadElement(coeff,nodes,elem,e);

    %
    % Assemble the elements
    %
    rows=[elem(e,1); elem(e,2); elem(e,3); elem(e,4)];
    colums= rows;
    K(rows,colums)=K(rows,colums)+Ke; %assembly
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end %end for elements

%
% We save a copy of the initial F array for the postprocess step
Kini= K;
Fini= F;

%Booundary Conditions
fixedNodes= [indNodBot', indNodTop'];      %fixed Nodes (global numbering)
freeNodes= setdiff(1:numNod,fixedNodes);   %free Nodes (global numbering)

%------------- Convetion BC
indCV=indNodLeft';
[K,Q]=applyConvQuad(indCV,betaLeft,TinfLeft,K,Q,nodes,elem);
indCV=indNodRight';
[K,Q]=applyConvQuad(indCV,betaRight,TinfRight,K,Q,nodes,elem);

% ------------ Essential BC
u=zeros(numNod,1); %initialize uu vector
u(indNodTop)=tempTop;
u(indNodBot)=tempBot;
Fm=F(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);%here u can be 
                                                      %different from zero 
                                                      %only for fixed nodes
%Reduced system
Km=K(freeNodes,freeNodes);
Fm=Fm+Q(freeNodes);

%Compute the solution
%solve the System
format short e; %just to a better view of the numbers
um=Km\Fm;
u(freeNodes)=um;

%PostProcess: Compute secondary variables and plot results
Q=Kini*u-Fini;
titol='Equation solution';
colorScale='jet';
if showPlt == 'y' || showPlt == 'Y'
    plotContourSolution(nodes,elem,u,titol,colorScale);
end

clc
fprintf('%8s%5s%9s%9s%12s\n','Num.Nod','X','Y','T','Q')
fprintf('%5d%10.4f%9.4f%12.4e%12.4e\n',[(1:numNod)',nodes,u,Q]')

%Compute the maximum temperature at the nodes of region B
nodesHeatGen = elem(elemHeatGen, :);
nodesHeatGen = nodesHeatGen(:);
nodesHeatGen = unique(nodesHeatGen);
nodesWithoutHeatGen = setdiff(1:numNod, nodesHeatGen);
[maxTempB, idx] = max(u(nodesWithoutHeatGen));
nodMaxTempB = nodesWithoutHeatGen(idx); 

% Part (a)
fprintf(['\n',...
    '(a) Temp Max. in B is T= %.4eºC, at node %3d\n',...
    '    Hint. The temperature at node %d is, T= %.4e%cC\n'],...
    maxTempB,nodMaxTempB,indNodHint,u(indNodHint),char(176))

% Part (b)
%Compute the interpolated temperature at point p

for e=1:numElem
    nods = elem(e,:);
    vertexs= nodes(nods,:);
    [alphas,isInside] = baryCoordQuad(vertexs,interpPointP);
    if (isInside > 0)
        interpTempPointP = alphas * u(nods);
        fprintf(['\n',...
          '(b) Point P = (%.5f, %.5f) belongs to element number: %d\n',...
          '    Indexs of the nodes of elem %d: %d, %d, %d, %d\n', ...
          '    Interpolated temperature at point P, T= %.4e%cC\n'],...
          interpPointP,e,e,nods,interpTempPointP,char(176))
        break;
    end
end

xx = nodes(indNodDiag,1);
numNodsOnDiag = length(xx);
tempDiag = u(indNodDiag);
cf = polyfit(xx, tempDiag, numNodsOnDiag - 1);

interpTempPointR = polyval(cf, interpPointR(1,1));
averagedTempDiag = sum(tempDiag)/numNodsOnDiag;

% Part (c)
fprintf(['\n',...
    '(c) Interpolated temperature at point R = (%.5f, %.5f)\n', ...
    '    using nodes at the diagonal & 1D interpolation T = %.4e%cC\n', ...
    '    Hint. Averged temperature of the nodes at the diagonal, <T> = %.4e%cC\n'], ...
    interpPointR,interpTempPointR,char(176),averagedTempDiag,char(176))

%Plot the domain, the regions A and B, and the points P and R
if showPlt == 'y' || showPlt == 'Y'
    figure
    numbering = 0;
    %plotElements(nodes, elem, numbering);
    plotElementsOld(nodes, elem, numbering);
    hold on

    numElemHeatGen = length(elemHeatGen);
    X = nodes(elem(elemHeatGen,:)',1); X=reshape(X,4,numElemHeatGen);
    Y = nodes(elem(elemHeatGen,:)',2); Y=reshape(Y,4,numElemHeatGen);
    fill(X,Y,'yellow')

    plot(nodes(indNodTop,1),nodes(indNodTop,2), ...
        'ok','lineWidth',2,'markerFaceColor','red','markerSize',6)
    plot(nodes(indNodBot,1),nodes(indNodBot,2), ...
        'ok','lineWidth',2,'markerFaceColor','magenta','markerSize',6)
    plot(nodes(indNodLeft,1),nodes(indNodLeft,2), ...
        'ok','lineWidth',2,'markerFaceColor','blue','markerSize',6)
    plot(nodes(indNodRight,1),nodes(indNodRight,2), ...
        'ok','lineWidth',2,'markerFaceColor','blue','markerSize',6)     
    plot(nodes(indNodDiag,1),nodes(indNodDiag,2), ...
        'ok','lineWidth',2,'markerFaceColor','green','markerSize',6)

    plot(interpPointP(1,1),interpPointP(1,2), ...
     'ok','lineWidth',2,'markerFaceColor','red','markerSize',6)
    text(interpPointP(1,1)-0.02, interpPointP(1,2)+0.05, '$P$',...
     'interpreter','LaTeX','fontSize', 14)

    plot(interpPointR(1,1),interpPointR(1,2), ...
        'ok','lineWidth',2,'markerFaceColor','blue','markerSize',6)
    text(interpPointR(1,1)-0.02, interpPointR(1,2)+0.05, '$R$',...
        'interpreter','LaTeX','fontSize', 14)
    text(0.6, 0.35, '$A$', 'Interpreter', 'LaTeX', 'FontSize', 30, ...
    'color', 'black')
    hold off 
end


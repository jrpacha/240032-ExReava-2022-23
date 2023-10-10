function [stress,VonMisses]=computeQuadStrainStressVM(nodes, elem, u, C)
%
% For *quadrilateral elements* 
%
numNodes=size(nodes,1);


strain=zeros(numNodes,3);
stress=zeros(numNodes,3);
numElem=size(elem,1);
    gaussP1D=[-1/sqrt(3),1/sqrt(3)];
    order=[1, 4; 2, 3]; %align gauss points with element corners
    % Shape functions
    Psi1=@(x,y)(1-x).*(1-y)/4;
    Psi2=@(x,y)(1+x).*(1-y)/4;
    Psi3=@(x,y)(1+x).*(1+y)/4;
    Psi4=@(x,y)(1-x).*(1+y)/4;
    shapeFunctions = @(x,y)[Psi1(x,y),Psi2(x,y),Psi3(x,y),Psi4(x,y)];        
    % Shape function derivatives
    dPsi11=@(x,y) -(1-y)/4;
    dPsi21=@(x,y) (1-y)/4;
    dPsi31=@(x,y) (1+y)/4;
    dPsi41=@(x,y) -(1+y)/4;
    dPsi12=@(x,y) -(1-x)/4;
    dPsi22=@(x,y) -(1+x)/4;
    dPsi32=@(x,y) (1+x)/4;
    dPsi42=@(x,y) (1-x)/4;
    % Derivative matrix 2x4
    deriv =@(x,y) [dPsi11(x,y),dPsi21(x,y),dPsi31(x,y),dPsi41(x,y);...
                   dPsi12(x,y),dPsi22(x,y),dPsi32(x,y),dPsi42(x,y)];

stressNod=zeros(numNodes,3);
for e=1:numElem
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    v4=nodes(elem(e,4),:); 
    vertices=[v1;v2;v3;v4];
    rows=[elem(e,1)*2-1; elem(e,1)*2; elem(e,2)*2-1; elem(e,2)*2; ...
         elem(e,3)*2-1; elem(e,3)*2; elem(e,4)*2-1; elem(e,4)*2];
    ue=u(rows,:); %displacements for nodes in this element
    
    for i=1:2
        for j=1:2
            locDeriv=deriv(gaussP1D(i),gaussP1D(j));
            jacob=locDeriv*vertices;
            cartDeriv=jacob\locDeriv;
            B=[];
            for k=1:4
                B=[B,[cartDeriv(1,k), 0; 0,cartDeriv(2,k); cartDeriv(2,k),cartDeriv(1,k)]];
            end
            stressGauss(:,order(i,j)) = C*B*ue;
            nodParam = [1/gaussP1D(i), 1/gaussP1D(i)];
            extrapolAtVertex(order(i,j),:) = shapeFunctions(nodParam(1), nodParam(2));
        end
    end
    
    stressElem=extrapolAtVertex*stressGauss';
    stressNod(elem(e,:),:)=stressNod(elem(e,:),:)+stressElem;
end
nodInElem=zeros(numNodes,1);
for i=1:numNodes
    elements = elemContainNod(i, elem);
    nodInElem(i) = size(elements,2); %how many elements the node ith belongs
end
stress=stressNod./nodInElem;
for i=1:numNodes
    sTx=stress(i,1);
    sTy=stress(i,2);
    sTxy=stress(i,3);
    VonMisses(i)=sqrt(sTx^2+sTy^2-sTx*sTy+3*sTxy^2);
end
  

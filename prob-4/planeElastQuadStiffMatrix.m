function Ke=planeElastQuadStiffMatrix(nodes,elem,e,C,h)
    N=2; %number of GaussPoints
    [w,ptGaussRef]=gaussValues2DQuad(N);
    %
    % First compute Ke, Fe for each element
    %    
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    v4=nodes(elem(e,4),:); 
    vertices=[v1;v2;v3;v4];
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
    Jacob =@(x,y) [dPsi11(x,y), dPsi21(x,y),dPsi31(x,y),dPsi41(x,y);...
                   dPsi12(x,y), dPsi22(x,y),dPsi32(x,y),dPsi42(x,y)];
    % Compute the corresponding Gaussian points on the domain
    % evaluate Shape functions on Gaussian reference points
    xx = ptGaussRef(:,1);
    yy = ptGaussRef(:,2);
    % evaluate Jacobian contribution for each point 
    % We use a Matlab cell variable in order to load each matrix 
    numPtG=size(xx,1);  
    for i=1:numPtG
        Jtilde{i}=inv(Jacob(xx(i),yy(i))*vertices);
        evalDetJacob(i) = det(Jacob(xx(i),yy(i))*vertices);
        Jacobia{i}=Jacob(xx(i),yy(i)); %derivatives of the shape functions
        shapeFun{i}=shapeFunctions(xx(i),yy(i));
    end
    %
    % element stiff matrix
    %
        ik1=1; ik2=1; %index of the derivatives here x-x
        Kxx=compKijIntegQuad(ik1,ik2,Jacobia,w,evalDetJacob,Jtilde,numPtG);
        ik1=1; ik2=2; 
        Kxy=compKijIntegQuad(ik1,ik2,Jacobia,w,evalDetJacob,Jtilde,numPtG);
        ik1=2; ik2=1; 
        Kyx=compKijIntegQuad(ik1,ik2,Jacobia,w,evalDetJacob,Jtilde,numPtG);
        ik1=2; ik2=2; 
        Kyy=compKijIntegQuad(ik1,ik2,Jacobia,w,evalDetJacob,Jtilde,numPtG);
     KK={};   
     for i=1:4
         for j=1:4
           KK{i,j}=[C(1,1)*Kxx(i,j)+C(3,3)*Kyy(i,j), C(1,2)*Kxy(i,j)+C(3,3)*Kyx(i,j);
                    C(3,3)*Kxy(i,j)+C(2,1)*Kyx(i,j), C(3,3)*Kxx(i,j)+C(2,2)*Kyy(i,j)]; 
         end
     end
     
       Ke =h*[KK{1,1}, KK{1,2}, KK{1,3}, KK{1,4};
          KK{2,1}, KK{2,2}, KK{2,3}, KK{2,4};
          KK{3,1}, KK{3,2}, KK{3,3}, KK{3,4};
          KK{4,1}, KK{4,2}, KK{4,3}, KK{4,4}];
end
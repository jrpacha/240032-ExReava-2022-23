function Kij=compKijIntegQuad(ik1,ik2,Jacobia,w,evalDetJacob,Jtilde,numPtG)
Kij=zeros(4);
for i=1:4
    for j=1:4
        suma=0;
        for ptG=1:numPtG
            first=Jacobia{ptG}(1,i)*Jtilde{ptG}(ik1,1)+Jacobia{ptG}(2,i)*Jtilde{ptG}(ik1,2);
            second=Jacobia{ptG}(1,j)*Jtilde{ptG}(ik2,1)+Jacobia{ptG}(2,j)*Jtilde{ptG}(ik2,2);
            suma=suma+w(ptG)*first*second*evalDetJacob(ptG);
        end
        Kij(i,j)=suma;
    end
end
end
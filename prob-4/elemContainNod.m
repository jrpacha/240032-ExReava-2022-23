function elements = elemContainNod(indNod, elem)
nNodEl=size(elem,2);
i=indNod;
if (nNodEl==3)
    nElem1 = find(elem(:,1)==i);
    nElem2 = find(elem(:,2)==i);
    nElem3 = find(elem(:,3)==i);
    elements=[nElem1',nElem2',nElem3']
elseif (nNodEl==4)
    nElem1 = find(elem(:,1)==i);
    nElem2 = find(elem(:,2)==i);
    nElem3 = find(elem(:,3)==i);
    nElem4 = find(elem(:,4)==i);
    elements=[nElem1',nElem2',nElem3',nElem4'];
end   

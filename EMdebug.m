Y = [1 0 1 1 1;1 1 0 0 1;1 1 1 0 0]
readA= Y(:,1)/sum(Y(:,1));
readB= Y(:,2)/sum(Y(:,2));
readC= Y(:,3)/sum(Y(:,3));
readD= Y(:,4)/sum(Y(:,4));
readE= Y(:,5)/sum(Y(:,5));
INTreadA= Y(:,1)/sum(Y(:,1));
INTreadB= Y(:,2)/sum(Y(:,2));
INTreadC= Y(:,3)/sum(Y(:,3));
INTreadD= Y(:,4)/sum(Y(:,4));
INTreadE= Y(:,5)/sum(Y(:,5));
count = 0;
dParameter1 = 1;
%% Section (1A):
while dParameter1 > 0.01
    iterA = readA;
    iterB = readB;
    iterC = readC;
    iterD = readD;
    iterE = readE;
    %Single iteration prep (M step)
    totalChance = (readA(1,:)+readB(1,:)+readC(1,:)+readD(1,:)+readE(1,:))+(readA(2,:)+readB(2,:)+readC(2,:)+readD(2,:)+readE(2,:))+(readA(3,:)+readB(3,:)+readC(3,:)+readD(3,:)+readE(3,:));
    pred = (readA(1,:)+readB(1,:)+readC(1,:)+readD(1,:)+readE(1,:))/totalChance;
    pgreen = (readA(2,:)+readB(2,:)+readC(2,:)+readD(2,:)+readE(2,:))/totalChance;
    pblue = (readA(3,:)+readB(3,:)+readC(3,:)+readD(3,:)+readE(3,:))/totalChance;
    %Set new read percentag
    readA(1,:) = INTreadA(1,:)*pred;
    readB(1,:) = INTreadB(1,:)*pred;
    readC(1,:) = INTreadC(1,:)*pred;
    readD(1,:) = INTreadD(1,:)*pred;
    readE(1,:) = INTreadE(1,:)*pred;
    readA(2,:) = INTreadA(2,:)*pgreen;
    readB(2,:) = INTreadB(2,:)*pgreen;
    readC(2,:) = INTreadC(2,:)*pgreen;
    readD(2,:) = INTreadD(2,:)*pgreen;
    readE(2,:) = INTreadE(2,:)*pgreen;
    readA(3,:) = INTreadA(3,:)*pblue;
    readB(3,:) = INTreadB(3,:)*pblue;
    readC(3,:) = INTreadC(3,:)*pblue;
    readD(3,:) = INTreadD(3,:)*pblue;
    readE(3,:) = INTreadE(3,:)*pblue;
    %Normalize read percentages
    readA = readA/sum(readA);
    readB = readB/sum(readB);
    readC = readC/sum(readC);
    readD = readD/sum(readD);
    readE = readE/sum(readE);
    %change in parameter
    dParaMatrix = zeros(3,5);
    dParaMatrix(1,1) = abs(readA(1,:)-iterA(1,:));
    dParaMatrix(2,1) = abs(readA(2,:)-iterA(2,:));
    dParaMatrix(3,1) = abs(readA(3,:)-iterA(3,:));
    dParaMatrix(1,2) = abs(readB(1,:)-iterB(1,:));
    dParaMatrix(2,2) = abs(readB(2,:)-iterB(2,:));
    dParaMatrix(3,2) = abs(readB(3,:)-iterB(3,:));
    dParaMatrix(1,3) = abs(readC(1,:)-iterC(1,:));
    dParaMatrix(2,3) = abs(readC(2,:)-iterC(2,:));
    dParaMatrix(3,3) = abs(readC(3,:)-iterC(3,:));
    dParaMatrix(1,4) = abs(readD(1,:)-iterD(1,:));
    dParaMatrix(2,4) = abs(readD(2,:)-iterD(2,:));
    dParaMatrix(3,4) = abs(readD(3,:)-iterD(3,:));
    dParaMatrix(1,5) = abs(readE(1,:)-iterE(1,:));
    dParaMatrix(2,5) = abs(readE(2,:)-iterE(2,:));
    dParaMatrix(3,5) = abs(readE(3,:)-iterE(3,:));
    dParameterRow = max(dParaMatrix);
    dParameter1 = max(dParameterRow);
    count = count + 1;
end
readA
readB
readC
readD
readE
pred
pgreen
pblue
count
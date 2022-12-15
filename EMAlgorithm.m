%% Problem 1 - RNA-seq
% Gabriel Ketron 916571733 Bim 154 Hw3
% (Sorry for length, much of the code is copied over instead of functions.)
%Initialization (initial E step)
Y = [1 0 1 1 1;1 1 0 0 1;1 1 1 0 0];
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
dParameter1;
pred;
pgreen;
pblue;
count;
% This algorithm requires 6 iteratons to convergence. This is the same
% solution from class. 
%% Section (1B)
Y = [1 0 1 1 1;1 1 0 0 1;1 1 1 0 0];
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
while dParameter1 > 0.001
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
dParameter1;
pred;
pgreen;
pblue;
count;
% 9 iterations are required.
%% Section (1C)
%New Initialization 1: 
Y = [0 0 1 1 1;1 1 0 0 1;1 1 1 0 0];
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
dParameter1;
pred;
pgreen;
pblue;
X = [pred pgreen pblue];
labels = {'pred','pgreen','pblue'};
figure
pie(X,labels)
snapnow
count;
Y = [0 0 1 1 1;1 1 0 0 1;1 1 1 0 0];
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
while dParameter1 > 0.001
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
    if count < 5
        X = [pred pgreen pblue];
        labels = {'pred','pgreen','pblue'};
        figure
        pie(X,labels)
    end
end
dParameter1;
pred;
pgreen;
pblue;
X = [pred pgreen pblue];
labels = {'pred','pgreen','pblue'};
figure
pie(X,labels)
snapnow
count;
%This set converges for both 0.01 and 0.001 cutoffs.
Y = [0 0 1 1 1;0 1 0 0 1;1 1 1 0 0];
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
dParameter1;
pred;
pgreen;
pblue;
X = [pred pgreen pblue];
labels = {'pred','pgreen','pblue'};
figure
pie(X,labels)
snapnow
count;
Y = [0 0 1 1 1;0 1 0 0 1;1 1 1 0 0];
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
while dParameter1 > 0.001
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
    if count < 5
        X = [pred pgreen pblue];
        labels = {'pred','pgreen','pblue'};
        figure
        pie(X,labels)
        snapnow
    end
end
dParameter1;
pred;
pgreen;
pblue;
X = [pred pgreen pblue];
labels = {'pred','pgreen','pblue'};
figure
pie(X,labels)
snapnow
count;
%This initialization converges to be without any green. 
%% Section (1D)
% New likelihood function: L(P) = (pred/sum(p)+pgreen/sum(p)+pblue/sum(p))^2*(pgreen/sum(p)+pblue/sum(p))*((pred/sum(p)+pblue/sum(p))*(pred/sum(p))*(pred/sum(p)+pgreen/sum(p))*(pgreen/sum(p))
% log L(P) = log(pgreen + pblue) + log(pred + pblue)+ log(pred) + log(pred+pgreen)+log(pgreen).
% The new likelihood function is still concave since it converges to a
% local maximum. 
%% Section(1E)
Y = [1 0 1 1 1 1 0;1 1 0 0 1 1 1;1 1 1 0 0 1 0];
readA= Y(:,1)/sum(Y(:,1));
readB= Y(:,2)/sum(Y(:,2));
readC= Y(:,3)/sum(Y(:,3));
readD= Y(:,4)/sum(Y(:,4));
readE= Y(:,5)/sum(Y(:,5));
readF= Y(:,6)/sum(Y(:,6));
readG= Y(:,7)/sum(Y(:,7));
INTreadA= Y(:,1)/sum(Y(:,1));
INTreadB= Y(:,2)/sum(Y(:,2));
INTreadC= Y(:,3)/sum(Y(:,3));
INTreadD= Y(:,4)/sum(Y(:,4));
INTreadE= Y(:,5)/sum(Y(:,5));
INTreadF= Y(:,6)/sum(Y(:,6));
INTreadG= Y(:,7)/sum(Y(:,7));
count = 0;
dParameter1 = 1;
while dParameter1 > 0.001
    iterA = readA;
    iterB = readB;
    iterC = readC;
    iterD = readD;
    iterE = readE;
    iterF = readF;
    iterG = readG;
    %Single iteration prep (M step)
    totalChance = (readA(1,:)+readB(1,:)+readC(1,:)+readD(1,:)+readE(1,:)+readF(1,:)+readG(1,:))+(readA(2,:)+readB(2,:)+readC(2,:)+readD(2,:)+readE(2,:)+readF(2,:)+readG(2,:))+(readA(3,:)+readB(3,:)+readC(3,:)+readD(3,:)+readE(3,:)+readF(3,:)+readG(3,:));
    pred = (readA(1,:)+readB(1,:)+readC(1,:)+readD(1,:)+readE(1,:)+readF(1,:)+readG(1,:))/totalChance;
    pgreen = (readA(2,:)+readB(2,:)+readC(2,:)+readD(2,:)+readE(2,:)+readF(2,:)+readG(2,:))/totalChance;
    pblue = (readA(3,:)+readB(3,:)+readC(3,:)+readD(3,:)+readE(3,:)+readF(3,:)+readG(3,:))/totalChance;
    %Set new read percentag
    readA(1,:) = INTreadA(1,:)*pred;
    readB(1,:) = INTreadB(1,:)*pred;
    readC(1,:) = INTreadC(1,:)*pred;
    readD(1,:) = INTreadD(1,:)*pred;
    readE(1,:) = INTreadE(1,:)*pred;
    readF(1,:) = INTreadF(1,:)*pred;
    readG(1,:) = INTreadG(1,:)*pred;
    readA(2,:) = INTreadA(2,:)*pgreen;
    readB(2,:) = INTreadB(2,:)*pgreen;
    readC(2,:) = INTreadC(2,:)*pgreen;
    readD(2,:) = INTreadD(2,:)*pgreen;
    readE(2,:) = INTreadE(2,:)*pgreen;
    readF(2,:) = INTreadF(2,:)*pgreen;
    readG(2,:) = INTreadG(2,:)*pgreen;
    readA(3,:) = INTreadA(3,:)*pblue;
    readB(3,:) = INTreadB(3,:)*pblue;
    readC(3,:) = INTreadC(3,:)*pblue;
    readD(3,:) = INTreadD(3,:)*pblue;
    readE(3,:) = INTreadE(3,:)*pblue;
    readF(3,:) = INTreadF(3,:)*pblue;
    readG(3,:) = INTreadG(3,:)*pblue;
    %Normalize read percentages
    readA = readA/sum(readA);
    readB = readB/sum(readB);
    readC = readC/sum(readC);
    readD = readD/sum(readD);
    readE = readE/sum(readE);
    readF = readF/sum(readF);
    readG = readG/sum(readG);
    %change in parameter
    dParaMatrix = zeros(3,7);
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
    dParaMatrix(1,6) = abs(readF(1,:)-iterF(1,:));
    dParaMatrix(2,6) = abs(readF(2,:)-iterF(2,:));
    dParaMatrix(3,6) = abs(readF(3,:)-iterF(3,:));
    dParaMatrix(1,7) = abs(readG(1,:)-iterG(1,:));
    dParaMatrix(2,7) = abs(readG(2,:)-iterG(2,:));
    dParaMatrix(3,7) = abs(readG(3,:)-iterG(3,:));
    dParameterRow = max(dParaMatrix);
    dParameter1 = max(dParameterRow);
    count = count + 1;
    if count < 5
        X = [pred pgreen pblue];
        labels = {'pred','pgreen','pblue'};
        figure
        pie(X,labels)
        snapnow
    end
end
dParameter1;
readA;
count;
%system converges from multiple intializations to minimize pblue. starting
%with green and blue at 50/50 instead of 3/3rds causes convergence to 10%
%blue instead of 0.2%. red/blue 50/50 similarly converges blue to only 1%.
%% Section (1F)
% i.Likelihood: L(P) = (pred/sum(p)+pgreen/sum(p)+pblue/sum(p))^100*(pgreen/sum(p)+pblue/sum(p))^150*((pred/sum(p)+pblue/sum(p))^200*(pred/sum(p))^250*(pred/sum(p)+pgreen/sum(p))^300*(pgreen/sum(p))^400
%Log Likelihood:
%150log(pgreen/sum(p)+pblue/sum(p))+200log(pred/sum(p)+pblue/sum(p))+250log(pred/sum(p))+300log(pred/sum(p)+pgreen/sum(p))+400log(pgreen/sum(p)).
aReads = [1; 1; 1];
aReads = repmat(aReads,1,100);
bReads = [0; 1; 1];
bReads = repmat(bReads,1,150);
cReads = [1; 0; 1];
cReads = repmat(cReads,1,200);
dReads = [1; 0; 0];
dReads = repmat(dReads,1,250);
eReads = [1; 1; 0];
eReads = repmat(eReads,1,300);
gReads = [0; 1; 0];
gReads = repmat(gReads,1,400);

Y = [aReads bReads cReads dReads eReads gReads];
dY = Y;
count = 0;
dP = 1;
while dP > 0.001
    iterY = dY;
    totalchance = sum(dY,'all');
    pred = sum(dY(1,:),'all')./totalchance;
    pgreen = sum(dY(2,:),'all')./totalchance;
    pblue = sum(dY(3,:),'all')./totalchance;
    pChart = [pred pgreen pblue];
    for k = 1:3
        dY(k,:) = Y(k,:).*pChart(k);
    end
    for i = 1:1400
            dY(:,i) = dY(:,i)./sum(dY(:,i),'all');
    end 
    dParaMatrix = zeros(3,1400);
    for k = 1:3
        for i = 1:1400
            dParaMatrix(k,i) = abs(dY(k,i)-iterY(k,i));
        end
    end
    dParameterRow = max(dParaMatrix);
    dP = max(dParameterRow);
    count = count + 1;
end
dY(:,1);
dP;
count;
pSave = pChart;
pSabund = pChart.*1400;
% ii. The results are similar as the single read sample. This likely comes about from a lack of blue only reads, so the likelihood function is drawn over times to more prevalent readings.
% iii. The number of reads does not increase operational complexity. While
% the structure of the code would give the appearance of going from O(n) to
% O(n^2), both were O(n^2) in order to make the parameter difference
% matrix. 
%% Section (1G)
% The derivation of E and M would change. The relative length effects the
% likelihood function by dividing the sampling probability. Since new
% locations exist for the read to attach to on the same transcript, these
% differences apply by substituting At for p(t)*l(t)/sum(p(r)l(r)) where t
% is the investigated transcript and r is a summation across all existing
% transcripts. This would effect the E step as the relative lengths would
% need to be multiplied into the top and bottom.

%% Problem 2 
DUMBHO = 170;
OOF = 7;

%% Section (2A)

r = normrnd(170,7,100,1);
heightAVG= mean(r,'all')
histfit(r)
title('Normal pdf and Histogram')
xlabel('Observation')
ylabel('Occurances')
legend('Histogram','Normal Distribution')

%% Section (2B)

y = randsample(r,100,true);

%i.Probability of selecting same dataset: (1/100)^100
%ii. Probabiltiy of finding the same dataset of n observations in x attempts is
%P(X=<x) = 1-(1-(1/n)^n)^x
%iii. P(X=<5) = (1-(1-(1/1000)^1000)^5 = 0 (effectively).

%Z% Section (2c) 

x = [];
for index = 1:1000
    omniset = datasample(r,100);
    omnisetavg = mean(omniset);
    x(end+1) = omnisetavg;
end
x;
xavg = mean(x);
figure
hist2= histogram(x)

%% Section (2D)
Xstd = std(x)
StandardErrorM = Xstd/sqrt(length(x));
ts = tinv([0.025 0.975],length(x)-1);
CI = xavg + ts*StandardErrorM
xprime = x';
pd = fitdist(xprime,'Normal');
ci = paramci(pd)
%CI values agree

%% Section (2E)
x = normrnd(170,7,1000,1);
x;
xavg = mean(x);
figure
hist3= histogram(x)
hold off
Xstd = std(x);
StandardErrorM = Xstd/sqrt(length(x));
ts = tinv([0.025 0.975],length(x)-1);
CI = xavg + ts*StandardErrorM
xprime = x';
pd = fitdist(x,'Normal');
ci = paramci(pd)

%% Section (2F)
r = normrnd(170,7,1000,1);
x = [];
for index = 1:50
    omniset = datasample(r,1000);
    omnisetavg = mean(omniset);
    x(end+1) = omnisetavg;
end
x;
xavg = mean(x);
figure
hist4= histogram(x)
Xstd = std(x);
StandardErrorM = Xstd/sqrt(length(x));
ts = tinv([0.025 0.975],length(x)-1);
CI = xavg + ts*StandardErrorM
xprime = x';
pd = fitdist(xprime,'Normal');
ci = paramci(pd)
%the CI for the 50 iteration version is much tighter then the 1000
%original samples/

%% Section (2g)
%i.
x = [];
for index = 1:1000
    r = randsample(10,1);
    if r == 1|2|3
        val = normrnd(0,1,100,1);
    end
    if r == 4|5
        val = normrnd(3,0.1,100,1);
    end
    if r > 5
        val = normrnd(10,1,100,1);
    end
    x = [x val];
end
xi = [];
for index = 1:1000
    r = randsample(100,1);
    r = x(:, r);
    omniset = datasample(r,100);
    omnisetavg = mean(omniset);
    xi(end+1) = omnisetavg;
end
xavg = mean(xi);
figure
hist5= histogram(xi)
title('Hist 5')
hold on
Xstd = std(xi);
StandardErrorM = Xstd/sqrt(length(xi));
ts = tinv([0.025 0.975],length(xi)-1);
CI = xavg + ts*StandardErrorM
xprime = xi';
pd = fitdist(xprime,'Normal');
ci = paramci(pd)
%ii.
xi = [];
for index = 1:1000
    r = randsample(100,1);
    r = x(:, r);
    r = r';
    omniset = datasample(r,100);
    omnisetavg = mean(omniset);
    xi(end+1) = omnisetavg;
end
xi;
xavg = mean(xi);
figure
hist6= histogram(xi)
title('hist 6')
hold off
%iii. Since the bootstrap is based of of the different percentage
%allocations of the mixng coefficients, it also reflects the random
%selection of distribution found in the original sample. 
%% Problem 3 Section (3A)

%Starting with the read matrix Y, we can generated simulated samples of
%1400 reads by randomly selecting with replacement 1400 new reads. These
%will form new datasets that could then be used to create a set of EM
%aquired abundances.

%% Section (3B)
pRuns = [];
for j = 1:10
    sampledY = Y;
    for i = 1:1400
        r = randsample(1400,1);
        sampledY(:,i)=Y(:,r);
    end
    dY = sampledY;
    count = 0;
    dP = 1;
    while dP > 0.001
        iterY = dY;
        totalchance = sum(dY,'all');
        pred = sum(dY(1,:),'all')./totalchance;
        pgreen = sum(dY(2,:),'all')./totalchance;
        pblue = sum(dY(3,:),'all')./totalchance;
        pChart = [pred; pgreen; pblue];
        for k = 1:3
            dY(k,:) = Y(k,:).*pChart(k);
        end
        for i = 1:1400
                dY(:,i) = dY(:,i)./sum(dY(:,i),'all');
        end 
        dParaMatrix = zeros(3,1400);
        for k = 1:3
            for i = 1:1400
                dParaMatrix(k,i) = abs(dY(k,i)-iterY(k,i));
            end
        end
        dParameterRow = max(dParaMatrix);
        dP = max(dParameterRow);
        count = count + 1;
    end
    pRuns = [pRuns pChart];
end
pRuns

%% Section (3C)
pSave
pRunR =(pRuns(1,:));
pRunG =(pRuns(2,:));
pRunB =(pRuns(3,:));

diffR = pRunR - pSave(1);
diffG = pRunG - pSave(2);
diffB = pRunB - pSave(3);

dSqR = diffR.^2;
dSqG = diffG.^2;
dSqB = diffB.^2;

sqSumR = sum(dSqR,'all')/9;
sqSumG = sum(dSqG,'all')/9;
sqSumB = sum(dSqB,'all')/9;

stdRed = sqrt(sqSumR);
stdGreen = sqrt(sqSumG);
stdBlue = sqrt(sqSumB);
stderror = [stdRed, stdGreen, stdBlue];
l = 1:3;
figure
bar(l,pSave)

hold on

er = errorbar(l,pSave,stderror);
er.Color = [0 0 0];
er.LineStyle = 'none';

hold off
%Error bars are small because there is almost no variation. 

%% Section (3D)
pAbund = [];
for j = 1:10
    sampledY = Y;
    for i = 1:1400
        r = randsample(1400,1);
        sampledY(:,i)=Y(:,r);
    end
    dY = sampledY;
    count = 0;
    dP = 1;
    while dP > 0.001
        iterY = dY;
        totalchance = sum(dY,'all');
        pred = sum(dY(1,:),'all')./totalchance;
        pgreen = sum(dY(2,:),'all')./totalchance;
        pblue = sum(dY(3,:),'all')./totalchance;
        pChart = [pred; pgreen; pblue];
        for k = 1:3
            dY(k,:) = Y(k,:).*pChart(k);
        end
        for i = 1:1400
                dY(:,i) = dY(:,i)./sum(dY(:,i),'all');
        end 
        dParaMatrix = zeros(3,1400);
        for k = 1:3
            for i = 1:1400
                dParaMatrix(k,i) = abs(dY(k,i)-iterY(k,i));
            end
        end
        dParameterRow = max(dParaMatrix);
        dP = max(dParameterRow);
        count = count + 1;
    end
    pCount = (pChart.*1400);
    pAbund = [pAbund pCount];
end
pAbund

pRunR =(pAbund(1,:));
pRunG =(pAbund(2,:));
pRunB =(pAbund(3,:));

diffR = pRunR - pSabund(1);
diffG = pRunG - pSabund(2);
diffB = pRunB - pSabund(3);

dSqR = diffR.^2;
dSqG = diffG.^2;
dSqB = diffB.^2;

sqSumR = sum(dSqR,'all')/9;
sqSumG = sum(dSqG,'all')/9;
sqSumB = sum(dSqB,'all')/9;

stdRed = sqrt(sqSumR);
stdGreen = sqrt(sqSumG);
stdBlue = sqrt(sqSumB);
stderror = [stdRed, stdGreen, stdBlue];

l = 1:3;
figure
bar(l,pSabund)

hold on

er = errorbar(l,pSave,stderror);
er.Color = [0 0 0];
er.LineStyle = 'none';

hold off
%Error bars are small because there is almost no variation. 

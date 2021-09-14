%% QPSK Without Noise for Checking Viterbi Algorithm

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Bits and Sending Them
numOfBits = 10^4;

%Creating bits:
bits = (2 * randi([0,1], [numOfBits,1])-1)/sqrt(2) + 1i * (2 * randi([0,1], [numOfBits,1])-1)/sqrt(2);

%Creating uk based on part a:
u = zeros(numOfBits, 1);

for i=1:numOfBits
    if i == 1
        u(i) = (5/6) * bits(i) + (0.5) * bits(i+1) + (1/12) * bits(i+2);
    elseif i == 2
        u(i) = (5/6) * bits(i) + (0.5) * bits(i+1) + (1/12) * bits(i+2) + (0.5) * bits(i-1);
    elseif i == numOfBits-1
        u(i) = (5/6) * bits(i) + (0.5) * bits(i+1) + (1/12) * bits(i-2) + (0.5) * bits(i-1);
    elseif i == numOfBits
        u(i) = (5/6) * bits(i) + (1/12) * bits(i-2) + (0.5) * bits(i-1);
    else
        u(i) = (5/6) * bits(i) + (0.5) * bits(i+1) + (1/12) * bits(i+2) + (0.5) * bits(i-1) + (1/12) * bits(i-2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Detecting:

%Viterbi Algorithm and Detecting Bits:
c = zeros(4, numOfBits);
I = [(1/sqrt(2))+1i*(1/sqrt(2)); (-1/sqrt(2))+1i*(1/sqrt(2)); (1/sqrt(2))-1i*(1/sqrt(2)); (-1/sqrt(2))-1i*(1/sqrt(2))];
foundBits = zeros(numOfBits, 1);

for i=1:12
    if i==1
        c(1,i) = real( I(1)' * ( 2*u(i)-(5/6)*I(1) ) );
        c(2,i) = real( I(2)' * ( 2*u(i)-(5/6)*I(2) ) );
        c(3,i) = real( I(3)' * ( 2*u(i)-(5/6)*I(3) ) );
        c(4,i) = real( I(4)' * ( 2*u(i)-(5/6)*I(4) ) );
        
        foundPath = [1; 2; 3; 4];
        
    elseif i==2
        
        cCalc1 = c([1;2;3;4], i-1) + real( I(1)' * ( 2*u(i)-(5/6)*I(1)-I(foundPath) ) );
        [c(1,i), maxIndex1] = max(cCalc1);
        
        cCalc2 = c([1;2;3;4], i-1) + real( I(2)' * ( 2*u(i)-(5/6)*I(2)-I(foundPath) ) );
        [c(2,i), maxIndex2] = max(cCalc2);
        
        cCalc3 = c([1;2;3;4], i-1) + real( I(3)' * ( 2*u(i)-(5/6)*I(3)-I(foundPath) ) );
        [c(3,i), maxIndex3] = max(cCalc3);
        
        cCalc4 = c([1;2;3;4], i-1) + real( I(4)' * ( 2*u(i)-(5/6)*I(4)-I(foundPath) ) );
        [c(4,i), maxIndex4] = max(cCalc4);
        
        pathTemp(1,:) = foundPath(maxIndex1,:);
        pathTemp(1,i) = 1;
        pathTemp(2,:) = foundPath(maxIndex2,:);
        pathTemp(2,i) = 2;
        pathTemp(3,:) = foundPath(maxIndex3,:);
        pathTemp(3,i) = 3;
        pathTemp(4,:) = foundPath(maxIndex4,:);
        pathTemp(4,i) = 4;
        
        foundPath = pathTemp;
        pathTemp = [];
        
      
    else
        cCalc1 = c([1;2;3;4], i-1) + real( I(1)' * ( 2*u(i)-(5/6)*I(1)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(1,i), maxIndex1] = max(cCalc1);
        
        cCalc2 = c([1;2;3;4], i-1) + real( I(2)' * ( 2*u(i)-(5/6)*I(2)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(2,i), maxIndex2] = max(cCalc2);
        
        cCalc3 = c([1;2;3;4], i-1) + real( I(3)' * ( 2*u(i)-(5/6)*I(3)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(3,i), maxIndex3] = max(cCalc3);
        
        cCalc4 = c([1;2;3;4], i-1) + real( I(4)' * ( 2*u(i)-(5/6)*I(4)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(4,i), maxIndex4] = max(cCalc4);
        
        pathTemp(1,:) = foundPath(maxIndex1,:);
        pathTemp(1,i) = 1;
        pathTemp(2,1:(i-1)) = foundPath(maxIndex2,:);
        pathTemp(2,i) = 2;
        pathTemp(3,1:(i-1)) = foundPath(maxIndex3,:);
        pathTemp(3,i) = 3;
        pathTemp(4,1:(i-1)) = foundPath(maxIndex4,:);
        pathTemp(4,i) = 4;
        
        foundPath = pathTemp;
        pathTemp = [];
        
    end
end
[maxCost, maxCostIndex] = max(c(:,12));
foundBits(1,1) = I(foundPath(maxCostIndex,1));

for i=13:numOfBits
    
    cCalc1 = c([1;2;3;4], i-1) + real( I(1)' * ( 2*u(i)-(5/6)*I(1)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
    [c(1,i), maxIndex1] = max(cCalc1);

    cCalc2 = c([1;2;3;4], i-1) + real( I(2)' * ( 2*u(i)-(5/6)*I(2)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
    [c(2,i), maxIndex2] = max(cCalc2);

    cCalc3 = c([1;2;3;4], i-1) + real( I(3)' * ( 2*u(i)-(5/6)*I(3)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
    [c(3,i), maxIndex3] = max(cCalc3);

    cCalc4 = c([1;2;3;4], i-1) + real( I(4)' * ( 2*u(i)-(5/6)*I(4)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
    [c(4,i), maxIndex4] = max(cCalc4);

    pathTemp(1,:) = foundPath(maxIndex1,:);
    pathTemp(1,i) = 1;
    pathTemp(2,1:(i-1)) = foundPath(maxIndex2,:);
    pathTemp(2,i) = 2;
    pathTemp(3,1:(i-1)) = foundPath(maxIndex3,:);
    pathTemp(3,i) = 3;
    pathTemp(4,1:(i-1)) = foundPath(maxIndex4,:);
    pathTemp(4,i) = 4;

    foundPath = pathTemp;
    pathTemp = [];
    
    [maxCost, maxCostIndex] = max(c(:,i));
    foundBits(i-11,1) = I(foundPath(maxCostIndex,i-11));
end
pe = nnz(foundBits(1:end-11,1)-bits(1:end-11,1))/(numOfBits-11);
fprintf('Probability of error of QPSK modulation without noise is %d as expected. \n' ,pe);


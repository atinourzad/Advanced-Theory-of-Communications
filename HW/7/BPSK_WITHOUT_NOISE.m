%% BPSK Without Noise for Checking Viterbi Algorithm: First Way

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Bits and Sending Them
numOfBits = 10^4;

%Creating -1 and +1 bits:
bits = 2 * randi([0,1], [numOfBits,1])-1;

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
c = zeros(2, numOfBits);
foundBits = zeros(numOfBits,1);
I = [1;-1];

for i=1:12
    if i==1
        c(1,i) = real( I(1)' * ( 2*u(i)-(5/6)*I(1) ) );
        c(2,i) = real( I(2)' * ( 2*u(i)-(5/6)*I(2) ) );

        foundPath = [1; 2];

    elseif i==2
        cCalc1 = c([1;2], i-1) + real( I(1)' * ( 2*u(i)-(5/6)*I(1)-I(foundPath) ) );
        [c(1,i), maxIndex1] = max(cCalc1);

        cCalc2 = c([1;2], i-1) + real( I(2)' * ( 2*u(i)-(5/6)*I(2)-I(foundPath) ) );
        [c(2,i), maxIndex2] = max(cCalc2);

        pathTemp(1,:) = foundPath(maxIndex1,:);
        pathTemp(1,i) = 1;
        pathTemp(2,:) = foundPath(maxIndex2,:);
        pathTemp(2,i) = 2;

        foundPath = pathTemp;
        pathTemp = [];

    else
        cCalc1 = c([1;2], i-1) + real( I(1)' * ( 2*u(i)-(5/6)*I(1)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(1,i), maxIndex1] = max(cCalc1);

        cCalc2 = c([1;2], i-1) + real( I(2)' * ( 2*u(i)-(5/6)*I(2)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(2,i), maxIndex2] = max(cCalc2);

        pathTemp(1,:) = foundPath(maxIndex1,:);
        pathTemp(1,i) = 1;
        pathTemp(2,1:(i-1)) = foundPath(maxIndex2,:);
        pathTemp(2,i) = 2;

        foundPath = pathTemp;
        pathTemp = [];
    end
end
[maxCost, maxCostIndex] = max(c(:,12));
foundBits(1,1) = I(foundPath(maxCostIndex,1));

for i=13:numOfBits
    cCalc1 = c([1;2], i-1) + real( I(1)' * ( 2*u(i)-(5/6)*I(1)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
    [c(1,i), maxIndex1] = max(cCalc1);

    cCalc2 = c([1;2], i-1) + real( I(2)' * ( 2*u(i)-(5/6)*I(2)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
    [c(2,i), maxIndex2] = max(cCalc2);

    pathTemp(1,:) = foundPath(maxIndex1,:);
    pathTemp(1,i) = 1;
    pathTemp(2,1:(i-1)) = foundPath(maxIndex2,:);
    pathTemp(2,i) = 2;

    foundPath = pathTemp;
    pathTemp = [];

    [maxCost, maxCostIndex] = max(c(:,i));
    foundBits(i-11,1) = I(foundPath(maxCostIndex,i-11));
end
pe = nnz(foundBits(1:end-11)-bits(1:end-11))/(numOfBits-11);
fprintf('Probability of error of BPSK modulation without noise is %d as expected. \n' ,pe);


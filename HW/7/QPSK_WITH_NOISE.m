clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Bits and Sending Them
numOfBits = 10^4;

%Creating -1 and +1 bits:
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


%Simulating h:
step = 0.01;
f = @(t) (heaviside(t-2)-heaviside(t-3)).*((3-t)/sqrt(2)) + (heaviside(t-1)-heaviside(t-2))*(1/sqrt(2)) + (heaviside(t)-heaviside(t-1)).*((t)/sqrt(2));
h = f(0:step:3);

%Creating noise:
nMean = 0;
nSD = sqrt(0.5);
n = normrnd(nMean, nSD, [numOfBits*1/step,1]) + 1i * normrnd(nMean, nSD, [numOfBits*1/step,1]); %White Gaussian Noise

%White Gaussian Noise Goes Through channel h(-t):
v = conv(n, h)*step*10;
v = v(round((1:numOfBits)/step+1));
fprintf("Since variance of lowpass noise is 1, we expect variance of v be x_0 = 5/6, which is %f. So it's correct. \n", std(v)^2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Detecting:

%Viterbi Algorithm and Detecting Bits:
c = zeros(4, numOfBits);
I = [(1/sqrt(2))+1i*(1/sqrt(2)); (-1/sqrt(2))+1i*(1/sqrt(2)); (1/sqrt(2))-1i*(1/sqrt(2)); (-1/sqrt(2))-1i*(1/sqrt(2))];
foundBits = zeros(numOfBits, 1);
pe = zeros(21,1);

for k = 0:20
    %Creating yk:
    y = u + sqrt(0.72216*10^(-k/10))*v;
    
    for i=1:12
        if i==1
            c(1,i) = real( I(1)' * ( 2*y(i)-(5/6)*I(1) ) );
            c(2,i) = real( I(2)' * ( 2*y(i)-(5/6)*I(2) ) );
            c(3,i) = real( I(3)' * ( 2*y(i)-(5/6)*I(3) ) );
            c(4,i) = real( I(4)' * ( 2*y(i)-(5/6)*I(4) ) );

            foundPath = [1; 2; 3; 4];

        elseif i==2
            cCalc1 = c([1;2;3;4], i-1) + real( I(1)' * ( 2*y(i)-(5/6)*I(1)-I(foundPath) ) );
            [c(1,i), maxIndex1] = max(cCalc1);

            cCalc2 = c([1;2;3;4], i-1) + real( I(2)' * ( 2*y(i)-(5/6)*I(2)-I(foundPath) ) );
            [c(2,i), maxIndex2] = max(cCalc2);

            cCalc3 = c([1;2;3;4], i-1) + real( I(3)' * ( 2*y(i)-(5/6)*I(3)-I(foundPath) ) );
            [c(3,i), maxIndex3] = max(cCalc3);

            cCalc4 = c([1;2;3;4], i-1) + real( I(4)' * ( 2*y(i)-(5/6)*I(4)-I(foundPath) ) );
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
            cCalc1 = c([1;2;3;4], i-1) + real( I(1)' * ( 2*y(i)-(5/6)*I(1)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
            [c(1,i), maxIndex1] = max(cCalc1);

            cCalc2 = c([1;2;3;4], i-1) + real( I(2)' * ( 2*y(i)-(5/6)*I(2)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
            [c(2,i), maxIndex2] = max(cCalc2);

            cCalc3 = c([1;2;3;4], i-1) + real( I(3)' * ( 2*y(i)-(5/6)*I(3)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
            [c(3,i), maxIndex3] = max(cCalc3);

            cCalc4 = c([1;2;3;4], i-1) + real( I(4)' * ( 2*y(i)-(5/6)*I(4)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
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
        cCalc1 = c([1;2;3;4], i-1) + real( I(1)' * ( 2*y(i)-(5/6)*I(1)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(1,i), maxIndex1] = max(cCalc1);

        cCalc2 = c([1;2;3;4], i-1) + real( I(2)' * ( 2*y(i)-(5/6)*I(2)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(2,i), maxIndex2] = max(cCalc2);

        cCalc3 = c([1;2;3;4], i-1) + real( I(3)' * ( 2*y(i)-(5/6)*I(3)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
        [c(3,i), maxIndex3] = max(cCalc3);

        cCalc4 = c([1;2;3;4], i-1) + real( I(4)' * ( 2*y(i)-(5/6)*I(4)-I(foundPath(:,i-1))-(1/6)*I(foundPath(:,i-2)) ) );
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
    pe(k+1,1) = nnz(foundBits(1:end-11)-bits(1:end-11))/(numOfBits-11);
end

semilogy(0:20, pe, 'Color', [19/255, 206/255, 188/255]);
xlim([0 20]);
xlabel('SNR (dB)');
ylabel('p_e');
title('Probability of Error for Detectiong QPSK Noisy Bits Using Viterbi Algorithm');
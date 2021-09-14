%% MMSE Equalizer:

clc
clear
close all

%Creating Bits:
numOfBits = 10^4;

%Creating -1 and +1 bits:
bits = 2 * randi([0,1], [numOfBits,1])-1;

%Simulating h:
step = 0.01;
f = @(t) (heaviside(t-2)-heaviside(t-3)).*((3-t)/sqrt(2)) + (heaviside(t-1)-heaviside(t-2))*(1/sqrt(2)) + (heaviside(t)-heaviside(t-1)).*((t)/sqrt(2));
h = f(0:step:3);

%Creating noise:
nMean = 0;
nSD = sqrt(0.5);
n1 = normrnd(nMean, nSD, [numOfBits*1/step,1]) + 1i * normrnd(nMean, nSD, [numOfBits*1/step,1]); %White Gaussian Noise
n2 = normrnd(nMean, nSD, [numOfBits,1]) + 1i * normrnd(nMean, nSD, [numOfBits,1]); %White Gaussian Noise

%White Gaussian Noise Goes through channel h(-t):
v = conv(n1, h)*step*10;
v = v(round((1:numOfBits)/step+1));

%Bits Go through Channel x(t):
x = [1/12 1/2 5/6 1/2 1/12];
uConv = conv(bits,x);
uConv = uConv(3:end-2);

%% For 5 Taps:

%Given Parameters:
numTaps = 5;
l = (length(x)-1)/2;
k = (numTaps-1)/2;

%Finding Ry(m) and RIy(m) Based on SNR:

%RIy(m):
e = zeros(numTaps,1);
e(k-l+1:k-l+length(x),1) = x';

%Ry(m):
s = conv(x,x);
s = s(length(x):end);
Ry = toeplitz(s);
gamma = zeros(numTaps,numTaps);
gamma(1:length(x), 1:length(x)) = Ry;
for i=1:numTaps
   for j=1:numTaps
       if abs(i-j) < length(x) && gamma(i,j) == 0 
            gamma(i,j) = s(abs(i-j)+1);
       end
   end
end

xTemp = x(ceil(length(x)/2):end);
SNRLen = 16;

pe5_1 = zeros(SNRLen,1);
pe5_2 = zeros(SNRLen,1);
peTheory5 = zeros(SNRLen,1);
SINR5 = zeros(SNRLen,1);

for i = 0:SNRLen-1
    
    %Calculating SNR:
    SNR = 10^(i/10);
    
    %Calculating Gamma Matrix:
    gammaUp = zeros(numTaps,numTaps);
    for m=1:numTaps
       for j=1:numTaps
           if abs(m-j) < length(xTemp)
               gammaUp(m,j) = gamma(m,j) + xTemp(abs(m-j)+1)/SNR;
           else
               gammaUp(m,j) = gamma(m,j);
           end
       end
    end
    
    %Calculating d:
    d = gammaUp\e;
    
    %Calculating q:
    q = conv(d, x);
    qLen = length(q);
    
    %Creating New Bits by Bits Going Through q:
    uQ = conv(bits, q);
    wantedIndex = floor(qLen/2);
    uQ = uQ(wantedIndex+1:end-wantedIndex);
    
    %Calculating Noise Variance:
    noiseFilter = conv(q, d);
    N0 = sum(abs(q).^2)/(2*SNR);
    sigmaEta = 2 * N0 * noiseFilter(ceil(length(noiseFilter)/2));
    eta = sqrt(sigmaEta) * n2;
    
    %Recieved Signal before Equalizer:
    y1 = uConv + v * sqrt(N0);
    %Recieced Signal After Equalizer:
    y1 = conv(y1,d);
    y1 = y1(1+k:numOfBits+k);
    
    %Recieved Signal After q:
    y2 = uQ + eta;
    
    %Detecting Recieved Signal:
    foundBits1 = 2 * (real(y1)>0) - 1;
    foundBits2 = 2 * (real(y2)>0) -1;
    
    %Calculating Probability of Error:
    pe5_1(i+1,1) = nnz(foundBits1-bits)/(numOfBits);
    pe5_2(i+1,1) = nnz(foundBits2-bits)/(numOfBits);
    
    %Caclculating Probability of Error Based on Theory:
    q0 = q(ceil(qLen/2));
    bitCoeff = q(q ~= q0);
    for j = 1:2^(length(bitCoeff))-1
        peTheory5(i+1,1) = peTheory5(i+1,1) + qfunc(sqrt((q0 + sum((2*de2bi(j, length(bitCoeff))-1).*bitCoeff))^2 /(sigmaEta/2)));
    end
    peTheory5(i+1,1) = 1/2^(length(bitCoeff))* peTheory5(i+1,1);
    
    %Calculating SINR:
    SINR5(i+1,1) = q0/(1-q0);
end

%Probability of Error Curves:
figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe5_1, '-*', 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, pe5_2, '-o', 'Color', [19/255, 206/255, 188/255]);
semilogy(SNRdB, peTheory5, '--', 'Color', [128/255, 128/255, 0/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
legend('Simulation 1', 'Simulation 2', 'Theory');
title('Probability of Error Comparision after Using MMSE Equalizer for 5 Taps');

%SINR curve:
figure
semilogy(SNRdB, SINR5, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('SINR');
title('SINR after Using MMSE Equalizer for 5 Taps');

%% For 9 taps:

%Given Parameters:
numTaps = 9;
l = (length(x)-1)/2;
k = (numTaps-1)/2;

%Finding Ry(m) and RIy(m) Based on SNR:

%RIy(m):
e = zeros(numTaps,1);
e(k-l+1:k-l+length(x),1) = x';

%Ry(m):
s = conv(x,x);
s = s(length(x):end);
Ry = toeplitz(s);
gamma = zeros(numTaps,numTaps);
gamma(1:length(x), 1:length(x)) = Ry;
for i=1:numTaps
   for j=1:numTaps
       if abs(i-j) < length(x) && gamma(i,j) == 0 
            gamma(i,j) = s(abs(i-j)+1);
       end
   end
end

xTemp = x(ceil(length(x)/2):end);
SNRLen = 16;

pe9_1 = zeros(SNRLen,1);
pe9_2 = zeros(SNRLen,1);
peTheory9 = zeros(SNRLen,1);
SINR9 = zeros(SNRLen,1);
for i = 0:SNRLen-1
    
    %Calculating SNR:
    SNR = 10^(i/10);
    
    %Calculating Gamma Matrix:
    gammaUp = zeros(numTaps,numTaps);
    for m=1:numTaps
       for j=1:numTaps
           if abs(m-j) < length(xTemp)
               gammaUp(m,j) = gamma(m,j) + xTemp(abs(m-j)+1)/SNR;
           else
               gammaUp(m,j) = gamma(m,j);
           end
       end
    end
    
    %Calculating d:
    d = gammaUp\e;
    
    %Calculating q:
    q = conv(d, x);
    q = q.';
    qLen = length(q);
    
    %Creating New Bits by Bits Going Through q:
    uQ = conv(bits, q);
    wantedIndex = floor(qLen/2);
    uQ = uQ(wantedIndex+1:end-wantedIndex);
    
    %Calculating Noise Variance:
    noiseFilter = conv(q, d);
    N0 = sum(abs(q).^2)/(2*SNR);
    sigmaEta = 2 * N0 * noiseFilter(ceil(length(noiseFilter)/2));
    eta = sqrt(sigmaEta) * n2;
    
    %Recieved Signal before Equalizer:
    y1 = uConv + v * sqrt(N0);
    %Recieced Signal After Equalizer:
    y1 = conv(y1,d);
    y1 = y1(1+k:numOfBits+k);
    
    %Recieved Signal After q:
    y2 = uQ + eta;
    
    %Detecting Recieved Signal:
    foundBits1 = 2 * (real(y1)>0) - 1;
    foundBits2 = 2 * (real(y2)>0) -1;
    
    %Calculating Probability of Error:
    pe9_1(i+1,1) = nnz(foundBits1-bits)/(numOfBits);
    pe9_2(i+1,1) = nnz(foundBits2-bits)/(numOfBits);
    
    %Caclculating Probability of Error Based on Theory:
    q0 = q(ceil(qLen/2));
    bitCoeff = q(q ~= q0);
    for j = 1:2^(length(bitCoeff))-1
        peTheory9(i+1,1) = peTheory9(i+1,1) + qfunc(sqrt((q0 + sum((2*de2bi(j, length(bitCoeff))-1).*bitCoeff))^2 /(sigmaEta/2)));
    end
    peTheory9(i+1,1) = 1/2^(length(bitCoeff))* peTheory9(i+1,1);
    
    %Calculating SINR:
    SINR9(i+1,1) = q0/(1-q0);
end

%Probability of Error Curves:
figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe9_1, '-*', 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, pe9_2, '-o', 'Color', [19/255, 206/255, 188/255]);
semilogy(SNRdB, peTheory9, '--', 'Color', [128/255, 128/255, 0/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
legend('Simulation 1', 'Simulation 2', 'Theory');
title('Probability of Error Comparision after Using MMSE Equalizer for 9 Taps');

%SINR curve:
figure
semilogy(SNRdB, SINR9, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('SINR');
title('SINR after Using MMSE Equalizer for 9 Taps');

%% For 13 taps:

%Given Parameters:
numTaps = 13;
l = (length(x)-1)/2;
k = (numTaps-1)/2;

%Finding Ry(m) and RIy(m) Based on SNR:

%RIy(m):
e = zeros(numTaps,1);
e(k-l+1:k-l+length(x),1) = x';

%Ry(m):
s = conv(x,x);
s = s(length(x):end);
Ry = toeplitz(s);
gamma = zeros(numTaps,numTaps);
gamma(1:length(x), 1:length(x)) = Ry;
for i=1:numTaps
   for j=1:numTaps
       if abs(i-j) < length(x) && gamma(i,j) == 0 
            gamma(i,j) = s(abs(i-j)+1);
       end
   end
end

xTemp = x(ceil(length(x)/2):end);
SNRLen = 16;

pe13_1 = zeros(SNRLen,1);
pe13_2 = zeros(SNRLen,1);
peTheory13 = zeros(SNRLen,1);
SINR13 = zeros(SNRLen,1);

for i = 0:SNRLen-1
    
    %Calculating SNR:
    SNR = 10^(i/10);
    
    %Calculating Gamma Matrix:
    gammaUp = zeros(numTaps,numTaps);
    for m=1:numTaps
       for j=1:numTaps
           if abs(m-j) < length(xTemp)
               gammaUp(m,j) = gamma(m,j) + xTemp(abs(m-j)+1)/SNR;
           else
               gammaUp(m,j) = gamma(m,j);
           end
       end
    end
    
    %Calculating d:
    d = gammaUp\e;
    
    %Calculating q:
    q = conv(d, x);
    q = q.';
    qLen = length(q);
    
    %Creating New Bits by Bits Going Through q:
    uQ = conv(bits, q);
    wantedIndex = floor(qLen/2);
    uQ = uQ(wantedIndex+1:end-wantedIndex);
    
    %Calculating Noise Variance:
    noiseFilter = conv(q, d);
    N0 = sum(abs(q).^2)/(2*SNR);
    sigmaEta = 2 * N0 * noiseFilter(ceil(length(noiseFilter)/2));
    eta = sqrt(sigmaEta) * n2;
    
    %Recieved Signal before Equalizer:
    y1 = uConv + v * sqrt(N0);
    %Recieced Signal After Equalizer:
    y1 = conv(y1,d);
    y1 = y1(1+k:numOfBits+k);
    
    %Recieved Signal After q:
    y2 = uQ + eta;
    
    %Detecting Recieved Signal:
    foundBits1 = 2 * (real(y1)>0) - 1;
    foundBits2 = 2 * (real(y2)>0) -1;
    
    %Calculating Probability of Error:
    pe13_1(i+1,1) = nnz(foundBits1-bits)/(numOfBits);
    pe13_2(i+1,1) = nnz(foundBits2-bits)/(numOfBits);
    
    %Caclculating Probability of Error Based on Theory:
    q0 = q(ceil(qLen/2));
    bitCoeff = q(q ~= q0);
    for j = 1:2^(length(bitCoeff))-1
        peTheory13(i+1,1) = peTheory13(i+1,1) + qfunc(sqrt((q0 + sum((2*de2bi(j, length(bitCoeff))-1).*bitCoeff))^2 /(sigmaEta/2)));
    end
    peTheory13(i+1,1) = 1/2^(length(bitCoeff))* peTheory13(i+1,1);

    %SINR Calculation:
    SINR13(i+1,1) = q0/(1-q0);
end

%Probability of Error Curves:
figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe13_1, '-*', 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, pe13_2, '-o', 'Color', [19/255, 206/255, 188/255]);
semilogy(SNRdB, peTheory13, '--', 'Color', [128/255, 128/255, 0/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
legend('Simulation 1', 'Simulation 2', 'Theory');
title('Probability of Error Comparision after Using MMSE Equalizer for 13 Taps');

%SINR curve:
figure
semilogy(SNRdB, SINR13, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('SINR');
title('SINR after Using MMSE Equalizer for 13 Taps');

%% Comparing Different Taps Effects:

figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe5_1, '-*', 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, pe9_1, '-o', 'Color', [19/255, 206/255, 188/255]);
semilogy(SNRdB, pe13_1, '--', 'Color', [128/255, 128/255, 0/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
legend('5 Taps', '9 Taps', '13 Taps');
title('Probability of Error Comparision after Using MMSE Equalizer for Different Taps Using Way1');

figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe5_2, '-*', 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, pe9_2, '-o', 'Color', [19/255, 206/255, 188/255]);
semilogy(SNRdB, pe13_2, '--', 'Color', [128/255, 128/255, 0/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
legend('5 Taps', '9 Taps', '13 Taps');
title('Probability of Error Comparision after Using MMSE Equalizer for Different Taps Using Way2');

figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, SINR5, '-*', 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, SINR9, '-o', 'Color', [19/255, 206/255, 188/255]);
semilogy(SNRdB, SINR13, '--', 'Color', [128/255, 128/255, 0/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
legend('5 Taps', '9 Taps', '13 Taps');
title('SINR Comparision after Using MMSE Equalizer for Different Taps');

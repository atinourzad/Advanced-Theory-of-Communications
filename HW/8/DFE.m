%% With DFE Equalizer:

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
u = conv(bits,x);
u = u(3:end-2);

%% For 5 Taps:

%Given Parameters:
numTaps = 5;
l = (length(x)-1)/2;
k = (numTaps-1)/2;

%Finding e and gamma Based on SNR:

%e:
e = zeros(numTaps,1);
e(k-l+1:k+1,1) = x(1:l+1)';

%gamma:

%Ry:
s = conv(x,x);
s = s(length(x):end);
Ry = toeplitz(s);
endIndexR = min(k+1, length(s));
endIndexRMat = min(k+1, length(x));
RyMat = zeros(k+1,k+1);
RyMat(1:endIndexRMat,1:endIndexRMat) = Ry(1:endIndexR,1:endIndexR);
for i=1:k+1
   for j=1:k+1
       if abs(i-j) < length(x) && RyMat(i,j) == 0 
            RyMat(i,j) = s(abs(i-j)+1);
       end
   end
end

%RIy:
RIyMat = zeros(k+1,k);
xCut = x(ceil(length(x)/2):end);
for i=1:k+1
   for j=1:k
       if k+1-i+j < length(xCut)
           RIyMat(i,j) = xCut(k+2-i+j);
       end
   end
end

%RI:
RIMat = eye(k,k);

SNRLen = 16;

pe5 = zeros(SNRLen,1);
SINR5 = zeros(SNRLen,1);

for i = 0:SNRLen-1
    
    %Calculating SNR:
    SNR = 10^(i/10);

    %Calculating Gamma Matrix:
    RyMatUp = zeros(k+1,k+1);
    for m=1:k+1
       for j=1:k+1
           if abs(m-j) < length(xCut)
               RyMatUp(m,j) = RyMat(m,j) + xCut(abs(m-j)+1)/SNR;
           else
               RyMatUp(m,j) = RyMat(m,j);
           end
       end
    end
    gamma = zeros(numTaps,numTaps);
    gamma(1:k+1,1:k+1) = RyMatUp;
    gamma(1:k+1,k+2:end) = RIyMat;
    gamma(k+2:end,1:k+1) = RIyMat';
    gamma(k+2:end,k+2:end) = RIMat;
    
    %Calculating d:
    d = gamma\e;
    dFF = d(1:k+1);%FeedForward 
    dFB = d(k+2:end);%FeedBack
    
    %Calculating q:
    q = conv(x, dFF);
    qLen = length(q);
    
    %Creating Bits going through Feedback:
    uF = conv(bits, dFB);
    uF = uF(1:numOfBits-1);
    
    %Calculating Noise Variance:
    N0 = sum(abs(q).^2)/(2*SNR);
    
    %Recieved Signal:
    y = u + v * sqrt(N0);
    y = conv(y, dFF);
    y = y(1+k:numOfBits+k);

    result = zeros(numOfBits,1);
    result(1)=y(1);
    result(2:end) = y(2:end) + uF;
    y = result;
    
    %Detecting Recieved Signal:
    foundBits = 2 * (real(y)>0) -1;
    
    %Calculating Probability of Error:
    pe5(i+1,1) = nnz(foundBits-bits)/(numOfBits);
    
    %Calculating SINR:
    q0 = q(k+3);
    SINR5(i+1,1) = q0/(1-q0);
end

figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe5, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
title('Probability of Error after Using DFE Equalizer for 5 Taps');

%SINR curve:
figure
semilogy(SNRdB, SINR5, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('SINR');
title('SINR after Using DFE Equalizer for 5 Taps');

%% For 9 Taps:

%Given Parameters:
numTaps = 9;
l = (length(x)-1)/2;
k = (numTaps-1)/2;

%Finding e and gamma Based on SNR:

%e:
e = zeros(numTaps,1);
e(k-l+1:k+1,1) = x(1:l+1)';

%gamma:

%Ry:
s = conv(x,x);
s = s(length(x):end);
Ry = toeplitz(s);
endIndexR = min(k+1, length(s));
endIndexRMat = min(k+1, length(x));
RyMat = zeros(k+1,k+1);
RyMat(1:endIndexRMat,1:endIndexRMat) = Ry(1:endIndexR,1:endIndexR);
for i=1:k+1
   for j=1:k+1
       if abs(i-j) < length(x) && RyMat(i,j) == 0 
            RyMat(i,j) = s(abs(i-j)+1);
       end
   end
end

%RIy:
RIyMat = zeros(k+1,k);
xCut = x(ceil(length(x)/2):end);
for i=1:k+1
   for j=1:k
       if k+1-i+j < length(xCut)
           RIyMat(i,j) = xCut(k+2-i+j);
       end
   end
end

%RI:
RIMat = eye(k,k);

SNRLen = 16;

pe9 = zeros(SNRLen,1);
SINR9 = zeros(SNRLen,1);

for i = 0:SNRLen-1
    
    %Calculating SNR:
    SNR = 10^(i/10);

    %Calculating Gamma Matrix:
    RyMatUp = zeros(k+1,k+1);
    for m=1:k+1
       for j=1:k+1
           if abs(m-j) < length(xCut)
               RyMatUp(m,j) = RyMat(m,j) + xCut(abs(m-j)+1)/SNR;
           else
               RyMatUp(m,j) = RyMat(m,j);
           end
       end
    end
    gamma = zeros(numTaps,numTaps);
    gamma(1:k+1,1:k+1) = RyMatUp;
    gamma(1:k+1,k+2:end) = RIyMat;
    gamma(k+2:end,1:k+1) = RIyMat';
    gamma(k+2:end,k+2:end) = RIMat;
    
    %Calculating d:
    d = gamma\e;
    dFF = d(1:k+1);%FeedForward 
    dFB = d(k+2:end);%FeedBack
    
    %Calculating q:
    q = conv(x, dFF);
    qLen = length(q);
    
    %Creating Bits going through Feedback:
    uF = conv(bits, dFB);
    uF = uF(1:numOfBits-1);
    
    %Calculating Noise Variance:
    N0 = sum(abs(q).^2)/(2*SNR);
    
    %Recieved Signal:
    y = u + v * sqrt(N0);
    y = conv(y, dFF);
    y = y(1+k:numOfBits+k);

    result = zeros(numOfBits,1);
    result(1)=y(1);
    result(2:end) = y(2:end) + uF;
    y = result;
    
    %Detecting Recieved Signal:
    foundBits = 2 * (real(y)>0) -1;
    
    %Calculating Probability of Error:
    pe9(i+1,1) = nnz(foundBits-bits)/(numOfBits);
    
    %Calculating SINR:
    q0 = q(k+3);
    SINR9(i+1,1) = q0/(1-q0);
end

figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe9, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
title('Probability of Error after Using DFE Equalizer for 9 Taps');

%SINR curve:
figure
semilogy(SNRdB, SINR9, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('SINR');
title('SINR after Using DFE Equalizer for 9 Taps');

%% For 13 Taps:

%Given Parameters:
numTaps = 13;
l = (length(x)-1)/2;
k = (numTaps-1)/2;

%Finding e and gamma Based on SNR:

%e:
e = zeros(numTaps,1);
e(k-l+1:k+1,1) = x(1:l+1)';

%gamma:

%Ry:
s = conv(x,x);
s = s(length(x):end);
Ry = toeplitz(s);
endIndexR = min(k+1, length(s));
endIndexRMat = min(k+1, length(x));
RyMat = zeros(k+1,k+1);
RyMat(1:endIndexRMat,1:endIndexRMat) = Ry(1:endIndexR,1:endIndexR);
for i=1:k+1
   for j=1:k+1
       if abs(i-j) < length(x) && RyMat(i,j) == 0 
            RyMat(i,j) = s(abs(i-j)+1);
       end
   end
end

%RIy:
RIyMat = zeros(k+1,k);
xCut = x(ceil(length(x)/2):end);
for i=1:k+1
   for j=1:k
       if k+1-i+j < length(xCut)
           RIyMat(i,j) = xCut(k+2-i+j);
       end
   end
end

%RI:
RIMat = eye(k,k);

SNRLen = 16;

pe13 = zeros(SNRLen,1);
SINR13 = zeros(SNRLen,1);

for i = 0:SNRLen-1
    
    %Calculating SNR:
    SNR = 10^(i/10);

    %Calculating Gamma Matrix:
    RyMatUp = zeros(k+1,k+1);
    for m=1:k+1
       for j=1:k+1
           if abs(m-j) < length(xCut)
               RyMatUp(m,j) = RyMat(m,j) + xCut(abs(m-j)+1)/SNR;
           else
               RyMatUp(m,j) = RyMat(m,j);
           end
       end
    end
    gamma = zeros(numTaps,numTaps);
    gamma(1:k+1,1:k+1) = RyMatUp;
    gamma(1:k+1,k+2:end) = RIyMat;
    gamma(k+2:end,1:k+1) = RIyMat';
    gamma(k+2:end,k+2:end) = RIMat;
    
    %Calculating d:
    d = gamma\e;
    dFF = d(1:k+1);%FeedForward 
    dFB = d(k+2:end);%FeedBack
    
    %Calculating q:
    q = conv(x, dFF);
    qLen = length(q);
    
    %Creating Bits going through Feedback:
    uF = conv(bits, dFB);
    uF = uF(1:numOfBits-1);
    
    %Calculating Noise Variance:
    N0 = sum(abs(q).^2)/(2*SNR);
    
    %Recieved Signal:
    y = u + v * sqrt(N0);
    y = conv(y, dFF);
    y = y(1+k:numOfBits+k);

    result = zeros(numOfBits,1);
    result(1)=y(1);
    result(2:end) = y(2:end) + uF;
    y = result;
    
    %Detecting Recieved Signal:
    foundBits = 2 * (real(y)>0) -1;
    
    %Calculating Probability of Error:
    pe13(i+1,1) = nnz(foundBits-bits)/(numOfBits);
    
    %Calculating SINR:
    q0 = q(k+3);
    SINR13(i+1,1) = q0/(1-q0);
end

figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe13, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
title('Probability of Error after Using DFE Equalizer for 13 Taps');

%SINR curve:
figure
semilogy(SNRdB, SINR13, 'Color', [19/255, 206/255, 188/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('SINR');
title('SINR after Using DFE Equalizer for 13 Taps');

%% Comparing Different Taps Effects:

figure
SNRdB = 0:SNRLen-1;
semilogy(SNRdB, pe5, '-*', 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, pe9, '-o', 'Color', [19/255, 206/255, 188/255]);
semilogy(SNRdB, pe13, '--', 'Color', [128/255, 128/255, 0/255]);
xlim([0 SNRLen-1]);
xlabel('SNR (dB)');
ylabel('p_e');
legend('5 Taps', '9 Taps', '13 Taps');
title('Probability of Error Comparision after Using DFE Equalizer for Different Taps');

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
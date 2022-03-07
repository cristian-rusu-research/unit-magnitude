close all
clear
clc

%%% sizes
Nt = 128;
Lt = 10;
Nrf = 10;

%%% generate random X
[Q, ~] = qr(randn(Nt)+1i*randn(Nt));
Q = Q(:, 1:Nrf);

%%% to normalize, or not to normalize, that is the question
normalize = 0;

%%% run the two methods
tic; [Frf, Fbb, error] = hd_lsr(Q, Lt, normalize); time = toc;
hold on; plot(0:length(error)-1, error*100, '--ro');
tic; [Frf_extended, Fbb_extended, error_extended] = hd_lsr_extended(Q, Lt, normalize); time_extended = toc;
hold on; plot(0:length(error_extended)-1, error_extended*100, '--bx');
grid on; box on;
xlabel('iteration'); ylabel('objective function error');

LL = min(Lt, Nrf);

%%% check if results is correct, see Remark 1 from the paper
QQ = custom_product(Frf_extended(:, 1:LL) + conj(Frf_extended(:, LL+1:end)), Fbb_extended(1:LL, :), abs(Frf_extended(:, LL+1:end)) == zeros(Nt, LL)); % norm(Frf6*Fbb6 - QQ, 'fro') will be approx. zero

%%% using the notation from the paper
Z1 = Frf_extended(:, 1:LL);
Z2 = Frf_extended(:, LL+1:end);
Y = Fbb_extended(1:LL, :);
% norm(Z1*Y + conj(conj(Z2)*Y) - Frf6*Fbb6, 'fro') will be approx. zero

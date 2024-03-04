load("./5Vpp.mat");

clearvars -except signal_acq;
close all;
clc;

%% Dati

% Limit signal_acq length
%signal_acq = signal_acq(1:2^14);

% frequenza campionamento
fs = 10^5;
% Stima frequenza di partenza attraverso l'FFt del segnale acquisito (signal_acq)
N = length(signal_acq);
frequencies = (0:N-1)*(fs/N);
fft_signal_acq = fft(signal_acq);
magnitude_spectrum = abs(fft_signal_acq(1:N/2+1));
[~, dominant_index] = max(magnitude_spectrum);
f0 = frequencies(dominant_index); % valore iniziale di f0 basato sulla frequenza dominante

% vettore tempo
t_vector = 0:1/fs:(length(signal_acq)-1)/fs;

%% Stima a 3 fattori iniziale basata sulla conoscenza di f0
D0 = [cos(2*pi*f0*t_vector)', sin(2*pi*f0*t_vector)', ones(length(signal_acq),1)];

s0 = D0\signal_acq;
A0_est = s0(1);
B0_est = s0(2);
C0_est = s0(3);

%A_est = sqrt(A0_est^2 + B0_est^2);
%phi_est = atan2(B0_est, A0_est);

signal_est_3params = A0_est*cos(2*pi*f0*t_vector) + B0_est*sin(2*pi*f0*t_vector) + C0_est;

%% Stima a 4 fattori basata sulla conoscenza iniziale di A0 ,B0, C0 ed f0
iterations = 5;

%delta_fi = 0;
% A0
% B0
% C0
% fi= zeros(iterations,1)
fi = f0;
%delta_fi = zeros(iterations, 1);
%dec = zeros(iterations, 1);

for i = 1:iterations
    % Matrice osservazioni
    Di = [cos(2*pi*fi*t_vector)', sin(2*pi*fi*t_vector)', ones(length(signal_acq),1), ...
        (-A0_est*t_vector.*sin(2*pi*fi*t_vector)+B0_est*t_vector.*cos(2*pi*fi*t_vector))'];

    % least squares
    si = Di\signal_acq;
    
    % Aggiornamento parametri
    A0_est = si(1);
    B0_est = si(2);
    C0_est = si(3);
    delta_fi = si(4);
    
    % Aggiornamento della frequenza incrementando progressivamente il numero di cifre significatve 
    fi = fi + delta_fi * 0.1;
    % Calcolo errore
    sig = Di * si;

end

% Il calcolo dell'errore con il metodo sotto genere una mole di dati
% maggiore rispetto a sqrt(sum((signal_acq - sig).^2));
%error = sqrt(sum((signal_acq - (A0_est*cos(2*pi*f0*t_vector) + B0_est*sin(2*pi*f0*t_vector) + C0_est )).^2));
% error = sqrt(sum((signal_acq - sig).^2));
% rss = sum(signal_acq - sig).^2;
% tss = sum(signal_acq - mean(sig)).^2;
% r2 = 1 - rss/tss;

%display valori stimati
disp(['Frequenza stimata: ', num2str(fi), ' Hz']);
% disp(['Coefficiente di determinazione: ', num2str(r2*100), '%']);

A_est = sqrt(A0_est^2 + B0_est^2);
phi_est = atan2(B0_est, A0_est);

signal_est_4params = A0_est*cos(2*pi*fi*t_vector) + B0_est*sin(2*pi*fi*t_vector) + C0_est;

%% Plots
figure(1);
plot(t_vector, signal_acq, t_vector, signal_est_4params);
xlim([0 0.1]);
grid on;
legend("acquired", "estimated 4 params");

% figure(2);
% plot(sig-signal_acq);
% grid on;
% legend("error");
% 
% figure(3);
% plot(db(fft(sig-signal_acq)));
% grid on;
% legend("residual signal");

%% Probabilità di errore

% Calculate the error between the acquired signal and the estimated signal
error_signal = sig - signal_acq;

% Create probability plot
figure;
probplot(error_signal);
title('Probability Plot of Errors in Comparison to Gaussian');

% Additional plot for visualization
figure;
subplot(2, 1, 1);
plot(error_signal);
title('Error Signal');
xlabel('Sample');
ylabel('Error');

subplot(2, 1, 2);
plot(db(fft(error_signal)));
title('FFT of Error Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
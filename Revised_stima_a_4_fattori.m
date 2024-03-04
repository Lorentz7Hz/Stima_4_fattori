load("./5Vpp.mat");
%load("./5Vpp_5s.mat");
clearvars -except signal_acq;
close all;
clc;

%% Dati

% frequenza campionamento
fs = 10^5;
% Stima frequenza di partenza attraverso l'FFt del segnale acquisito (signal_acq)
N = length(signal_acq);
frequencies = (0:N-1)*(fs/N);
fft_signal_acq = fft(signal_acq);
magnitude_spectrum = abs(fft_signal_acq(1:N/2+1));
[~, dominant_index] = max(magnitude_spectrum);
f0 = frequencies(dominant_index); % valore iniziale di f0 basato sulla frequenza dominante

best_error = Inf;  % Initialize with a large value
best_params = [];  % Store the best parameters

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

iterations = 6;
delta_fi = 0;
fi = f0;

% A0= zeros(iterations,1);
% B0= zeros(iterations,1);
% C0= zeros(iterations,1);
% fi= zeros(iterations,1);
%delta_fi = zeros(iterations, 1);
%Delta_fi = zeros(iterations, 1);

for i = 1:iterations
    
    % Aggiornamento della frequenza incrementando progressivamente il numero di cifre significatve
    fi = fi + delta_fi*0.1;

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
    
    sig = Di * si;
    error = sqrt(sum((signal_acq - sig).^2)); 
    
    if error < best_error
        best_error = error;
        best_params = [A0_est, B0_est, C0_est, fi];  
    end

end

% After the loop, use the best_params to compute the estimated signal and other results
A0_est = best_params(1);
B0_est = best_params(2);
C0_est = best_params(3);
fi = best_params(4);

% Il calcolo dell'errore con il metodo sotto genere una mole di dati
% maggiore rispetto a sqrt(sum((signal_acq - sig).^2));
%error = sqrt(sum((signal_acq - (A0_est*cos(2*pi*f0*t_vector) + B0_est*sin(2*pi*f0*t_vector) + C0_est )).^2));
%error = sqrt(sum((signal_acq - sig).^2));   

%% Calcolo errore

error_signal = sig - signal_acq;
% rss = sum(signal_acq - sig).^2;
% tss = sum(signal_acq - mean(sig)).^2;
csvwrite('signal_acq.csv', signal_acq);
csvwrite('sig.csv', sig);
% r2 =  1 - rss/tss;

A_est = sqrt(A0_est^2 + B0_est^2);
phi_est = atan2(B0_est, A0_est);
signal_est_4params = A0_est*cos(2*pi*fi*t_vector) + B0_est*sin(2*pi*fi*t_vector) + C0_est;

%display valori stimati
disp(['Ampiezza stimata: ', num2str(A_est)]);
disp(['fase stimata: ', num2str(phi_est)]);
disp(['Frequenza stimata: ', num2str(fi), ' Hz']);
%disp(['Coefficiente di determinazione: ', num2str(r2*100), '%']);

%% Plots

figure;
plot(t_vector, signal_acq, t_vector, signal_est_4params);
xlim([0 0.1]);
grid on;
legend("signal acq", "stima 4 parametri");

figure;
probplot(error_signal);
title('Probability Plot degli errori rispetto ad una gaussiana');

figure;
subplot(2, 1, 1);
plot(error_signal);
title('Segnale di errore');
xlabel('Campioni');
ylabel('Errore');

subplot(2, 1, 2);
plot(db(fft(error_signal)));
title('FFT of Error Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
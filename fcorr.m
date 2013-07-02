function [ris] = fcorr(w)

% author L.Scorzato (2003)
% calcola la funzione di correlazione (sottraendo il valormedio);
% tratta il vettore ciclicamente , usa la trasf. di fourier.
% Sulle matrici esegue la operazione suddetta su ciascuna colonna.
% ris(1,:) contiene le varianze.

row = size(w,1);
dw = w - ones(row,1) * mean(w);
f = fft(dw);
m2f = abs(f).^2;
temp = ifft(m2f);
ris = real(temp)./ row;


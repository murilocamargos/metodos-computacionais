clear
clc
format long

addpath('../../tools');

%% Defini��es iniciais
F = @(x) 1./(x.^4+x.^2+0.9);
%F = @(x) sqrt(abs(x+0.5));
N = 22;

[I, ar] = chebyshev(F, -1, 1, N);

%% Plota a fun��o real e a aproximada
w = [.5 ones(1, N-1) .5];                   % Vetor de pesos
x = -1:0.05:1;                              % Intervalo
yr = F(x);                                  % Fun��o verdadeira
yt = cos((0:N)'*acos(x))' * (w' .* ar);     % Aproxima��o
hold on;
plot(x,yr,'bx');
plot(x,yt,'ro');
grid on;
hold off;
legend('Fun��o', 'Aproxima��o');

%% Calcula o erro quadr�tico m�dio da aproxima��o
I
rmse = sqrt(mean((yr-yt').^2))
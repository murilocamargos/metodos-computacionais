clear
clc
format long

addpath('../../tools');

%% Definições iniciais
F = @(x) 1./(x.^4+x.^2+0.9);
%F = @(x) sqrt(abs(x+0.5));
N = 22;

[I, ar] = chebyshev(F, -1, 1, N);

%% Plota a função real e a aproximada
w = [.5 ones(1, N-1) .5];                   % Vetor de pesos
x = -1:0.05:1;                              % Intervalo
yr = F(x);                                  % Função verdadeira
yt = cos((0:N)'*acos(x))' * (w' .* ar);     % Aproximação
hold on;
plot(x,yr,'bx');
plot(x,yt,'ro');
grid on;
hold off;
legend('Função', 'Aproximação');

%% Calcula o erro quadrático médio da aproximação
I
rmse = sqrt(mean((yr-yt').^2))
clear
clc
format long

%% Definições iniciais
F = @(x) 1./(x.^4+x.^2+0.9);
%F = @(x) sqrt(abs(x+0.5));
N = 22;

%% Matriz de pesos da soma finita com primeiro e último pesos iguais a 0.5
w = [.5 ones(1, N-1) .5];

%% Cria iteradores de 0 a N
[s,r] = meshgrid(0:N, 0:N);

%% Calcula os valores dos coeficientes a_r
ar = F(cos(s * pi / N)) .* cos(pi * s .* r / N);
ar = (2/N) * ar*w';

%% Calcula os valores dos coeficientes b_r
r = 2:N;
br = (ar(r-1) - ar(r+1))./(2*(r-1)');

%% Calcula o valor da integral
I = 2*sum(br(1:2:end))

%% Plota a função real e a aproximada
x = -1:0.05:1;
yr = F(x);
yt = cos((0:N)'*acos(x))' * (w' .* ar);
hold on;
plot(x,yr,'bx');
plot(x,yt,'ro');
grid on;
hold off;
legend('Função', 'Aproximação');

%% Calcula o erro quadrático mádio da aproximação
rmse = sqrt(mean((yr-yt').^2))
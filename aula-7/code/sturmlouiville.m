% Soluções do problema simples de Sturm-Louiville com diferentes condições
% de contorno.

L = pi;
x = 0:L/100:L;
k_1 = 1;
k_2 = 1;
n = 1:5;

%% Dirichlet, Dirichlet
lambda = (n*pi/L).^2;
u = k_2 * sin(x' * sqrt(lambda));
grafico(x, u);

%% Dirichlet, Neumman
lambda = ((2*n-1)*pi/L).^2;
u = k_2 * sin(x' * sqrt(lambda));
grafico(x, u);

%% Dirichlet, Robin
lambda = ((2*n-1)*pi/L).^2;
u = k_2 * sin(x' * sqrt(lambda));
grafico(x, u);

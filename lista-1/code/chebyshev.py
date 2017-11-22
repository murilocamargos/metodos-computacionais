import numpy as np
import matplotlib.pyplot as plt

# Definicoes iniciais
N = 22
#F = lambda x: 1/(np.power(x, 4) + np.power(x, 2) + 0.9)
F = lambda x: np.sqrt(np.abs(x + 0.5))

# Matriz de pesos da soma finita com primeiro e ultimos pesos iguais a 0.5
w = np.array([0.5] + [1 for i in range(N-1)] + [0.5])

# Cria iteradores de 0 a N
s, r = np.meshgrid(np.arange(N+1), np.arange(N+1))

# Calcula os valores dos coeficientes a_r
ar = F(np.cos(s * np.pi / N)) * np.cos(np.pi * s * r / N)
ar = (2/N) * np.dot(ar, w)

# Calcula os valores dos coeficientes b_r
r = np.arange(1, N)
br = (ar[r-1] - ar[r+1]) / (2*r)

# Calcula o valor da integral
I = 2 * np.sum( br[ np.arange(0, len(br), 2) ] )


# Plota a funcao real e a aproximada
x = np.linspace(-1,1,41)
yr = F(x)
yt = np.cos(np.dot(np.arange(N+1).reshape((N+1, 1)), np.arccos(x).reshape(1, x.shape[0]))).T
yt = np.dot(yt, (ar * w).reshape((N+1, 1)))

# Calcula o erro quadratico medio da aproximacao
rmse = np.mean((yr - yt.T)**2)**0.5

print('Integral:', I)
print('RMSE:', rmse)

# Plota a funcao real e a aproximada
plt.scatter(x, yr, edgecolors='b', marker='o', facecolors='none', label='Função')
plt.scatter(x, yt, color='r', marker='.', label='Aproximação')
plt.grid()
plt.legend()
plt.show()
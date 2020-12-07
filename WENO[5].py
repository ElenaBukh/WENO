# https://pyweno.readthedocs.io/en/latest/tutorial.html !!!

#для решения уравнения переноса используется схема WENO 5-го порядка
import numpy as np
import matplotlib.pyplot as plt
import pylab
import time

N = 1000 #число рассматриваемых точек
S = 3  #число используемых шаблонов (для 5-го порядка используется 3-хточечная схема)
cur = 0.25 # коэффициент Куранта, должен быть меньше 1
x_left = 0
x_mid = 0.2
x_right = 0.4
v = 1
tau = 0
u1 = 1
u2 = 0
dx = (x_right - x_left)/N
dt = cur*dx/v
x = np.zeros(N, dtype = np.float64)
b = np.zeros((N,S), dtype=np.float64)
u = np.zeros((N, S), dtype=np.float64)
eps = 10**-6
# переменные для присваивания значений массивов для отрисовки графиков
u0 = np.zeros(N, dtype = np.float64)
uk = np.zeros(N, dtype = np.float64)

d = np.zeros(S, dtype=np.float64)
w = np.zeros(S, dtype=np.float64)
a = np.zeros(S, dtype=np.float64)
u_center = np.zeros((N,S), dtype=np.float64)
u_sol = np.zeros((N, S), dtype=np.float64)
u_0_1 = np.zeros(N, dtype=np.float64)
u_0_2 = np.zeros(N, dtype=np.float64)
u_0_3 = np.zeros(N, dtype=np.float64)
u_1_1 = np.zeros(N, dtype=np.float64)
u_1_2 = np.zeros(N, dtype=np.float64)
u_1_3 = np.zeros(N, dtype=np.float64)
u_2_1 = np.zeros(N, dtype=np.float64)
u_2_2 = np.zeros(N, dtype=np.float64)
u_2_3 = np.zeros(N, dtype=np.float64)

file = open('weno_solutions.txt', 'w')
starttime = time.time()
#начальные условия 
x[0] = x_left 
for s in range (S):
    for j in range (N-1):
        x[j+1] = x[j] + dx
        if x[j] <= x_mid:
              u[j, s] = u1
        else:
              u[j, s] = u2   
               
#_________________________________________________________________________________________________
# сетка для характеристик, направленных вправо
# индикаторы гладкости
for s in range (S - 1):
    for i in range(N-2):
        if s == 0:
            b[i, 0] = 13/12*((u[i, s] - 2*u[i+1, s] + u[i+2, s]))**2 + 1/4*(3*u[i, s] - 4*u[i+1, s] + u[i+2, s])**2
        if s == 1:
            b[i, 1] = 13/12*((u[i-1, s] - 2*u[i, s] + u[i+1, s]))**2 + 1/4*(u[i-1, s] - u[i+1, s])**2
        if s == 2:
            b[i, 2] = 13/12*((u[i-2, s] - 2*u[i-1, s] + u[i, s]))**2 + 1/4*(u[i-2, s] - 4*u[i-1, s] + 3*u[i, s])**2 

# числовые потоки, значения в середине ячеек
for s in range (S - 1):    
    for i in range (N-1):
        u_0_1[i] = -1/6*u[i, s] 
        u_1_1[i] = 1/3*u[i+1, s] + 5/6*u[i, s] 
        u_2_1[i] = 11/6*u[i, s]

for s in range (S - 1):
    for j in range(N-1):
        u_0_2[j] = 5/6*u[j, s]
        u_1_2[j] = - 1/6*u[j, s]
        u_2_2[j] = - 7/6*u[j, s]

for s in range (S - 1):
    for k in range(N-1):
        u_0_3[k] = 1/3*u[k, s]
        u_1_3[k] = 0
        u_2_3[k] = 1/3*u[k, s]

for z in range (N - 1):            
    u_center[z, 0] = u_0_1[z] + u_0_2[z] + u_0_3[z]
    u_center[z, 1] = u_1_1[z] + u_1_2[z] + u_1_3[z]
    u_center[z, 2] = u_2_1[z] + u_2_2[z] + u_2_3[z]   

# оптимальные коэффициенты
d[0] = 0.3
d[1] = 0.6
d[2] = 0.1

# веса
for s in range(S - 1):
    for i in range(N):
        a[s] = d[s]/(eps + b[i, s])**2

for s in range(S - 1):
    w[s] = a[s]/(a[0] + a[1] + a[2])

# решения для характеристик
for s in range (S - 1):
    for i in range(N):
            u_sol[i] = w[0]*u_center[i, 0] + w[1]*u_center[i, 1] + w[2]*u_center[i, 2]
            file.write(str(round(u_sol[i, s], 5)) + "\n")

#изменение графика в зависимости от временного слоя
tau = cur*((x_right - x_left)/v)
u_sol[0, 0] = u[0, 0]
for s in range (S - 1):
    for i in range (1, N):
        u_sol[i, s] = u_sol[i, s-1] - tau*v*((u_sol[i, s-1] - u_sol[i-1, s-1])/(x_right - x_left))
        x[i] = x[i-1] + dx

#эти циклы для того, чтобы в легенде названия графиков отображались по одному разу
for s in range (S -1):
    for i in range (N):
        u0[i] = u[i, s]

for s in range (S-1):
    for i in range (N):
        uk[i] = u_sol[i, s]

print(time.time() - starttime)

#вывод графика и запись полученных значений в текстовый файл  
pylab.xlim(-0.1, 0.5)
pylab.ylim(-0.1, 1.1)      
line1 = plt.plot (x, uk, 'b-o', label='WENO 5')
line1 = plt.plot (x, u_sol, 'b-o')
line2 = plt.plot (x, u0, 'r-h', label = 'original')
pylab.legend(loc = 'upper right')

file.close()
plt.title('Transport equation (WENO[5] scheme)')
plt.xlabel('Taken points')
plt.ylabel('Solutions of scheme (concentration)')
plt.grid(True)
plt.show()
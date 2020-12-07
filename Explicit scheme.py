#для решения уравнения переноса используется явная схема
import numpy as np
import matplotlib.pyplot as plt
import pylab
import time

N = 1000 #число рассматриваемых точек
T = 3 #временные слои
u1 = 1
a = 1 #скорость изменения функции
u2 = 0
xmin = 0
xmid = 0.2
xmax = 0.4
tau = 0 # изменение по времени
u = np.zeros((N, T), dtype=np.float64)
u_0 = np.zeros((N, T), dtype=np.float64)
x = np.zeros(N, dtype = np.float64)
dx = (xmax - xmin)/N
u0 = np.zeros(N, dtype = np.float64)
uk = np.zeros(N, dtype = np.float64)

#цикл для вычисления входящих значений и заполнения ими массивов
starttime = time.time()
x[0] = xmin
for j in range (N-1):
    x[j] = x[j-1] + dx
    if x[j] <= xmid:
          u[j] = u1
          u_0[j] = u[j]  
    else:
          u[j] = u2       
          u_0[j] = u[j] 

tau = 0.25*((xmax - xmin)/a)

#вычисление значений по явной схеме в зависимости от временного слоя
u[0] = u_0[0]
for t in range (T):
    for i in range (1, N):
        u[i, t] = u[i, t-1] - tau*a*((u[i, t-1]-u[i-1, t-1])/(xmax - xmin))
        x[i] = x[i-1] + dx

#эти циклы для того, чтобы в легенде названия графиков отображались по одному разу
for t in range (T-1):
    for i in range (N):
        u0[i] = u_0[i, t]

for t in range (T-1):
    for i in range (N):
        uk[i] = u[i, t]
        
print(time.time() - starttime)

#вывод графика и запись результатов в текстовый файл
pylab.xlim(-0.1, 0.5)
pylab.ylim(-0.1, 1.1)
plt.plot (x, uk, 'b-o', label='Explicit')
plt.plot (x, u, 'b-o')
plt.plot ( x, u0, 'r-h', label = 'original')
plt.plot ( x, u_0, 'r-h')
pylab.legend(loc = 'upper right')

plt.title('Transport equation (Explicit scheme)')
plt.xlabel('Taken points')
plt.ylabel('Solutions of scheme (concentration)')
plt.grid(True)
plt.show()
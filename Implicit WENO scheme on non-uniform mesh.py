# https://pyweno.readthedocs.io/en/latest/tutorial.html !!!

# for the solving of transport equation WENO5 scheme was used
import numpy as np
import matplotlib.pyplot as plt
import pylab
import time

# non-uniform mesh
w_f_0_Xf_x = 5E-3
x_f_x = 25 
L_x = 100 # length of considered section 
    
x_f_x_Inj1 = x_f_x               # [m] half length of fracture Xf1_Inj1  (on axis Ox)
x_f_x_Prod1 = x_f_x_Inj1         # [m] half length of fracture Xf1_Prod1 (on axis Ox)

x_f_x_Prod2 = x_f_x              # [m] half length of fracture Xf1_Prod2 (on axis Ox)
x_f_x_Prod3 = x_f_x_Prod2        # [m] half length of fracture Xf1_Prod3 (on axis Ox)

D_x = L_x-x_f_x_Prod1 - x_f_x_Prod2 # [m] bridge between fractures in producing wells row
 
Nx = 24*42  # multiple 24, quantity of poins

k100 = 100 # index of mesh smoothness

dx_RU_Prod1 = ((x_f_x_Prod1)-(x_f_x_Prod1/k100))/(6*Nx/24+1)
bb_Xf_x_Prod1 = (2*(x_f_x_Prod1/k100)/2/(Nx/24) - 4*(w_f_0_Xf_x/2))/((Nx/24+0)-1)

dx_RU_D_x = (D_x - 1*(D_x / k100))/(6*Nx/24+1+1)    
bb_D_x = (2*(D_x / k100) / 2/ (Nx / 24) - 4*(w_f_0_Xf_x/2))/((Nx/24+0)-1)

dx_RU_Prod2 = ((x_f_x_Prod2) - (x_f_x_Prod2 / k100))/(6*Nx / 24+1)
bb_Xf_x_Prod2 = (2*(x_f_x_Prod2/k100)/2/(Nx/24) - 4*(w_f_0_Xf_x/2))/((Nx/24+0)-1)

dx  = np.zeros(Nx-1, dtype = np.float64)
x  = np.zeros(Nx-1, dtype = np.float64)

for i in range (1, Nx-1):
    if ((i >= 1) & (i <= 1*Nx/24)):
        dx[i] = (4*w_f_0_Xf_x/2 + bb_Xf_x_Prod1*(i-1))/2*(i)

    elif ((i> Nx/24) & (i< 7*Nx/24+1)):
        dx[i] = dx_RU_Prod1

    elif ((i >= 7*Nx/24+1) & (i <= 8*Nx/24)):
        dx[i] = dx[np.int(i-2*(i-4*Nx/24)+1)]

    elif ((i > 8*Nx/24) & (i < 9*Nx/24)):
        dx[i] = (4*w_f_0_Xf_x/2 +bb_D_x * (i-0-8*Nx/24))/2*(i+1-8*Nx/24)
    
    elif ((i >= 9*Nx/24) & (i < 15*Nx/24+0)):    
        dx[i] = dx_RU_D_x

    elif ((i >= 15*Nx/24+0) & (i<=16*Nx/24)):
        dx[i] = dx[np.int(i-2*(i-12*Nx/24))]

    elif ((i > 16*Nx/24) & (i < 17*Nx/24)):
        dx[i] = (4*w_f_0_Xf_x/2 +bb_Xf_x_Prod2*(i-0-16*Nx/24))/2*(i+1-16*Nx/24)

    elif ((i >= 17*Nx/24) & (i <= 23*Nx/24-2)):    
        dx[i] = dx_RU_Prod2

    elif ((i > 23*Nx/24-2) & (i <= Nx-1)):
        dx[i] = dx[np.int(i-2*(i-20*Nx/24)-1)]

for i in range (1, Nx-1): 
    if i == 1:
        x[i] = 0
    
    elif i == 2:
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1]) + 2*w_f_0_Xf_x/2

    elif ((i > 2) & (i <= 1*Nx/24)): 
        x[i]= x[i-1] + np.abs(dx[i] - dx[i-1])

    elif ((i > 1*Nx/24) & (i <= 7*Nx/24+0)): # the middle
        x[i] = x[i-1] + dx[i]

    elif (i == 7*Nx/24+1):               # the middle + 1
        x[i] = x[i-1] + dx[i-1]

    elif ((i > 7*Nx/24+1) & (i <= 8*Nx/24-1)): 
        x[i] = x[i-1] - (dx[i] - dx[i-1])

    elif ((i == 8*Nx/24)):
        x[i] = x[i-1] + np.abs((dx[i] - dx[i-1])) + 2*w_f_0_Xf_x/2

    elif (i == 8*Nx/24+1):
        x[i] = x[i-1] + np.abs((dx[i] - dx[i-1])) + 2*w_f_0_Xf_x/2

    elif ((i > 8*Nx/24+1) & (i < 9*Nx/24)): 
        x[i] = x[i-1] + np.abs((dx[i] - dx[i-1]))

    elif ((i >= 9*Nx/24) & (i <= 15*Nx/24+0)): # the middle
        x[i] = x[i-1] + dx[i]    

    elif (i == 15*Nx/24+1):               # the middle + 1
        x[i] = x[i-1] + dx[i-1] 

    elif ((i > 15*Nx/24+1) & (i <= 16*Nx/24-1)): 
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1])

    elif (i == 16*Nx/24):
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1]) + 2*w_f_0_Xf_x/2

    elif (i == 16*Nx/24+1):
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1]) + 2*w_f_0_Xf_x/2

    elif ((i > 16*Nx/24+1) & (i < 17*Nx/24-0)): 
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1])

    elif ((i >= 17*Nx/24) & (i <= 23*Nx/24-1)):  # the middle 
        x[i] = x[i-1] + dx[i]

    elif (i == 23*Nx/24):              # the middle + 1
        x[i] = x[i-1] + dx[i-1]

    elif ((i > 23*Nx/24) & (i <= (Nx-1)-1)):
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1])

    elif (i == (Nx-1)):
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1]) + 2*w_f_0_Xf_x/2
# end of non-uniform mesh_______________________________________________________

# Implicit WENO scheme
N = Nx-1 # number of considered poins
S = 3  # number of used stensils (for 5-th order 3-point stensils)
eps = 10**-6
x_left = 0
x_right = L_x
v = 1 # speed of changing graph thouhg the time
u = np.zeros((N, S), dtype=np.float64)
t = np.zeros(N, dtype = np.float64)
t_max = 0
t_low = 0
u_max = 1
u_low = 0
b = np.zeros((N, S), dtype=np.float64)
# variables for assignment of array values for plotting
u0 = np.zeros(N, dtype = np.float64)
uk = np.zeros(N, dtype = np.float64)
# variables used in calculating values of WENO scheme
d = np.zeros(S, dtype=np.float64)
w = np.zeros(S, dtype=np.float64)
a = np.zeros(S, dtype=np.float64)
u_center = np.zeros((N,S), dtype=np.float64)
u_sol = np.zeros((N), dtype=np.float64)
u_0_1 = np.zeros(N, dtype=np.float64)
u_0_2 = np.zeros(N, dtype=np.float64)
u_0_3 = np.zeros(N, dtype=np.float64)
u_1_1 = np.zeros(N, dtype=np.float64)
u_1_2 = np.zeros(N, dtype=np.float64)
u_1_3 = np.zeros(N, dtype=np.float64)
u_2_1 = np.zeros(N, dtype=np.float64)
u_2_2 = np.zeros(N, dtype=np.float64)
u_2_3 = np.zeros(N, dtype=np.float64)
u_pred = np.zeros((N, S), dtype=np.float64)

'''
print('Введите левую границу по оси Х в метрах:')
x_left = float(input())
print('Введите правую границу по оси Х в метрах:')
x_right = float(input())
print('Введите размерность отрезка, на котором будет закачиваться трассер:')
x_size = float(input())
print('Введите значение скорости изменения функции в м/с:')
v = float(input()) 
'''
x_left = 20
x_size = 0.2

t_max = x_size/v
t_step = (t_max - t_low)/N
t[0] = t_low

# values calculation in dummy cells
'''
for s in range (S):
    for i in range (1, 2):
        u[i, s] = (10*(-i + 1))**10
        t[i] = t[i - 1] + t_step
'''
starttime = time.time()
# values calculation in inter cells
for s in range (S):
    for i in range (N - 1):
        t[i] = t[i-1] + t_step
        if x[i] < x_left:
                u[i, s] = u_max
        if x_left <= x[i] <= x_right:
                u[i, s] = u_low

# values calculation in dummy cells
for s in range (S):
    for i in range (N, N):
        u[i, s] = (10*(i + 1))**10
        t[i] = t[i - 1] + t_step
      
#_________________________________________________________________________________________________
# smothness monitors
for s in range (S - 1):
    for i in range(N-2):
        if s == 0:
            b[i, 0] = 13/12*((u[i-2, s] - 2*u[i-1, s] + u[i, s]))**2 + 1/4*(u[i-2, s] - 4*u[i-1, s] + 3*u[i, s])**2
        if s == 1:
            b[i, 1] = 13/12*(((u[i-1, s] + u[i, s]) - 4*u[i, s] + (u[i+1, s] + u[i, s])))**2 + 1/4*(-(u[i-1, s] + u[i, s]) + (u[i+1, s] + u[i, s]))**2
        if s == 2:
            b[i, 2] = 13/12*((u[i, s] - 2*u[i+1, s] + u[i+2, s]))**2 + 1/4*(3*u[i, s] - 4*u[i+1, s] + u[i+2, s])**2 

# numerical flow, values in the midde of the cells
for s in range (S - 1):    
    for i in range (N - 1):
        u_0_1[i] = 2/6*u[i-2, s] 
        u_1_1[i] = -1/6*(u[i-1, s] + u[i, s])
        u_2_1[i] = 2/6*u[i, s]

for s in range (S - 1):
    for j in range(N-1):
        u_0_2[j] = -7/6*u[j-1, s]
        u_1_2[j] = 2/3*u[j, s]
        u_2_2[j] = 5/6*u[j+1, s]

for s in range (S - 1):
    for k in range(N - 2):
        u_0_3[k] = 11/6*u[k, s]
        u_1_3[k] = 2/6*(u[k+1,s] + u[k, s])
        u_2_3[k] = -1/6*u[k+2, s]

for z in range (N - 1):            
    u_center[z, 0] = u_0_1[z] + u_0_2[z] + u_0_3[z]
    u_center[z, 1] = u_1_1[z] + u_1_2[z] + u_1_3[z]
    u_center[z, 2] = u_2_1[z] + u_2_2[z] + u_2_3[z]   

# optimal factors
d[0] = 0.3
d[1] = 0.6
d[2] = 0.1

# weights
for s in range(S - 1):
    for i in range(N):
        a[s] = d[s]/(eps + b[i, s])**2 

for s in range(S - 1):
    w[s] = a[s]/(a[0] + a[1] + a[2])

#_______________________________________________________________________________________________________________________________________________ 
# implicit scheme
for s in range (S - 1):
    for i in range (N - 1):
        u_pred[i, s] = u[i, s] + t_step*w[s]*u_center[i, s]

# smothness monitors
for s in range (S - 1):
    for i in range(N-2):
        if s == 0:
            b[i, 0] = 13/12*((u_pred[i-2, s] - 2*u_pred[i-1, s] + u_pred[i, s]))**2 + 1/4*(u_pred[i-2, s] - 4*u_pred[i-1, s] + 3*u_pred[i, s])**2
        if s == 1:
            b[i, 1] = 13/12*(((u_pred[i-1, s] + u_pred[i, s]) - 4*u_pred[i, s] + (u_pred[i+1, s] + u_pred[i, s])))**2 + 1/4*(-(u_pred[i-1, s] + u_pred[i, s]) + (u_pred[i+1, s] + u_pred[i, s]))**2
        if s == 2:
            b[i, 2] = 13/12*((u_pred[i, s] - 2*u_pred[i+1, s] + u_pred[i+2, s]))**2 + 1/4*(3*u_pred[i, s] - 4*u_pred[i+1, s] + u_pred[i+2, s])**2 

# numerical flow, values in the midde of the cells
for s in range (S - 1):    
    for i in range (N - 1):
        u_0_1[i] = 2/6*u[i-2, s] 
        u_1_1[i] = -1/6*(u[i-1, s] + u[i, s])
        u_2_1[i] = 2/6*u[i, s]

for s in range (S - 1):
    for j in range(N-1):
        u_0_2[j] = -7/6*u[j-1, s]
        u_1_2[j] = 2/3*u[j, s]
        u_2_2[j] = 5/6*u[j+1, s]

for s in range (S - 1):
    for k in range(N - 2):
        u_0_3[k] = 11/6*u[k, s]
        u_1_3[k] = 2/6*(u[k+1,s] + u[k, s])
        u_2_3[k] = -1/6*u[k+2, s]

for z in range (N - 1):            
    u_center[z, 0] = u_0_1[z] + u_0_2[z] + u_0_3[z]
    u_center[z, 1] = u_1_1[z] + u_1_2[z] + u_1_3[z]
    u_center[z, 2] = u_2_1[z] + u_2_2[z] + u_2_3[z]   

# weights
for s in range(S - 1):
    for i in range(N):
        a[s] = d[s]/(eps + b[i, s])**2 

for s in range(S - 1):
    w[s] = a[s]/(a[0] + a[1] + a[2])

u_sol[0] = u[0, 0]
for s in range (S-1):
    for i in range (N-1):
        u_sol[i+1] = u[i, s] + t_step*w[s]*u_center[i+1, s]
#_____________________________________________________________________________________________________________________________________________

# these cycles for displaying names of graphs in the legend once
for s in range (S - 1):
    for i in range (N):
        uk[i] = u_sol[i]

for s in range (S - 1):
    for i in range (N):
        u0[i] = u[i, s]
        
print(time.time() - starttime)

#print('Время, в течение которого рассматривалась закачка трассера:', (round(t_max, 3)), ' секунд.')
#print('Время от начала процесса до окончания закачки трассера:', (round((x_size + x_left)/v, 3)), ' секунд.')

# chart output  
pylab.xlim(- 1, 101)
pylab.ylim(-0.1, 1.2)      
line1 = plt.plot (x, uk, 'b-o', label='WENO 5')
line2 = plt.plot (x, u_sol, 'b-o')
line3 = plt.plot (x, u0, 'r-h', label = 'original')
pylab.legend(loc = 'upper right')
plt.title('Transport equation (WENO[5] scheme)')
plt.xlabel('Taken points (x)')
plt.ylabel('Solutions of scheme, concentration (u)')
plt.grid(True)
plt.show()
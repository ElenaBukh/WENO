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

    elif (i == 7*Nx/24+1):               # middle + 1
        x[i] = x[i-1] + dx[i-1]

    elif ((i > 7*Nx/24+1) & (i <= 8*Nx/24-1)): 
        x[i] = x[i-1] - (dx[i] - dx[i-1])

    elif ((i == 8*Nx/24)):
        x[i] = x[i-1] + np.abs((dx[i] - dx[i-1])) + 2*w_f_0_Xf_x/2

    elif (i == 8*Nx/24+1):
        x[i] = x[i-1] + np.abs((dx[i] - dx[i-1])) + 2*w_f_0_Xf_x/2

    elif ((i > 8*Nx/24+1) & (i < 9*Nx/24)): 
        x[i] = x[i-1] + np.abs((dx[i] - dx[i-1]))

    elif ((i >= 9*Nx/24) & (i <= 15*Nx/24+0)): # middle
        x[i] = x[i-1] + dx[i]    

    elif (i == 15*Nx/24+1):               # middle + 1
        x[i] = x[i-1] + dx[i-1] 

    elif ((i > 15*Nx/24+1) & (i <= 16*Nx/24-1)): 
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1])

    elif (i == 16*Nx/24):
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1]) + 2*w_f_0_Xf_x/2

    elif (i == 16*Nx/24+1):
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1]) + 2*w_f_0_Xf_x/2

    elif ((i > 16*Nx/24+1) & (i < 17*Nx/24-0)): 
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1])

    elif ((i >= 17*Nx/24) & (i <= 23*Nx/24-1)):  # middle 
        x[i] = x[i-1] + dx[i]

    elif (i == 23*Nx/24):              # middle + 1
        x[i] = x[i-1] + dx[i-1]

    elif ((i > 23*Nx/24) & (i <= (Nx-1)-1)):
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1])

    elif (i == (Nx-1)):
        x[i] = x[i-1] + np.abs(dx[i] - dx[i-1]) + 2*w_f_0_Xf_x/2
# end of non-uniform mesh_______________________________________________________

# WENO scheme___________________________________________________________________
starttime = time.time()
N = Nx-1 # number of considered poins
S = 3  # number of used stensils (for 5-th order 3-point stensils)
x_left = 0
x_size = 0
x_right = L_x
v = 1 # speed of changing graph thouhg the time
u = np.zeros(N, dtype=np.float64)
eps = 10**-6
t = np.zeros(N, dtype = np.float64)
t_max = 0
t_low = 0
u_max = 1
u_low = 0
u_sol = np.zeros((N), dtype=np.float64)
# additional factors for calculating final solution
u_sol2 = np.zeros(N, dtype=np.float64)
u_sol3 = np.zeros(N, dtype=np.float64)
u_sol4 = np.zeros(N, dtype=np.float64)
u_sol5 = np.zeros(N, dtype=np.float64)

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
for i in range (1, 2):
    u[i] = (10*(-i + 1))**10
    t[i] = t[i - 1] + t_step
'''
# values calculation in inter cells
for i in range (N - 1):
    t[i] = t[i-1] + t_step
    if x[i] < x_left:
            u[i] = u_max
    if x_left <= x[i] <= x_right:
            u[i] = u_low

# values calculation in dummy cells
for i in range (N, N):
    u[i] = (10*(i + 1))**10
    t[i] = t[i - 1] + t_step

# left flow
for i in range (N - 2):
    u_sol[i] = u[i] + ((x[i+1] - x[i])/2)*((-1/15)*((u[i] - u[i-2])/(x[i] - x[i-2] + eps)) + (11/30)*((u[i] - u[i-1])/(x[i] - x[i-1] + eps)) + (4/5)*((u[i+1] - u[i])/(x[i+1] - x[i] + eps)) - (1/10)*((u[i+2] - u[i])/(x[i+2] - x[i] + eps))) 

u_sol2[0] = u[0]
for i in range (1, N - 1):
    u_sol2[i] = u[i - 1] + ((x[i+1] - x[i])/2)*((-1/15)*((u[i] - u[i-3])/(x[i] - x[i-3] + eps)) + (11/30)*((u[i] - u[i-2])/(x[i] - x[i-2] + eps)) + (4/5)*((u[i] - u[i-1])/(x[i] - x[i-1] + eps)) - (1/10)*((u[i+1] - u[i])/(x[i+1] - x[i] + eps)))

# right flow
for i in range (N - 3):
    u_sol3[i] = u[i+1] - ((x[i+1] - x[i])/2)*((-1/10)*((u[i] - u[i-1])/(x[i] - x[i-1] + eps)) + (4/5)*((u[i+1] - u[i])/(x[i+1] - x[i] + eps)) + (11/30)*((u[i+2] - u[i])/(x[i+2] - x[i] + eps)) - (1/15)*((u[i+3] - u[i])/(x[i+3] - x[i] + eps))) 

for i in range (N - 2):
    u_sol4[i] = u[i] - ((x[i+1] - x[i])/2)*((-1/10)*((u[i] - u[i-2])/(x[i] - x[i-2] + eps)) + (4/5)*((u[i] - u[i-1])/(x[i] - x[i-1] + eps)) + (11/30)*((u[i+1] - u[i])/(x[i+1] - x[i] + eps)) - (1/15)*((u[i+2] - u[i])/(x[i+2] - x[i] + eps))) 

#final flow - summ of left and right flows    
for i in range (N):
    u_sol5[i] = ((u_sol[i] + u_sol2[i] + u_sol3[i] + u_sol4[i])/4)
    
print(time.time() - starttime)

# chart output  
#pylab.xlim(- 1, 101)
#pylab.ylim(-0.1, 0.9) 
line1 = plt.plot (x, u_sol5, 'b-o', label='WENO 5')
line2 = plt.plot (x, u, 'r-h', label = 'exact')
pylab.legend(loc = 'upper right')
plt.title('WENO5 scheme on non-uniform mesh')
plt.xlabel('Taken points')
plt.ylabel('Solutions of scheme, concentration (u)')
plt.grid(True)
plt.show()
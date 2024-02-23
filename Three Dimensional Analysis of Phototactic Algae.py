import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

######## 1. Configuring Data ########

#loads data in time, x, y, z format
data = [t, x, y, z] = np.loadtxt('10x_500fps_0.025msec_setup3B_mid_LED.txt', unpack = True)

#separates data into an array containing only the spatial components
xyz=np.array([data[1],data[2],data[3]])

#changes the frames into a measure of time
#by dividing the frame number by the frame rate
time=[]
a=0

while a<len(data[0]):
    time.append(data[0][a]/500)
    a += 1


#creates a boxcar average of 20 points of the x-direction data points
x_bc=[]
a=0
while a<len(data[1]):
    x_bc.append(sum(data[1][a-20:a+200])/40)
    a += 1


#creates a boxcar average of 20 points of the y-direction data points
y_bc=[]
a=0
while a<len(data[2]):
    y_bc.append(sum(data[2][a-20:a+200])/40)
    a += 1

#creates a boxcar average of 20 points of the z-direction data points
z_bc=[]
a=0
while a<len(data[3]):
    z_bc.append(sum(data[3][a-20:a+200])/40)
    a += 1

#separates the smoothed data into an array containing only the spatial components
xyz_bc=np.arrary([x_bc, y_bc, z_bc])

#turns xyz from a list into an array and generates the transpose of this array
array_xyz=np.array(xyz); transposed_xyz=array_xyz.T


######## 2. Defining Functions ########

#outputs ds from dP and calculates the modulus of dP
#defined in lab book page 26
def arclength(dP):
    return np.sqrt(dP[0]**2+dP[1]**2+dP[2]**2)

#calculates the modulus of dT/ds
#later used as a modulus function for other variables
def modulus(vector_N):
    return np.sqrt((vector_N[0]**2)+(vector_N[1]**2)+(vector_N[2]**2))

#generates the darboux vector using the curvature, torsion, and T and B vectors
def darboux(tau,T,kappa,B):
    return (tau*T)+(kappa*B)


######## 3. Calculating Curvature ########

#creates an empty list for the values of dP
#(the vector displacement between two adjacent points on the curve)
dP=[]

#while loop to produce dP
#dP = position of point n+1 - position of point n
i = 0
while i<len(transposed_xyz)-1:
    dP.append(transposed_xyz[i+1] - transposed_xyz[i])
    i += 1

array_xyz_bc=np.array(xyz_bc); transposed_xyz_bc=array_xyz_bc.T

dP_bc=[]
i = 0
while i<len(transposed_xyz_bc)-1:
    dP_bc.append(transposed_xyz_bc[i+1] - transposed_xyz_bc[i])
    i += 1

#empty list for the arclength ds
#(the scalar distance between two points)
ds = []

#For loop to produce arclengths, ds
for i in dP:
    ds.append(arclength(i))

ds_bc=[]
for i in dP_bc:
    ds_bc.append(arclength(i))

#empty list for T, the tangent unit vector
unit_vector_T=[]

#while loop to produce T, where T=dP/ds
j = 0; d = 0
while j<len(dP)-1:
    while d<len(ds)-1:
        unit_vector_T.append(dP[j]/ds[d])
        j += 1
        d += 1

unit_vector_T_bc=[]
j = 0; d = 0
while j<len(dP_bc)-1:
    while d<len(ds_bc)-1:
        unit_vector_T_bc.append(dP_bc[j]/ds_bc[d])
        j += 1
        d += 1


#empty list for the differential of T, wrt the distance between the points
#as vector and scalar
dT_ds=[]; mod_dT_ds=[]

#produces dT/ds from T and ds
l=0; m=0
while l<len(unit_vector_T)-1:
    while m<len(ds)-1:
        dT_ds.append((unit_vector_T[l-1] - unit_vector_T[l])/ds[m])
        l += 1
        m += 1


dT_ds_bc=[];
l=0; m=0
while l<len(unit_vector_T_bc)-1:
    while m<len(ds_bc)-1:
        dT_ds.append((unit_vector_T_bc[l-1] - unit_vector_T_bc[l])/ds_bc[m])
        l += 1
        m += 1

#produces the modulus of dT/ds
for i in dT_ds:
    mod_dT_ds.append(modulus(i))

mod_dT_ds_bc=[]
for i in dT_ds_bc:
    mod_dT_ds_bc.append(modulus(i))

#determines average curvature
avg_mod_dT_ds_bc=sum(mod_dT_ds_bc[20:-20])/len(mod_dT_ds_bc[20:-200])
print 'average curvature is ',avg_mod_dT_ds_bc

#calculates the variance in the average
sq_avg_K_bc=[]
var_avg_K_bc=[]

a=0
while a<len(mod_dT_ds_bc[20:-20])-1:
    sq_avg_K_bc.append((mod_dT_ds_bc[20:-20][a]-avg_mod_dT_ds_bc)**2)
    a += 1

var_avg_K_bc=sum(sq_avg_K_bc)/len(sq_avg_K_bc)

std_avg_K_bc=math.sqrt(var_avg_K_bc)
print 'variance and std respectively are',var_avg_K_bc,std_avg_K_bc

######## 4. Calculating Torsion ########

#empty list for N, the normal unit vector
unit_vector_N=[]

#produces N using N=(dT/ds)/|dT/ds| (page 38 lab book)
q=0; r=0
while q<len(dT_ds)-1:
    while r<len(mod_dT_ds)-1:
        unit_vector_N.append(dT_ds[q]/mod_dT_ds[r])
        q += 1
        r += 1

unit_vector_N_bc=[]
q=0; r=0
while q<len(dT_ds_bc)-1:
    while r<len(mod_dT_ds_bc)-1:
        unit_vector_N_bc.append(dT_ds_bc[q]/mod_dT_ds_bc[r])
        q += 1
        r += 1

#empty list for K, the curvature-kappa
K=[]

#produces list of curvature values at each point using N and dT/ds
#where curvature K=(dT/ds)/N
#used to compare with other method of determining K, K=|dT/ds|
u=0;w=0
while w<len(unit_vector_N)-1:
    while u<len(dT_ds)-1:
        K.append(dT_ds[u]/unit_vector_N[w])
        u += 1
        w += 1

#####as all values in K vector are the same, turns xyz from list to array
#generates transpose of this array
array_K=np.array(K); transposed_K=array_K.T

curvature=transposed_K[0]


######## Torsion ########
#produces list of binormal unit vectors as B=NxT
#the crossproduct of N and T, shown on page 39

array_N=np.array(unit_vector_N)
array_T=np.array(unit_vector_T)

#boxcar average
curv_avg=[]

a=0
while a<len(curvature):
    curv_avg.append(sum(curvature[a-6:a+6])/12)
    a += 1

#empty list for B, the binormal unit vector
unit_vector_B=[]

o=0;p=0
while o<len(unit_vector_N):
    while p<len(unit_vector_T)-1:
        unit_vector_B.append(np.cross(unit_vector_T[o],unit_vector_N[p]))
        o += 1
        p += 1


unit_vector_B_bc=[]

o=0;p=0
while o<len(unit_vector_N_bc):
    while p<len(unit_vector_T_bc)-1:
        unit_vector_B_bc.append(np.cross(unit_vector_T_bc[o],unit_vector_N_bc[p]))
        o += 1
        p += 1

#empty list for differential of B wrt the distance between the points
dB_ds=[]

#produces list of dB/ds values in order to calculate torsion
a=0;b=0
while a<len(unit_vector_B)-1:
    while b<len(ds)-2:
        dB_ds.append(unit_vector_B[a]/ds[b])
        a += 1
        b += 1

dB_ds2=[]
#produces dB/ds from B and ds
l=0;m=0
while l<len(unit_vector_B):
    while m<len(ds)-2:
        dB_ds2.append((unit_vector_B[l-1] - unit_vector_B[l])/ds[m])
        l += 1
        m += 1

dB_ds_bc=[]
#produces dB/ds from B and ds
l=0;m=0
while l<len(unit_vector_B_bc):
    while m<len(ds_bc)-2:
        dB_ds_bc.append((unit_vector_B_bc[l-1] - unit_vector_B_bc[l])/ds_bc[m])
        l += 1
        m += 1

mod_dB_ds_bc=[]
a=0
while a<len(dB_ds_bc)-1:
    mod_dB_ds_bc.append(modulus(dB_ds_bc[a]))
    a += 1

mod_dB_ds2=[]
a=0
while a<len(dB_ds2)-1:
    mod_dB_ds2.append(modulus(dB_ds2[a]))
    a += 1

torsion_avg=[]
a=0
while a<len(dB_ds2):
    torsion_avg.append(sum(dB_ds2[a-2:a+2])/4)
    a += 1

#empty list for torsion
Torsion=[]

#calculates torsion using -|dB/ds|/|N| (page 39)
c=0;e=0
while c<len(dB_ds)-1:
    while e<len(unit_vector_N)-1:
        Torsion.append((modulus(dB_ds[c])/modulus(unit_vector_N[e]))*(-1))
        c += 1
        e += 1

Torsion2=[]
c=0;e=0
while c<len(dB_ds2)-1:
    while e<len(unit_vector_N)-1:
        Torsion2.append((modulus(dB_ds2[c])/modulus(unit_vector_N[e]))*(-1))
        c += 1
        e += 1

mod_tors2=[]
a=0
while a<len(Torsion2)-1:
    mod_tors2.append(modulus(Torsion2[a]))
    a += 1

#determines average torsion
avg_mod_dB_ds_bc=sum(mod_dB_ds_bc[21:-20])/len(mod_dB_ds_bc[21:-20])

#calculates variance in average
sq_avg_T_bc=[]
var_avg_T_bc=[]

a=0
while a<len(mod_dB_ds_bc[21:-20])-1:
    sq_avg_T_bc.append((mod_dB_ds_bc[21:-20][a]-avg_mod_dB_ds_bc)**2)
    a += 1

var_avg_T_bc=sum(sq_avg_T_bc)/len(sq_avg_T_bc)

std_avg_T_bc=math.sqrt(var_avg_T_bc)
print 'variance and std respectively are',var_avg_T_bc,std_avg_T_bc


######## 5. Plotting Trajectory ########

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(data[1], data[2], data[3])

ax.set_title("Trajectory")
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')


plt.show()

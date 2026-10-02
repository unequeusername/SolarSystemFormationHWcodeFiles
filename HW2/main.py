#Homework 2 Origins of Solar Systems


import numpy as np
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt


#Problem 1

#Problem 2

def Area(R1,R2):
	return np.pi * abs( (R2**2)-(R1**2))

R = np.array([5.9,9.4,15,26.4])*5.32e10*u.cm


R1 = np.array([4.67, 7.45, 11.87,19.9])*5.32e10*u.cm
R2 = np.array([7.45,11.87,19.9,35.02])*5.32e10*u.cm

A = Area(R1,R2)

#A = np.array([A(R2[i],R1[i]) for i in range(len(R))])

M = np.array([8.93,4.80,14.8,10.8])*1e25*u.g

Maug = M/np.array([0.005,0.005,0.010,0.010])

Sigma =Maug/A

#print("Area:",A,"\n","Sigma:",Sigma)

x = np.log10(R.to(u.km).value )
y = np.log10(Sigma.to(u.g*u.cm**-2).value)

w = 1/(0.1**2)
W =w*np.identity(len(x))

X = np.transpose(np.array([np.ones(len(x)),x]))
XT =np.transpose(X)
Y = np.transpose(y)
YT =y
#print(X,"\n",Y)

#B = np.matmul(np.linalg.inv(np.matmul(np.transpose(X),X)),np.matmul(np.transpose(X),Y))


covB = np.linalg.inv(XT@W@X)
B= covB@XT@W@Y
#print(B)
p= B[1]
dp =covB[1,1]**0.5
#print("p=",p,"+-",dp)

Tab2 =np.transpose(np.array([x,Maug,A,y,Sigma]))

#print(Tab2)
#plt.plot(x,y, marker=".")
#plt.plot(x,B[0] +B[1]*x)
#plt.xlabel("log(R)")
#plt.ylabel("log($\Sigma$)")
#plt.show()


#Problem 3

for i in [const.M_sun,const.M_earth,const.G, u.AU]:
	print(i.cgs)

#Solar System Formation Homework 1 
#Liam Royle-Grimes
#09/15/2026


#problem 1

#from periodictable import H,He,O,C,Ne,N,Mg,Si,Fe,S,Ar,Al,Ca,Na,Ni,Cr,Mn,P,Cl,K,Ti,Co
import periodictable.core
import numpy as np
import matplotlib.pyplot as plt


def init(table, reload =False):
	if "abund" not in table.properties or reload:
		table.properties.append("abund")
		periodictable.core.Element.abund = 0
		for name,a in Elem.items():
			el = table.symbol(name)
			el.abund = a
	if "Tc" not in table.properties or reload:
		table.properties.append("Tc")
		periodictable.core.Element.Tc = 0
		for name,t in Tc.items():
			el = table.symbol(name)
			el.Tc = t
	if "Onumber" not in table.properties or reload:
		table.properties.append("Onumber")
		periodictable.core.Element.Tc = 0
		for name,t in oxy.items():
			el = table.symbol(name)
			el.Onumber = t


Elem=dict(H = 2.431e10,
He = 2.343e9,
O = 1.413e7,
C = 7.079e6,
Ne = 2.148e6,
N = 1.950e6,
Mg = 1.020e6,
Si = 1.000e6,
Fe = 8.380e5,
S = 4.449e5,
Ar = 1.025e5,
Al = 8.410e4,
Ca = 6.287e4,
Na = 5.751e4,
Ni = 4.780e4,
Cr = 1.29e4,
Mn = 9.17e3,
P = 8.37e3,
Cl = 5.24e3,
K = 3.70e3,
Ti = 2.42e3,
Co = 2.32e3)

Tc = dict(H=182,He =3,Li=1142,Be=1452,B=908,C=78,N=131,O=182,F=739,Ne=9.3,Na=958,Mg=1397,Al=1677,Si=1529,P=1248,S=704,Cl=954,Ar=48,K=1006,Ti=1593,Ca=1659,Cr=1296,Mn=1158,Fe=1357,Co=1352,Ni=1353)

oxy = dict(H=0.5,Li=0.5,Be=1,B=1.5,O=0,F=0.5,Na=0.5,Mg=1,Al=1.5,Si=2,Cl=0.5,K=0.5,Ca=1,Sc=1.5,Ti=2,V=0,Br=0.5,Rb=0.5,Sr=1,Y=1.5,Zr=2,Nb=0,I=0.5,Cs=0.5,Ba=1,Lu=1.5,Hf=2,Ta=0)

init(periodictable.core.default_table())
from periodictable import *

Elem ={H,He,O,C,Ne,N,Mg,Si,Fe,S,Ar,Al,Ca,Na,Ni,Cr,Mn,P,Cl,K,Ti,Co}
atmophile = {H,He,C,N,Ne,Ar,Kr,Xe}
calcophile = {S,Cu,Zn,Ga,Ge,As,Se,Ag,Cd,In,Sn,Sb,Te,Hg,Ti,Pb,Bi}
siderophile = {P,Mn,Fe,Co,Ni,Mo,Ru,Rh,Pd,W,Re,Os,Ir,Pt,Au,Cr}
lithophile = {Li,Be,B,O,F,Na,Mg,Al,Si,Cl,K,Ca,Sc,Ti,V,Br,Rb,Sr,Y,Zr,Nb,I,Cs,Ba,Lu,Hf,Ta}
volatiles=[i for i in Elem if i.Tc<(300)]


print("volatiles")
for i in Elem:
	if i in volatiles:
		print(i,i.Tc)
print("siderophile elements")

siderophile = [i for i in siderophile if i in Elem and i not in volatiles]

for i in siderophile:
	print(i,i.Tc)

print("lithophile elements")

lithophile = [i for i in lithophile if i in Elem and i not in volatiles]
for i in lithophile:
	print(i,i.Tc)

#lithophille condensation
for i in lithophile:
	print(i,i.Onumber, i.Onumber*i.abund)

solidoxused = np.sum([i.Onumber*i.abund for i in lithophile])

totalabund = np.sum([i.abund for i in Elem])

gassoxused = O.abund - solidoxused

totalmass = np.sum([i.mass*i.abund for i in Elem])
silicatemass = solidoxused * O.mass + np.sum([i.abund*i.mass for i in lithophile])
metalmass = np.sum([i.mass*i.abund for i in siderophile])
gasmass = totalmass -silicatemass - metalmass
solidmass = silicatemass+metalmass

print("gas ratio",gasmass/totalmass)
print("metal ratio",metalmass/totalmass)
print("lith ratio",silicatemass/totalmass)

totalCaTiO = Ti.abund
massCaTiO = (Ca.mass+Ti.mass+3*O.mass)*totalCaTiO

totalCaAlSi = (Ca.abund - totalCaTiO)/2
massCaAlSi = (2*Ca.mass + 2*Al.mass + Si.mass +7*O.mass)

totalAlO =  (Al.abund - 2*totalCaAlSi)/2
massCaAlO = (2*Al.mass +3*O.mass)*totalAlO

CAImass = massCaTiO+massCaAlSi+massCaAlO

print("CAIs", CAImass,solidmass)
print(CAImass/solidmass)

#problem 2
print("P2")
#a
x1 = 13.1
x2=61.37
y1=0.139097
y2=0.139794

m = (y1-y2)/(x1-x2)
b = y1-m*x1

print(m,b)
#b

x=[142.7,101.7,61]
dx=[18.0,17.4,5.2]
y=[0.980,0.637,0.457]
dy=[0.173,0.173,0.068]
dy1sig =[i/2 for i in dy]

S0 = np.sum([dy[i]**-2 for i in range(3)])
Sx = np.sum([x[i]/dy[i]**2 for i in range(3)])
Sy = np.sum([y[i]/dy[i]**2 for i in range(3)])
Sxx = np.sum([x[i]**2/dy[i]**2 for i in range(3)])
Sxy = np.sum([y[i]*x[i]/dy[i]**2 for i in range(3)])
Syy = np.sum([y[i]**2/dy[i]**2 for i in range(3)])
Delta = S0*Sxx-Sx*Sx
b = (Sy*Sxx-Sx*Sxy)/Delta
db= (Sxx/Delta)**0.5
m = (S0*Sxy-Sx*Sy)/Delta
dm = (S0/Delta)**1/2
MSWD = np.sum([((y[i] -m*x[i]-b)**2/(dy1sig[i])**2) for i in range(3)])
critval= 1+2*(2/(3-2))**0.5
print("section b")
print(S0,Sx,Sy,Sxx,Sxy,Syy,Delta,b,db,m,dm,MSWD,critval)
plt.plot(x,y)
plt.plot(x,[m*x[i]+b for i in range(3)])
plt.show()


#problem 3
V= (np.pi*(3.27**2- 2.05**2)*0.5)*(1.46e11)**3
N=1e6
n = N/V
d = n**-(1/3)

sig = np.pi*(0.5e3)**2
v= 5e3

R = n*sig*v

t = 1/R
print("P3")
print(V)
print(n)
print(d)
print("sigma",sig)
print("colision time")
print(R)
print(t)
print(t/3.154e7)
print(t/3.154e7/1e9)
print(R*4.5e9*3.15e7)
print(N*R*4.5e9*3.15e7/2)

















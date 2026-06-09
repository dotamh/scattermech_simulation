# -*- coding: utf-8 -*-
"""
filmcurve_simulator
"""

import matplotlib.pyplot as plt
from pySCATMECH.local import *
import pandas as pd
from pySCATMECH.fresnel import *
from pySCATMECH.integrate import *
import numpy as np

pi = np.pi
deg = pi/180

#User Input1：薄膜材料折射率
# User Input film type, default film index n + j*k as follwoing, user as the option to change those number in diag. window
Si = OpticalFunction(1.850864+4.431444j)    # n + k*j, 根据用户输入的n、k生成薄膜材料的折射率，此处为材料Si
Oxide = OpticalFunction(1.510800+0.000203j)
Nit = OpticalFunction(2.22161+0.00086281j)
Air = OpticalFunction(1)

#CO = OpticalFunction(1.799+1.767j)
#SiN = OpticalFunction(2.22161+0.00086281j)
#SiOCH = OpticalFunction(1.462 + 0.004j)
#SiO = OpticalFunction(1.46)
#Ti = OpticalFunction(3+3.17j)
#SiN = OpticalFunction(1.98+0.04j)
#AlO = OpticalFunction(1.62)

#User Input2：薄膜材料厚度，注意此处单位是um，1um = 10000A
# User Input film thickness, use only one single film layer. 
#stack = {None:'No_StackModel'} # for Bare Si
# 此处需要UI输入一个FilmStack的表格，包括每一次薄膜的（n、k）和Thickness；由靠近Substrate的film作为第一层

# film1 = Film(Oxide,thickness = 500/10000) #closer to substrate

# film2 = Film(Nit,thickness = 10000/10000)   #以此类推

# 可以一直叠加
# stack = FilmStack([film1,film2])
stack = {None:'No_StackModel'}



# User Input3：入射角，Oblique则选择75度，Normal则选择0度；
thetai = 75*deg # Oblique thetai = 75deg, Normal thetai = 0 deg

# User Input4：入射偏振 IncPol
incpol = Polarization('p') # Oblique case: p or s or c; Normal case: c only


#Config Input1：定义Silica小球和PSL小球的折射率
Silica = OpticalFunction(1.510800+0.000203j) #Silica Sphere
#Silica = OpticalFunction(1.4994);
PSL = OpticalFunction(1.743846) # PSL Sphere #1.743846, #2.0352

# 定义模型与模型参数部分
spherecoat = {None:'No_StackModel'}  #hard coded
parameters = {'lambda' : 0.266, #hard coded
              'substrate' : Si, #Si Substrate, hard coded
              'type' : 0,
              'density': 1 ,      #hard coded
              'sphere' : PSL, #UserInput5：小球材料类型
              'radius' : 0.05, # default setting, hard coded
              'spherecoat' : spherecoat,  #hard coded
              'stack' : stack,            #hard coded
              'delta' : 0,  #hard coded
              'lmax' : 0,   #hard coded
              'order' : -1,  #hard coded
              'Norm_Inc_Approx' : 0,  #hard coded
              'improve' : 3}  #hard coded

model = Local_BRDF_Model("Bobbert_Vlieger_BRDF_Model",parameters)  #hard coded
print(model)


#预定义机台的P1-P5和PU，for ch1
p5 = Polarization(angle = 0*deg) # 0deg, s
p5 = StokesVector(np.array(p5)*0.5) 

p1 = Polarization(angle = 22.5*deg)
p1 = StokesVector(np.array(p1)*0.5) #normalize to perserve total energy

p2 = Polarization(angle = 45*deg)
p2 = StokesVector(np.array(p2)*0.5) 

p3 = Polarization(angle = 67.5*deg)
p3 = StokesVector(np.array(p3)*0.5) 

p4 = Polarization(angle = 90*deg)  #90deg, p
p4 = StokesVector(np.array(p4)*0.5) 

p1r = Polarization(angle = -67.5*deg)
p1r = StokesVector(np.array(p1r)*0.5) 

p2r = Polarization(angle = -45*deg)
p2r = StokesVector(np.array(p2r)*0.5) 

p3r = Polarization(angle = -22.5*deg)
p3r = StokesVector(np.array(p3r)*0.5) 

p4r = Polarization(angle = 0*deg)
p4r = StokesVector(np.array(p4r)*0.5) 

p5r = Polarization(angle = 90*deg)
p5r = StokesVector(np.array(p5r)*0.5)  

pu = Polarization('u')


# 根据用户输入定义收集相关部分
# P-FULL-U - IDEAL
# full = CircularCone(theta=0*deg,phi = 0*deg,alpha = 75*deg,sensitivity=p2)
# halfmask = ProjectedPolygon([(-1,0),
#                              (-1,-1),
#                              (1,-1),
#                              (1,0)],sensitivity=p2)
#                        #  (-0.5,0.5),
#                        #  (0,0.2)])
# narrow = CircularCone(theta=0*deg,phi = 0*deg,alpha = 28*deg,sensitivity=p2)
# detector1 = full & halfmask & ~narrow
# ch1_integrator = Integrator(1*deg, detector1,type = 1)
# ch1_integrator.PlotSamplingPoints() # drawing


full = CircularCone(theta=0*deg,phi = 0*deg,alpha = 75*deg,sensitivity=pu)

halfmask = ProjectedPolygon([(-1,0),
                             (-1,-1),
                             (1,-1),
                             (1,0)],sensitivity=p2)

narrow = CircularCone(theta=0*deg,phi = 0*deg,alpha = 28*deg,sensitivity=pu)

# sector = ProjectedPolygon([(-np.cos(67.5*deg),np.sin(67.5*deg)),(-np.cos(56.25*deg),np.sin(56.25*deg)), # 定义一个扇形区域，用多边形去逼近，同时再对坐标值做一个判断，去除没有被多边形完全删除的点
#                          (-np.cos(45*deg),np.sin(45*deg)), (-np.cos(33.75*deg),np.sin(33.75*deg)),
#                          (-np.cos(22.5*deg),np.sin(22.5*deg)), (-np.cos(11.25*deg),np.sin(11.25*deg)),
#                          (-1,0), (-np.cos(11.25*deg),-np.sin(11.25*deg)), (-np.cos(22.5*deg),-np.sin(22.5*deg)),
#                          (-np.cos(33.75*deg),-np.sin(33.75*deg)), (-np.cos(45*deg),-np.sin(45*deg)),
#                          (-np.cos(56.25*deg),-np.sin(56.25*deg)), (-np.cos(67.5*deg),-np.sin(67.5*deg)), (0,0)], sensitivity=pu)

sector = ProjectedPolygon([(-1,0), (-np.cos(11.25*deg),-np.sin(11.25*deg)), (-np.cos(22.5*deg),-np.sin(22.5*deg)),
                         (-np.cos(33.75*deg),-np.sin(33.75*deg)), (-np.cos(45*deg),-np.sin(45*deg)),
                         (-np.cos(56.25*deg),-np.sin(56.25*deg)), (-np.cos(67.5*deg),-np.sin(67.5*deg)), (0,0)], sensitivity=pu)
# sector= ProjectedPolygon([(0, 0),
#                           (0,-1),
#                           (1, 0)])

slot = RectangularCone(theta=0, phi=0, alpha=(75*deg, 13.7659*deg), sensitivity=pu)
detector1 = full & halfmask & ~narrow & ~sector & ~slot
ch1_integrator = Integrator(1*deg, detector1,type = 3)
# detector1 = halfmask
# ch1_integrator = Integrator(1*deg, detector1,type =  3)
ch1_integrator.PlotSamplingPoints()


# P-FULL-U- Cut
block1 = ProjectedPolygon([(0.95,0.1),
                             (0.95,-0.1),
                             (1,-0.1),
                             (1,0.1)],sensitivity=p2)


block2 = ProjectedPolygon([(-0.95,0.1),
                             (-0.95,-0.1),
                             (-1,-0.1),
                             (-1,0.1)],sensitivity=p2)

block3 = ProjectedPolygon([(0.1,0.95),
                             (-0.1,0.95),
                             (-0.1,1),
                             (0.1,1)],sensitivity=p2)

block4 = ProjectedPolygon([(0.1,-0.95),
                             (-0.1,-0.95),
                             (-0.1,-1),
                             (0.1,-1)],sensitivity=p2)
block4 = ProjectedPolygon([(0.1,-0.95),
                             (-0.1,-0.95),
                             (-0.1,-1),
                             (0.1,-1)],sensitivity=p2)

detector2 = detector1 & ~block1  & ~block2 & ~block4
ch2_integrator = Integrator(1*deg, detector2,type = 1)
ch2_integrator.PlotSamplingPoints()

# 举例Ch#3,P#定义只有P、S、U
cone5 = CircularCone(theta=0*deg,phi = 0*deg,alpha = 25*deg,sensitivity=Polarization('u'))
detector3 =cone5 &  ~block3 & ~block2
ch3_integrator = Integrator(1*deg, detector3,type = 1)
ch3_integrator.PlotSamplingPoints()


# 根据用户输入的粒径计算ISC数值
#diameters = np.linspace(0.02,0.8,0.2); #
diameters = np.array([0.02,0.026,0.032,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.3,0.5,0.8])

# 计算对应通道的ISC数值
resCh1 = []
resCh2 = []
resCh3 = []

for d in diameters:
    model.setParameters(radius=d/2)
    resCh1.append(ch1_integrator.CrossSection(model, thetai, incpol) *2)
    resCh2.append(ch2_integrator.CrossSection(model, thetai, incpol) *2)
    resCh3.append(ch3_integrator.CrossSection(model, thetai, incpol) *2)

plt.figure()
plt.semilogy(diameters,resCh1,label='ch1-p#')
plt.semilogy(diameters,resCh2,label='ch2-p#r') 
plt.semilogy(diameters,resCh3,label='ch3-pu')
plt.xlabel("Diameter / $\mathrm{\mu m}$")
plt.ylabel("Cross section / $\mathrm{\mu m}^2$")
plt.title("Cross Section versus Diameter")
plt.legend()
plt.show()

a = np.vstack((diameters,resCh1,resCh2,resCh3))


# 储存对应通道的ISC数值
np.savetxt("fc.csv", a, delimiter = ",")

list=[]
list.append(diameters)
list.append(resCh1)
list.append(resCh2)
list.append(resCh3)
print(list)
column=['diameters','resCh1','resCh2','resCh3'] #列表头名称
test=pd.DataFrame(index=column,data=list)#将数据放进表格
test.to_csv('D:\Autofilm_Save_Data\Autofilm.csv') #数据存入csv,存储位置及文件名称

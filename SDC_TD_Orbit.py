import csv
import math
import numpy as np

β_min = 0 #β角の最小値
μ = 3.98603*(10**(5)) #ケプラー定数
Req = 6378.1 #地球の赤道半径
h = 200 # 軌道高度
R = Req + h # 半直弦

entry = 104.164
exit = 255.836

θ_en = math.radians(entry)
θ_ex = math.radians(exit)

β = β_min

T = (4*math.pi*math.pi*(R**3)/μ)**(1/2) #周期計算
n = 100 #timesの精度
dt = T/n

data = [] #times,x,y,z,roll,pitch,yow

counter_en = 0
counter_ex = 0

for i in np.arange(0,n):
  θ = 2*dt*i*math.pi/T

  if θ < θ_en:
    times = i*dt
    θ = 2*dt*i*math.pi/T
    x = R*math.cos(θ)*math.cos(math.radians(β))
    y = R*math.sin(θ)
    z = R*math.cos(θ)*math.sin(math.radians(β))
    data.append([times,x,y,z,0,0,0])

  if θ_ex > θ > θ_en:
    if counter_en == 0:
      θ_ = θ_en
      times = θ_*T/(2*math.pi)
      x = R*math.cos(θ_)*math.cos(math.radians(β))
      y = R*math.sin(θ_)
      z = R*math.cos(θ_)*math.sin(math.radians(β))
      counter_en =counter_en+1
      data.append([times,x,y,z,0,0,0])
    
    times = i*dt
    θ = 2*dt*i*math.pi/T
    x = R*math.cos(θ)*math.cos(math.radians(β))
    y = R*math.sin(θ)
    z = R*math.cos(θ)*math.sin(math.radians(β))
    ψ = math.degrees(math.pi/2+math.asin(y/R))

    data.append([times,x,y,z,0,0,-ψ])

  if θ > θ_ex:
    if counter_ex == 0:
      θ_ = θ_ex
      times = θ_*T/(2*math.pi)
      x = R*math.cos(θ_)*math.cos(math.radians(β))
      y = R*math.sin(θ_)
      z = R*math.cos(θ_)*math.sin(math.radians(β))
      counter_ex =counter_ex+1
      ψ = math.degrees(math.pi/2+math.asin(y/R))
      data.append([times,x,y,z,0,0,-ψ])
    
    times = i*dt
    θ = 2*dt*i*math.pi/T
    x = R*math.cos(θ)*math.cos(math.radians(β))
    y = R*math.sin(θ)
    z = R*math.cos(θ)*math.sin(math.radians(β))

    data.append([times,x,y,z,0,0,0])


with open('sdc_orbit','w') as f:
  writer = csv.writer(f)
  writer.writerows(data)

f.close



import numpy as np
from environment import LEOEnvironment
from satellite import Satellite
import time

import matplotlib.pyplot as plt


# 開始
start_time = time.process_time()

# 座標系について
# 宇宙機座標系は変数の後ろにBをつける
# Zを地球方向　Xを進行方向　それに合わせてY軸を表す．
# 
# 全体座標系(地心赤道面をXY平面とし太陽を+X軸方向の先に固定した座標系)

# LEO
en = LEOEnvironment(
    vec_sun_W = np.array([-1,0,0]),
    beta = 0,
    h = 400
)

# 衛星
sat = Satellite(
    filepath = "box.stl", 
    vec_sun_B = en.vec_sun_W
)

pcsas = []
for theta in np.linspace(0, 360, 10):
    print(theta)
    theta = np.radians(theta)
    # 影にある時
    if en.theta_shoku <= theta and theta <= 2*np.pi - en.theta_shoku:
        pcsas.append(0)
    else:
        # 衛星の位置
        vec_sat_r_W = en.calc_sat_position(theta)
        # ベクトルを全体座標系から宇宙機座標系にする
        vec_sun_B = sat.W2B(
                vec_sat_r_W,
                theta,
                en.beta,
                e = None # 向きたい方向のベクトルを与える
            ) @ en.vec_sun_W
        print(vec_sun_B)
        sat.pcsa_calc(vec_sun_B)
        pcsas.append(sat.pcsa)
print(pcsas)

# 終了
end_time = time.process_time()

# 経過時間を出力(秒)
elapsed_time = end_time - start_time
print(elapsed_time)

fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111) 
ax.set_xlabel('theta')
ax.set_ylabel('pcsa')

ax.grid()

ax.autoscale()
ax.plot(np.linspace(0, 360, 10),pcsas, color = "red")
plt.show()



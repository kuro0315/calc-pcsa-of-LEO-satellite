import numpy as np
from environment import LEOEnvironment
from satellite import Satellite

# 座標系
# 宇宙機座標系は変数の後ろにBをつける
# Zを地球方向　Xを進行方向　それに合わせてY軸を表す．
# 
# 全体座標系(地心赤道面をXY平面とし太陽をX軸に固定した座標系)

# LEO
en = LEOEnvironment(
    vec_sun_W = np.array([-1,0,0])
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
        continue
    else:
        # 衛星の位置
        vec_sat_r_W = en.calc_sat_position(theta)
        # ベクトルを全体座標系から宇宙機座標系にする
        vec_sun_B = sat.W2B(
                vec_sat_r_W,
                theta,
                en.beta
            ) @ en.vec_sun_W
        print(vec_sun_B)
        sat.pcsa_calc(vec_sun_B)
        pcsas.append(sat.pcsa)
print(pcsas)
    
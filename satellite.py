import numpy as np
from pcsa_calculator import PcsaCalculator

class Satellite(PcsaCalculator):
    def __init__(self, filepath, vec_sun_B):    
        super().__init__(filepath, vec_sun_B)    
    
    # 全体の座標系から宇宙機座標系へ
    def W2B(self,vec_sat_r_W, theta,beta):
        #原点
        origin_B = vec_sat_r_W
        # x軸方向の単位ベクトル rをθで微分したもの
        ex_B = np.array(
                [ -np.sin(theta)*np.cos(beta), np.cos(theta) ,-np.sin(theta)*np.sin(beta) ]
            )
        ez_B = -vec_sat_r_W / np.linalg.norm(vec_sat_r_W)
        ey_B= np.cross(ex_B, ez_B)
        
        #絶対座標から宇宙機座標に変換する行列
        R_GtoL = np.identity(4)
        R_GtoL[:3, :3] = np.array([ex_B, ey_B, ez_B])
        R_GtoL[:3, 3] = - np.dot(np.array([ex_B, ey_B, ez_B]), origin_B)
        
        return R_GtoL[:3, :3]
    
    def calc_sat_position():
        pass
    

    
    
import numpy as np
from pcsa_calculator import PcsaCalculator

class Satellite(PcsaCalculator):
    def __init__(self, filepath, vec_sun_B):    
        super().__init__(filepath, vec_sun_B)    
    
    # 全体の座標系から宇宙機座標系へ
    def W2B(self,vec_sat_r_W, theta,beta,e=None):
        
        if e == None:
            # x軸方向の単位ベクトル rをθで微分したもの
            ex_B = np.array(
                    [ -np.sin(theta)*np.cos(beta), np.cos(theta) ,-np.sin(theta)*np.sin(beta) ]
                )
            ez_B = -vec_sat_r_W / np.linalg.norm(vec_sat_r_W)
            ey_B= np.cross(ex_B, ez_B)
        else :
            ex_B, ey_B, ez_B = e
        
        #絶対座標から宇宙機座標に変換する行列
        R_GtoL = np.array([ex_B, ey_B, ez_B])
        
        return R_GtoL
    
    def calc_sat_position():
        pass
    

    
    
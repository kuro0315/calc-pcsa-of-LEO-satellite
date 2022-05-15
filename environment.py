from skyfield.api import load
from skyfield import api
import numpy as np
from datetime import datetime, timedelta
from pytz import timezone

class LEOEnvironment():
    def __init__(
            self, 
            vec_sun_W, 
            beta=30, 
            i=51.64, 
            re=6378.1 , 
            h=400
        ):
        # 太陽面軌道角
        self.beta = np.radians(beta)
        # 高度[km]
        self.h = h
        # 軌道傾斜角
        self.i = np.radians(i)
        # 地球の赤道半径
        self.re = re
        # 太陽方向のベクトル
        self.vec_sun_W = vec_sun_W
        # 半直弦
        self.R = self.re + self.h
        # 日影に入る角度
        self.theta_shoku = self.calc_shoku()
        
        # 地球中心から太陽方向と地球中心から衛星方向の角度
        self.theta = 0
        # 衛星の位置計算
        self.vec_sat_r_W = self.calc_sat_position(self.theta)

    def calc_sat_position(self,theta):
            self.theta = theta
            # 衛星の位置
            self.vec_sat_r_W = (self.re + self.h) * np.array(
                [ np.cos(theta)*np.cos(self.beta), np.sin(theta) ,np.cos(theta)*np.sin(self.beta) ]
            )
            return self.vec_sat_r_W
        
    def calc_shoku(self):
        theta_shoku = np.arcsin(
            (
                1.0/np.cos(self.beta)**2 * ( ( self.re / (self.re + self.h) )**2 - np.sin(self.beta)**2 )
            )**0.5
        )
        theta_shoku = np.pi - theta_shoku
        return theta_shoku
    
    
    
    
    
    
        
    
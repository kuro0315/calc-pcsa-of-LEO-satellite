# -*- coding: utf-8 -*-
from pyexpat import model
import numpy as np
import pyvista as pv
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import patches
import itertools
from pprint import pprint
from utils import *
import random
import copy

##
# @file pcsa_calculator.py
# @author Shuntarou Kuroiwa, Nishio Keita
# @date 2022/03/30
# @brief 展開構造をもつ衛星の投影面積を解析的に求めます
# @details 詳細な説明

##
# @class PcsaCalculator
# @brief 投影面積を求めるライブラリ
# @details cadファイルと投影面積を求めたい方向のベクトルから投影面積を計算するライブラリです． CADファイルのパスと投影面積を求めたい方向のベクトルを与えてください．
class PcsaCalculator():
    def __init__(self, filepath, vec_k):
        ## ファイルのパス
        self.filepath = filepath
        ## 衛星のCADモデル
        self.model = pv.read(filepath)
        ## 頂点
        self.vertices = self.model.points
        ## 面を構成する頂点 (1つの面(＝3つの頂点セット)でreshapeしている。 stlやめた時変更する必要あり)
        self.index_polygons = self.model.faces.reshape(-1, 4)[:, 1:]
        ## 面を構成する頂点 (1つの面(＝3つの頂点セット)でreshapeしている。 stlやめた時変更する必要あり)
        self.polygons = []
        for face in self.model.faces.reshape(-1, 4)[:, 1:]:
            buffer = []
            for index_vertice in face:
                buffer.append(self.vertices[index_vertice])
            self.polygons.append(buffer)
        # pcsaの計算!
        self.pcsa_calc(vec_k)
    ##
    # @brief 投影面積の計算
    # @param vec_k 宇宙機座標系での求める投影面積方向に沿ったベクトル
    # @return None 戻り値なし
    def pcsa_calc(self, vec_k):
        ## ポリゴンを構成する頂点が属する面のパラメータの配列
        self.planes = []
        for polygon in self.polygons:
            # 面を構成する点の座標
            vec_r1 = polygon[0]
            vec_r2 = polygon[1]
            vec_r3 = polygon[2]
            # 平面のパラメータを計算
            self.planes.append(self.calc_plane_param(vec_r1,vec_r2,vec_r3))
        
        ## Kベクトルに対してある垂直な面に投影した頂点
        self.vertices_proj = []
        for vertice in copy.deepcopy(self.vertices):
            self.vertices_proj.append(
                self.calc_projection_of_vertex(
                        vertice, vec_k
                    )
                )
        
        ## W座標系での投影した頂点
        self.vertices_proj_w = []
        # 投影した面を構成する頂点を全てkベクトルを基準とした座標系での頂点に座標変換する．
        for vertice in self.vertices_proj:
            vertice =  self.b2w_matrix(vec_k) @ vertice
            self.vertices_proj_w.append([
                vertice[0,0], vertice[0,1], vertice[0,2]
            ])
        self.vertices_proj_w = np.array(self.vertices_proj_w)
        
        ## W座標系での投影したポリゴン
        polygons_proj_w = []
        for i in self.index_polygons:
            buffer = []
            for j in i:
                buffer.append(self.vertices_proj_w[j])
            polygons_proj_w.append(buffer)
        polygons_proj_w = np.array(polygons_proj_w)
        
        # 線になってしまったポリゴンを省く
        polygons_proj_w = self.remove_line_from_polygons_proj_w(polygons_proj_w)
        
        # 1回目は普通のポリゴンと重なってできたポリゴンの二つを計算する必要がある
        first_flag = 1
        # 重なってできたポリゴン用
        overlapping_polygons_sorted = []
        # overlappin_poygonsがあるか一回目の時以外は続ける
        count = -1
        while(overlapping_polygons_sorted != [] or first_flag):
            ## 重なってできたポリゴン
            overlapping_polygons = []
            ## 面の組み合わせ
            comb_polygons = list(itertools.combinations((polygons_proj_w), 2))
            # 交差点と内部点と一致している点を求めて重なっているポリゴンを求める
            for index_polygons, comb_polygon in enumerate(comb_polygons):
                ## 重なってできたポリゴンの頂点
                overlapping_vertices = []

                ## 交差点
                intersection_vertices = self.identify_intersection_points(comb_polygon)

                ## 同じ点
                same_vertices = self.identify_same_points(comb_polygon)

                ## 内部点
                interior_vertices = self.identify_interior_points(comb_polygon, []) # 1つ目のポリゴンの頂点が2つ目のポリゴンに含まれるかどうか計算
                interior_vertices = self.identify_interior_points(comb_polygon[::-1], interior_vertices) # 2つ目のポリゴンの頂点が1つ目のポリゴンに含まれるかどうか計算

                ## 重なってできたの頂点
                overlapping_vertices = intersection_vertices + same_vertices + interior_vertices

                if len(overlapping_vertices) >= 3:
                    if not search_polygons_4_polygons_array(overlapping_polygons, overlapping_vertices):
                        overlapping_polygons.append(overlapping_vertices)

            ## 重なってできたポリゴンの頂点をソートしたもの
            overlapping_polygons_sorted = self.sort_overlapping_polygon_vertices(overlapping_polygons)
            # 面積計算
            #　1回目は元のポリゴンと重なってできたポリゴン
            if first_flag == 1:
                G = self.calc_pcsa_from_polygons(polygons_proj_w)
                over_G = self.calc_pcsa_from_polygons(overlapping_polygons_sorted)
                ## 投影面積
                self.pcsa = G - over_G
                first_flag = 0

            #　2回目以降はさらに重なってできたポリゴン
            else: 
                # 重なってできたポリゴン同士に重なりがない時
                if overlapping_polygons_sorted == []:
                    pass #終了
                # 重なってできたポリゴン同士に重なりがある，つまり，余分に引いてる
                else:
                    over_G = self.calc_pcsa_from_polygons(overlapping_polygons_sorted)
                    self.pcsa += over_G
                    
            polygons_proj_w = overlapping_polygons_sorted[:]
            
            count += 1

    ##
    # @brief 宇宙機座標系からW座標系に変換する行列生成する関数
    # @param vec_k 宇宙機座標系での求める投影面積方向に沿ったベクトル
    # @return None 宇宙機座標系からW座標系に変換する行列生成する行列
    def b2w_matrix(self, vec_k):
        ## 球面座標系でいう θ
        delta = np.arcsin(vec_k[2] / np.linalg.norm(vec_k))
        if np.isclose(delta, np.radians(90)):
            delta = 90
        elif np.isclose(delta, -np.radians(90)):
            delta = -90
        buffer = vec_k[0] / (np.linalg.norm(vec_k) * np.cos(delta))
        
        if np.isclose(buffer,1.0):
            buffer = 1.0
        elif np.isclose(buffer,-1.0):
            buffer = -1.0
            
        ## 球面座標系でいう alpha
        alpha = np.arccos(buffer)
        if vec_k[1] < 0:
            alpha = 2*np.pi - alpha


        ## 宇宙機座標系からkベクトルを基準とする座標系に変換する行列
        brw = np.matrix((
            (np.cos(delta-np.pi/2) * np.cos(alpha), np.cos(delta-np.pi/2) * np.sin(alpha), np.sin(delta-np.pi/2)),
            (-np.sin(alpha), np.cos(alpha), 0),
            (-np.sin(delta-np.pi/2)*np.cos(alpha), -np.sin(delta-np.pi/2)*np.sin(alpha), np.cos(delta-np.pi/2))
        ))
        return brw

    ##
    # @brief 平面のパラメータを計算する
    # @param vec_r1,vec_r2,vec_r3 平面に含まれる3点
    # @return None 平面のパラメータ
    def calc_plane_param(self,vec_r1,vec_r2,vec_r3):
        x1,y1,z1 = vec_r1
        x2,y2,z2 = vec_r2
        x3,y3,z3 = vec_r3
        A = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1)
        B = (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1)
        C = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
        D = -(A*x1+B*y1+C*z1)
        return [A,B,C,D]

    ##
    # @brief 投影した点を計算
    # @param arg argの説明
    # @return None 戻り値なし
    def calc_projection_of_vertex(self, vec_r, vec_k):
        ## 見たい方向に沿ったベクトルの0の数
        n = len(vec_k) - np.count_nonzero(vec_k)
        
        if n == 3:
            raise ValueError("見たい方向のベクトルがおかしいよ")

        # この辺無理やりかいたすまぬ
        if n == 0:
            vec_r_proj = np.linalg.inv(
                    np.array([
                        [vec_k[1], -vec_k[0], 0],
                        [0, vec_k[2], -vec_k[1]],
                        [vec_k[0], vec_k[1], vec_k[2]]
                    ])
                ) @ np.array([
                    vec_k[1]*vec_r[0] - vec_k[0]*vec_r[1],
                    vec_k[2]*vec_r[1] - vec_k[1]*vec_r[2],
                    0
                ])
            return vec_r_proj

        elif n == 1:
            buffer = []
            for i in range(len(vec_k)):
                if vec_k[i] != 0:
                    buffer.append(i)

            buffer_vec_r_proj = np.linalg.inv(
                    np.array([
                        [vec_k[buffer[1]], -vec_k[buffer[0]]],
                        [vec_k[buffer[0]],  vec_k[buffer[1]]]
                    ])
                ) @ np.array([
                    vec_k[ buffer[1] ]*vec_r[ buffer[0] ] - vec_k[ buffer[0] ]*vec_r[ buffer[1] ], 
                    0
                ])
            
            vec_r_proj = vec_r[:]
            for index,num in enumerate(buffer):
                vec_r_proj[num] = buffer_vec_r_proj[index]
                
            return vec_r_proj

        else: #n==2
            vec_r_proj = np.array([0,0,0])
            for i in range(len(vec_k)):
                if vec_k[i] != 0.0:
                    continue
                vec_r_proj[i] = vec_r[i]
            return vec_r_proj
    
    ##
    # @brief 二組のポリゴンに交差点があるか計算する
    # @param comb_polygon 二組のポリゴン(W座標系)
    # @return intersection_vertices ポリゴンの交差点
    def identify_intersection_points(self, comb_polygon):
        intersection_vertices = []
        # 2種類のポリゴン
        polygon1 = comb_polygon[0]
        polygon2 = comb_polygon[1]
    
        for i, vertices1 in enumerate(polygon1):
            vec_r0 = np.array(polygon1[i])
            vec_r1 = np.array(polygon1[np.mod(i+1, len(polygon1))])
            for j, vertices2 in enumerate(polygon2):
                vec_r2 = np.array(polygon2[j])
                vec_r3 = np.array(polygon2[np.mod(j+1, len(polygon2))])
                intersection_vertice_t_s = self.calc_intersection_vertice(vec_r0,vec_r1,vec_r2,vec_r3)
                isIntersection_vertice = self.judge_intersect(intersection_vertice_t_s)
                # intersectionがあった
                if isIntersection_vertice:
                    t = intersection_vertice_t_s[0]
                    intersection_vertice = vec_r0 + t*(vec_r1 - vec_r0)
                    # もう追加している点ではないか調べる
                    if not search_array_4_array(intersection_vertices, intersection_vertice):
                        intersection_vertices.append(intersection_vertice)

        return intersection_vertices

    #
    # @brief 二組のポリゴンに同じ点があるか計算する deprecated
    # @param comb_polygon 二組のポリゴン(W座標系)
    # @return same_vertices ポリゴンの同じ点
    def identify_same_points(self, comb_polygon):
        same_vertices = []
        for vertice0 in comb_polygon[0]:
            for vertice1 in comb_polygon[1]:
                c1 = np.isclose( round(vertice0[0], 5), round(vertice1[0], 5) )
                c2 = np.isclose( round(vertice0[1], 5), round(vertice1[1], 5) )
                if  c1 and c2:
                    same_vertices.append(vertice0)
        return same_vertices

    ##
    # @brief 二組のポリゴンに内部の点があるか計算する
    # @param comb_polygon 二組のポリゴン(W座標系)
    # @return interior_vertices ポリゴンの内部点
    def identify_interior_points(self, comb_polygon, interior_vertices):
        for vertice in comb_polygon[0]:
            # ある点から対象の面の各点へのベクトルを求めている
            vec_list = []
            # 始点
            start_vertice = np.array(vertice)
            for vertice_end in comb_polygon[1]:
                # 終点
                end_vertice = np.array(vertice_end)
                vec_list.append(end_vertice - start_vertice)

            # 頂点が一致してしまった時0を省きたい
            vec_list_new = []
            for vec in vec_list:
                flag = False
                for v in vec:
                    v = round(v, 6)
                    if not np.isclose(v,0):
                        flag = True
                if flag:
                    vec_list_new.append(vec)

            # 始点からあるポリゴンの頂点に向けたベクトル全てのなす角の合計
            sum_theta = 0
            # 最初のベクトルが左回りか右回りか  左回りはTrue(+)
            vec_z = True
            for index, vec in enumerate(vec_list_new):
                i = np.inner(vec_list_new[index], vec_list_new[np.mod(index+1, len(vec_list_new))])
                n = np.linalg.norm(vec_list_new[index]) * np.linalg.norm( vec_list_new[np.mod(index+1, len(vec_list_new))] ) 
                c = i / n
                if np.isclose(c,1):
                    c = 1
                elif np.isclose(c,-1):
                    c = -1
                theta = np.rad2deg(np.arccos(c))
                if np.isclose(theta, 180):
                    theta = -theta
                # 最初の基準となるベクトルがどっち回りなのか計算する
                #vec_i_n = np.cross(vec_list_new[index], vec_list_new[index+1])
                sum_theta += theta

            sum_theta = round(sum_theta, 6)
            if np.isclose( np.abs(sum_theta), 0.0 ):
                pass
            elif np.isclose( np.abs(sum_theta), 360.0 ):
                if not search_array_4_array(interior_vertices, vertice):
                    interior_vertices.append(vertice)
                    ## print("INTERNAL(360)", vertice)
            elif np.isclose( np.abs(sum_theta), 180.0 ):
                pass
            else :
                pass

        return interior_vertices
    
    ##
    # @brief 重なってできたポリゴンの頂点をソートする
    # @param overlapping_polygons 重なってできたポリゴン
    # @return overlapping_polygons_sorted 重なってできたポリゴンの頂点をソートしたもの
    def sort_overlapping_polygon_vertices(self, overlapping_polygons):
        overlapping_polygons_sorted = []
        for polygon in overlapping_polygons:
        
            if polygon == []:
                continue

            # 重なっているポリゴンの中心を求める
            polygon_center = np.array([0.,0.,0.])
            for vertice in polygon:
                polygon_center[0] += vertice[0]
                polygon_center[1] += vertice[1]
                polygon_center[2] += vertice[2]
            polygon_center[0] = polygon_center[0] / len(polygon)
            polygon_center[1] = polygon_center[1] / len(polygon)
            polygon_center[2] = polygon_center[2] / len(polygon)


            # 基準線 なんでもいいっぽい? 
            refelence_vec = polygon_center[:]

            # ポリゴンの中心から各点へのベクトル
            vec_list = []
            for vertice in polygon:
                vec = vertice[0:3] - polygon_center
                vec_list.append(vec)

            # 基準線に対して角度がどれだけあるのか
            theta_list = []

            # 基準線に対して左回りなのか右回りなのか
            vec_z = True
            for index, vec in enumerate(vec_list):
                i = np.inner(refelence_vec, vec)
                n = np.linalg.norm(refelence_vec * np.linalg.norm(vec))
                
                if np.isclose(i,0):
                    c = 0
                else:    
                    c = i/n
                    
                if np.isclose(c,1):
                    c = 1
                elif np.isclose(c,-1):
                    c = -1
                theta = np.rad2deg(np.arccos(c))

                # 最初の基準となるベクトルが外積を求めてどっち回りなのか計算する
                if index == 0:
                    vec_i_n = np.cross(refelence_vec, vec)
                    if vec_i_n[-1] < 0:
                        vec_z = False
                else: # 最初の基準となるベクトルに合わせる
                    vec_i_n = np.cross(refelence_vec, vec)
                    if vec_z: #基準は左回り
                        if vec_i_n[-1] < 0.0:
                            theta = 360 - theta
                    else: #基準は右回り
                        if vec_i_n[-1] > 0.0:
                            theta = 360 - theta
                theta_list.append(theta)

            # 角度の大きいものからソートする
            theta_list = np.array(theta_list)
            theta_index_list = np.argsort(theta_list)[::-1]
            #ソート後の重なってできたポリゴンの頂点
            overlapping_vertices_sorted = []
            for i in theta_index_list:
                overlapping_vertices_sorted.append(polygon[i])
            overlapping_polygons_sorted.append(overlapping_vertices_sorted)
        return overlapping_polygons_sorted
    
    ##
    # @brief 交差点を求める
    # @param vec_r0,vec_r1,vec_r2,vec_r3 2組のポリゴンのそれぞれ２つの頂点
    # @return intersection_vertice_t_s 交差点の位置を示すパラメータ
    def calc_intersection_vertice(self, vec_r0,vec_r1,vec_r2,vec_r3):
        # x,y座標のみ利用します
        x0,y0,z0 = vec_r0
        x1,y1,z1 = vec_r1
        x2,y2,z2 = vec_r2
        x3,y3,z3 = vec_r3
        
        # 交差点の計算
        det = np.linalg.det(np.array([
                [x1-x0, x2-x3],
                [y1-y0, y2-y3]
            ]))
        
        if np.isclose(det,0):
            return [1000,1000] # 正則ではない 範囲から出る値ならなんでもいい
        
        intersection_vertice_t_s = np.linalg.inv(np.array([
                [x1-x0, x2-x3],
                [y1-y0, y2-y3]
            ])) @ np.array([
                x2-x0, y2-y0
            ])
        return intersection_vertice_t_s
    
    ##
    # @brief 交差点の位置を示すパラメータから交差点の状態を判定する
    # @param True/False 交差点の有無
    # @return overlapping_polygons_sorted 重なってできたポリゴンの頂点をソートしたもの
    def judge_intersect(self, intersection_vertice_t_s):
        t = round(intersection_vertice_t_s[0], 6)
        s = round(intersection_vertice_t_s[1], 6)

        # i 延長しても交点なし(並行)
        if (t < 0.0 or 1.0 < t) and (s < 0.0 or 1.0 < s):
            return False
        # ii 交点なし(延長すると交わる)
        elif (0.0 <= t and t <= 1.0) and (s < 0.0 or 1.0 < s):
            return False
        elif (t < 0.0 or 1.0 < t) and (0.0 <= s and s <= 1.0):
            return False
        # iii 片方の線上に交点をもつ
        elif (0.0 < t and t < 1.0) and ( np.isclose(s,0) or np.isclose(s,1) ):
            return True
        elif ( np.isclose(t,0) or np.isclose(t,1) ) and (0.0 < s and s < 1.0):
            return True
        # iv 線分の片方の点が一致する
        elif ( np.isclose(t,0) or np.isclose(t,1) ) and ( np.isclose(s,0) or np.isclose(s,1) ):
            return False
        # v 交点を持つ
        elif (0.0 < t and t < 1.0) and (0.0 < s and s < 1.0):
            return True
        else:
            return False

    ##
    # @brief 複数のポリゴンの面積を求める
    # @param polygons 複数のポリゴンで構成された面
    # @return sum_G 複数のポリゴンで構成された面積
    def calc_pcsa_from_polygons(self, polygons):
        sum_G = 0
        for vertices in polygons:
            G = 0
            for index, vertice in enumerate(vertices):       
                x0 = vertices[index][0] - vertices[0][0]
                y0 = vertices[index][1] - vertices[0][1]
                if index + 1 == len(vertices):
                    x1 = vertices[0][0] - vertices[0][0]
                    y1 = vertices[0][1] - vertices[0][1]
                else:
                    x1 = vertices[index+1][0] - vertices[0][0]
                    y1 = vertices[index+1][1] - vertices[0][1]
                G += (
                    x0 * y1 - x1 * y0
                )
            sum_G += 1/2 * np.abs(G)
        return sum_G
    
    ##
    # @brief 投影したポリゴンの点から線になった要素を省いたもの
    # @param polygons_proj_w 投影したポリゴン
    # @return sum_G polygons_proj_w_ 投影したポリゴンの点から線になった要素を省いたもの
    def remove_line_from_polygons_proj_w(self, polygons_proj_w):
        polygons_proj_w_ = []
        for i, polygon in enumerate(polygons_proj_w):
            flag = []
            for j, vertice in enumerate(polygon):
                if j <= 1:
                    continue
                vec1 = np.array(polygon[j][0:2]) - np.array(polygon[0][0:2])
                vec2 = np.array(polygon[1][0:2]) - np.array(polygon[0][0:2])
                
                n = np.cross(vec1, vec2)
                
                if np.isclose(n, 0.0):
                    flag.append(False)
                else :
                    flag.append(True)
                    
            if np.all(flag):
                polygons_proj_w_.append(polygon)
        return polygons_proj_w_
    
    ##
    # @brief PCSA関連の画像を出力する
    # @param arg argの説明
    # @return None 戻り値なし
    def pcsa_viewer(self):
        # Figureを追加
        fig = plt.figure(figsize=(8,4))
        # 3DAxesを追加
        ax = fig.add_subplot(111) 
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        x = []
        y = []
        z = []
        
        for vertice in self.vertices_proj_w:
            x.append(vertice[0])
            y.append(vertice[1])
            z.append(vertice[2])
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

        ax.grid()
        #ax.set_zlabel('z')

        for index_vertices in self.index_polygons:
            poly = []
            poly.append(self.vertices_proj_w[index_vertices[0]][0:2])
            poly.append(self.vertices_proj_w[index_vertices[1]][0:2])
            poly.append(self.vertices_proj_w[index_vertices[2]][0:2])
            patch = patches.Polygon(xy=poly,closed=True,alpha=0.3,ec='black',fc="#00B0FF")
            ax.add_patch(patch)
        ax.autoscale()
        ax.scatter(x,y, color = "red")
        plt.show()


if __name__ == '__main__':
    pc = PcsaCalculator('../cad/smp_v4.stl',[2,4,2])
    print("PCSA = ", pc.pcsa)
    pc = PcsaCalculator('../cad/box_solid.stl',[0.84143585,0.0091384 ,0.54027975])
    print("PCSA = ", pc.pcsa)
    pc = PcsaCalculator('../cad/box_solid.stl',[0.,0.0091384,0.99995824])
    print("PCSA = ", pc.pcsa)
    pc.pcsa_calc([0.84143585,0.0091384 ,0.54027975])
    print("PCSA = ", pc.pcsa)
    
    
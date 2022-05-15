import math
from pprint import pprint
from matplotlib.pyplot import flag
import numpy as np
##
# @file utils.py
# @author Shuntarou Kuroiwa
# @date 2022/04/05
# @brief pcsa_calculatorを作る上で必要となった関数をまとめてます
# @details 詳細な説明

##
# @brief 配列を要素とする配列 に ある配列 が 含まれているかどうかチェックする関数
# @param main_array sub_array 配列を要素とする配列　　ある配列(チェックしたい要素)
# @return True/False 含まれていたらTrue 含まれていなかったらFalse
def search_array_4_array(main_array, sub_array):
    num_sub_array = len(sub_array)
    for array in main_array:
        # そもそも配列の要素数が違ったら違う
        num_array = len(array)
        if num_array != num_sub_array:
            continue
        
        # 要素が一致しているのかどうかチェックする配列
        flag_array_checker = []
        for i in range(num_array):
            flag_array_checker.append(
                np.isclose(
                    round(array[i], 6),
                    round(sub_array[i], 6)
                )
            )

        # 要素全て一緒の場合　含まれている
        if all(flag_array_checker):
            return True

    return False

##
# @brief ポリゴン(配列を要素とする配列)の配列 に ある配列 が 含まれているかどうかチェックする関数
# @param overlapping_polygons overlapping_vertices ポリゴンの配列　ある配列(チェックしたい要素)
# @return True/False 含まれていたらTrue 含まれていなかったらFalse
def search_polygons_4_polygons_array(overlapping_polygons, overlapping_vertices):
    # 何角形のポリゴン
    num_overlapping_vertices = len(overlapping_vertices)
    for polygon in overlapping_polygons:
        num_polygon = len(polygon)
        if not num_polygon == num_overlapping_vertices:
            continue
        flag_array_checker = []
        for vertice in overlapping_vertices:
            if search_array_4_array(polygon, vertice):
                flag_array_checker.append(True)
            else:
                flag_array_checker.append(False)
                break
        # 含まれている
        if all(flag_array_checker):
            return True
    return False
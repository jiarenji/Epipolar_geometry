from coordinate_conversion import *

  # 角度制转弧度制
PI = 3.1415926;
DEG2RAD = PI / 180.0;
RAD2DEG = 180.0 / PI;

def linalg_solve(R,uv):
    b = np.array([uv[0],uv[1],1]).T
    A = np.mat('{},{},{};{},{},{};{},{},{}'.format(R[0][0],R[0][1],R[0][2],R[1][0],R[1][1],R[1][2],R[2][0],R[2][1],R[2][2]))
    XYZ = np.linalg.solve(A,b)
    print("解得此方程")
    print("X/Z,",XYZ[0])
    print("Y/Z,",XYZ[1])
    return  XYZ


def solve_XYZ(e2x1,e2x2,XYZ1,XYZ2,enu_relative):

    b = np.array([0, 0, (e2x2[0][0]-e2x2[2][0]*XYZ2[0])*enu_relative[0] + (e2x2[0][1]-e2x2[2][1]*XYZ2[0])*enu_relative[1] +  (e2x2[0][2]-e2x2[2][2]*XYZ2[0])*enu_relative[2]]).T
    A = np.mat(
        '{},{},{};{},{},{};{},{},{}'.format(e2x1[0][0]-e2x1[2][0]*XYZ1[0], e2x1[0][1]-e2x1[2][1]*XYZ1[0], e2x1[0][2]-e2x1[2][2]*XYZ1[0],
                                            e2x1[1][0]-e2x1[2][0]*XYZ1[1], e2x1[1][1]-e2x1[2][1]*XYZ1[1], e2x1[1][2]-e2x1[2][2]*XYZ1[1],
                                            e2x2[0][0]-e2x2[2][0]*XYZ2[0], e2x2[0][1]-e2x2[2][1]*XYZ2[0], e2x2[0][2]-e2x2[2][2]*XYZ2[0]))

    XYZ_enu = np.linalg.solve(A, b)
    print(XYZ_enu)
    return XYZ_enu



def main():
     #相机内参矩阵
    K = np.array([[43800 / 15, 0.0, 320.0], [0.0, 43800 / 15, 204.8], [0.0, 0.0, 1.0e+00]])
    print("K:",K)

    #前后帧目标经纬高
    target_blh = [[84.483139,37.329095,10000],[84.483139,37.329095,10000]]
    #前后帧探测器经纬高
    dd_blh = [[84.517328,37.343778,10708],[84.5068,37.342108,10539]]
    #前后帧探测器-航向-俯仰-滚转
    dd_y_p_r = [[123.32,-7.34,2.85],[122.26,-6.98,5.27]]
    #前后帧成像坐标，u：x轴；v:y轴
    uv = [[584,413],[229,468]]
    uv = [[584.42271367,413.03002426],[229,468]]
    
    print('******************************************************************************')
    enu_now = ralative_ENU(dd_blh[0], target_blh[0])
    print(enu_now)
    enu_now = [i / enu_now[2] for i in enu_now]
    print(enu_now)


    
    
    print('******************************************************************************')
    
    yaw_last_frame = dd_y_p_r[0][0];
    pitch_last_frame = dd_y_p_r[0][1];
    roll_last_frame = dd_y_p_r[0][2];
    ypl_last = [yaw_last_frame * DEG2RAD, pitch_last_frame * DEG2RAD, roll_last_frame * DEG2RAD]  # 角度转弧度
    yaw_now_frame = dd_y_p_r[1][0];
    pitch_now_frame = dd_y_p_r[1][1];
    roll_now_frame = dd_y_p_r[1][2];
    yp2_now = [yaw_now_frame * DEG2RAD, pitch_now_frame * DEG2RAD, roll_now_frame * DEG2RAD]  #
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #第一张图所得结果
    e2x1 = ENU2XYZ(ypl_last)
    print('e2x1',e2x1)

    XYZ1 =  UV2XYZ([584,413], K)
    print(XYZ1)

    # 第一张图所得结果
    e2x2 = ENU2XYZ(yp2_now)
    print('e2x2', e2x2)

    XYZ2 = UV2XYZ([229,468], K)
    print(XYZ2)

    enu_relative = ralative_ENU(dd_blh[0], dd_blh[1])
    print(enu_relative)

    XYZ_ENU =  solve_XYZ(e2x1, e2x2, XYZ1, XYZ2, enu_relative)
    print(XYZ_ENU)




    print('******************************************************************************')
    
    
    
    # print("目标经纬高:", target_blh[0])  
    # print("探测器1经纬高:",dd_blh[0])
    # print("探测器-航向-俯仰-滚转",dd_y_p_r[0])
    # print("图像坐标",uv[0])
    
    # 求探测器1,2东北天
    # 经纬高转东北天坐标
    # dd_blh[0]为上一帧探测器经纬高，记为东北天（相对系）原点； 
    # dd_blh[1]为当前帧探测器经纬高，计算当前帧探测器相对于上一帧探测器的东北天坐标
    # enu_now = ralative_ENU(dd_blh[0], dd_blh[1])
    
    #求R_1
    # 前一帧探测器-航向-俯仰-滚转
    yaw_last_frame = dd_y_p_r[0][0];
    pitch_last_frame = dd_y_p_r[0][1];
    roll_last_frame = dd_y_p_r[0][2];
    ypl_last = [yaw_last_frame * DEG2RAD, pitch_last_frame * DEG2RAD, roll_last_frame * DEG2RAD]  # 角度转弧度
    R_1 = ENU2UV(K,ypl_last)
    print('R1',R_1)
    #上一帧目标在图像中的坐标
    uv_last = [584.42271367,413.03002426]
    print("上一帧")
    print(uv_last)
    print(R_1)
    XY_1 =  linalg_solve_22(R_1, uv_last)
    print(XY_1)
    
    U_A3 = R_1[0][0] * 4.280822275315625 + R_1[0][1] * 2.301487038916918
    print("U_A3", U_A3)
    U_A31 = uv_last[0] - R_1[0][2] 
    print("U_A31", U_A31)
    
    print('******************************************************************************')
    
    
    
    
    
    # R_2 
    # 当前帧探测器-航向-俯仰-滚转
    yaw_now_frame = dd_y_p_r[1][0];
    pitch_now_frame = dd_y_p_r[1][1];
    roll_now_frame = dd_y_p_r[1][2];
    
    yp2_now = [yaw_now_frame * DEG2RAD, pitch_now_frame * DEG2RAD, roll_now_frame * DEG2RAD]  # 
    print(yp2_now)
    R_2 = ENU2UV(K,yp2_now)
    uv_now = [229,468]
    print("当前帧")
    XY_2 =  linalg_solve_22(R_2, uv_now)
    print(XY_2)
    
    
    U_A3 = R_2[0][0] * XY_2[0] + R_2[0][1] * XY_2[1]  
    print("U_A3", U_A3)
    U_A31 = uv_now[0] - R_2[0][2] 
    print("U_A31", U_A31)
    
    print('******************************************************************************')
    
    
    
        
    enu_ddnow = ralative_ENU(dd_blh[0], dd_blh[1])   
    print(enu_ddnow)
    print('******************************************************************************')
    
    z = (-XY_2[0]*enu_ddnow[2] + enu_ddnow[0])/(XY_1[0] - XY_2[0])
    x = XY_1[0] * z
    print('z',z)
    print('x',x)
    print('******************************************************************************')
    
    
    
    
    
    # 求目标的东北天坐标
    A = np.array([[1,1,-(XY_1[0]+XY_1[1])],
                  [1,0,-XY_2[0]],
                  [0,1,-XY_2[1]]])
    
    b = np.array([0,
                  enu_now[0]-XY_2[0]*enu_now[2],
                  enu_now[1]-XY_2[1]*enu_now[2]])
    
    xyz = np.linalg.solve(A,b)
    
    print("线性方程组的解 x, y, z:")
    print(xyz)
    



if __name__ == '__main__':
    main()
  

   
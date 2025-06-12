#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
@Project ：Epipolar-Geometry-master
@File    ：main.py
@IDE     ：PyCharm
@Author  ：jijiaren
@Date    ：2024-04-16 9:26
'''


import numpy as np
import cv2 as cv
from matplotlib import pyplot as plt
import math


def get_distance_point2line(point, line):
    """
    Args:
        point: [x0, y0]
        line: [x1, y1, x2, y2]
    """
    line_point1, line_point2 = np.array(line[0:2]), np.array(line[2:])
    vec1 = line_point1 - point
    vec2 = line_point2 - point
    distance = np.abs(np.cross(vec1, vec2)) / np.linalg.norm(line_point1 - line_point2)
    return distance



def drawlines(img1,lines, pts,pts2):
    if len(img1.shape) > 2:
        r,c,_ = img1.shape
    else:
        r, c= img1.shape
    for r, pt1 in zip(lines, pts):
        color=(0,255,255)
        x0,y0 = map(int, [0, -r[2]/r[1] ])
        x1,y1 = map(int, [c, -(r[2]+r[0]*c)/r[1] ])
        img1 = cv.line(img1.copy(), (x0,y0), (x1,y1), color ,4)
        img1 = cv.circle(img1, tuple(pt1), 5, color, 5)
        print("当前图像目标真实坐标:",pt1)
        print("估计直线",[x0,y0,x1,y1])
        distance = get_distance_point2line(pt1, [x0,y0,x1,y1])
        print("真实坐标距离估计直线距离error",distance)

    # for pt2 in  pts2:
    #     img1 = cv.circle(img1, tuple(pt2), 5, color, 5)
    return img1

def drawcircle(img1, pts):
    if len(img1.shape) > 2:
        r,c,_ = img1.shape
    else:
        r, c= img1.shape

    for  pt1 in  pts:
        color=(0,255,255)
        img1 = cv.circle(img1, tuple(pt1), 5, color, 5)
    return img1


"""
经纬高坐标系转地心坐标系
"""
def llh2ecef(LBH):
    """
    LBH:经纬高
    """

    a = 6378137.0;
    b = 6356752.314;

    lon = LBH[0] * 3.1415926 / 180.0;

    lat = LBH[1] * 3.1415926 / 180.0;

    alt = LBH[2] ;

    n = a * a / math.sqrt(a * a * math.cos(lat) * math.cos(lat) + b * b * math.sin(lat) * math.sin(lat));
    Rx = (n + alt) * math.cos(lat) * math.cos(lon);
    Ry = (n + alt) * math.cos(lat) * math.sin(lon);
    Rz = (b * b / (a * a) * n + alt) * math.sin(lat);


    Rxyz = [Rx,Ry,Rz]

    return Rxyz

"""
地心坐标系转经纬高坐标系   
"""
def ecef2llh(Rxyz):

    pi = 3.1415926;

    x = Rxyz[0];
    y = Rxyz[1];
    z = Rxyz[2];

    x2 = pow(x, 2);
    y2 = pow(y, 2);
    z2 = pow(z, 2);

    a = 6378137.0000;
    b = 6356752.3142;
    e = math.sqrt(1 - (b / a) * (b / a));

    b2 = b * b;
    e2 = e * e;
    ep = e * (a / b);
    r = math.sqrt(x2 + y2);
    r2 = r * r;
    E2 = a * a - b * b;
    F = 54 * b2 * z2;
    G = r2 + (1 - e2) * z2 - e2 * E2;
    c = (e2 * e2 * F * r2) / (G * G * G);
    s = (1 + c + math.sqrt(c * c + 2 * c));
    s = pow(s, 1 / 3);

    P = F / (3 * ((s + 1 / s + 1) * (s + 1 / s + 1)) * G * G);
    Q = math.sqrt(1 + 2 * e2 * e2 * P);
    ro = -(P * e2 * r) / (1 + Q) + math.sqrt((a * a / 2) * (1 + 1 / Q) - (P * (1 - e2) * z2) / (Q * (1 + Q)) - P * r2 / 2);
    tmp = (r - e2 * ro) * (r - e2 * ro);
    U = math.sqrt(tmp + z2);
    V = math.sqrt(tmp + (1 - e2) * z2);
    zo = (b2 * z) / (a * V);

    height = U * (1 - b2 / (a * V));
    lat = math.atan((z + ep * ep * zo) / r);
    temp = math.atan(y / x);

    if (x >= 0):
        long_ = temp;
    elif ((x < 0) and (y >= 0)):
        long_ = pi + temp;
    else:
        long_ = temp - pi;
    llh = [(long_) * (180 / pi),(lat) * (180 / pi),height]

    return llh;


"""
地心坐标系转东北天坐标系ENU（X,Y,Z）
"""
def ecef2enu(LBH,Rxyz):

    pi = 3.1415926;
    DEG2RAD = pi / 180.0;
    RAD2DEG = 180.0 / pi;

    x = Rxyz[0];
    y = Rxyz[1];
    z = Rxyz[2];

    Rxyz_1 = llh2ecef(LBH)
    ox = Rxyz_1[0];
    oy = Rxyz_1[1];
    oz = Rxyz_1[2];

    dx = x - ox;
    dy = y - oy;
    dz = z - oz;


    lonDeg = LBH[0];
    latDeg = LBH[1];
    lon = lonDeg * DEG2RAD;
    lat = latDeg * DEG2RAD;

    enu_0 = -math.sin(lon) * dx + math.cos(lon) * dy;
    enu_1 = -math.sin(lat) * math.cos(lon) * dx - math.sin(lat) * math.sin(lon) * dy + math.cos(lat) * dz;
    enu_2 = math.cos(lat) * math.cos(lon) * dx + math.cos(lat) * math.sin(lon) * dy + math.sin(lat) * dz;

    enu = [enu_0,enu_1,enu_2]
    return enu;


def vec_inv(vec):

    vecinv = [[0 , -vec[2],  vec[1]],[vec[2],  0,  -vec[0]],[-vec[1], vec[0],  0]]

    return np.array(vecinv)


def calculate_projected_line(F,pts1,pts2,img1_path,img2_path):

    pts1 = np.array(pts1)
    pts2 = np.array(pts2)
    img1 = plt.imread(img1_path)
    img2 = plt.imread(img2_path)


    ones = np.ones((1,1))
    pts1_final = np.append(pts1, ones, axis=1)

    lines1 = np.dot(F, pts1_final.transpose())


    # ones = np.ones((1,1))
    # pts2_final = np.append(pts2, ones, axis=1)
    #
    # lines2 = np.dot(pts2_final, F)

    img3 = drawlines(img2, lines1.transpose(), pts2,pts1)
    plt.imshow(img3)
    plt.savefig("{}.png".format(img2_path.split(".")[0]+'_new'))
    # plt.imshow(img3)
    # plt.show()
    # img4 = drawcircle(img1,pts1)
    # plt.imshow(img4)
    # plt.savefig("{}.png".format(img1_path.split(".")[0]+'_new'))
    # plt.imshow(img4)
    # plt.show()

    # img4 = drawlines(img1, lines2, pts1)
    # plt.imshow(img4)
    # plt.show()



def isRotationMatrix(R):
    should_identity = R * R.transpose();
    I = np.array([[1,0,0],[0,1,0],[0,0,1]])
    err = 1e-3;

    return (R - I)


def computrdUVPos(XYZ,FOVX,FOVY):
    """
    计算3D空间点(X,Y,Z)在给定FOVX和FOVY的视野下，在2D图像上的投影位置(u,v)。
    
    Args:
        XYZ (list[float]): 包含三个元素的列表，分别表示3D空间点的X、Y、Z坐标。
        FOVX (float): 视野的X轴角度大小，单位是度。
        FOVY (float): 视野的Y轴角度大小，单位是度。
    
    Returns:
        list[float]: 包含两个元素的列表，分别表示2D图像上的投影位置u和v。
    
    """
    pi = 3.1415926535
    halfsizex = math.tan(FOVX * 0.5 * pi / 180)

    X,Y,Z = XYZ[0],XYZ[1],XYZ[2]
    u = (X/Y)/halfsizex

    if Z < 0:
        halfsizey = math.tan(FOVY * 0.6 * pi / 180)
        v = (Z/Y)/halfsizey
    else:
        halfsizey = math.tan(FOVY * 0.4 * pi / 180)
        v = (Z/Y)/halfsizey

    u = 320 + u*320
    if v < 0:
        v = 204.8 - v * 307.2
    else :
        v = 204.8 - v * 204.8

    return [u,v]



def vaild_point(EulerAngles,LBH_dd,LBH_plane,uv):
    # 知道目标经纬高，求uv
    pi = 3.1415926;

    DEG2RAD = pi / 180.0;
    print(f"探测器经纬高:{LBH_dd}")
    print(f"目标经纬高:{LBH_plane}")


    EulerAngles = [i*DEG2RAD for i in EulerAngles]  # 探测器的欧拉角
    print("欧拉角弧度",EulerAngles[0],EulerAngles[1],EulerAngles[2])
    Zyaw = np.array([[math.cos(EulerAngles[0]),math.sin(EulerAngles[0]),0],[-math.sin(EulerAngles[0]),math.cos(EulerAngles[0]),0],[0,0,1]]);
    Xpitch = np.array([[1,0,0],[0,math.cos(EulerAngles[1]),math.sin(EulerAngles[1])],[0,-math.sin(EulerAngles[1]),math.cos(EulerAngles[1])]]);
    Yroll = np.array([[math.cos(EulerAngles[2]),0,-math.sin(EulerAngles[2])],[0,1,0],[math.sin(EulerAngles[2]),0,math.cos(EulerAngles[2])]]);
    R_yanzheng = np.matmul(np.matmul(Yroll,Xpitch),Zyaw)  # 探测器的欧拉
    # R_yanzheng_wai = np.matmul(np.matmul(Zyaw,Xpitch),Yroll)

    Rxyz_last = llh2ecef(LBH_dd)  # 探测器经纬高坐标系转地心坐标系
    enu_last = ecef2enu(LBH_dd, Rxyz_last)  # 探测器地心坐标系转东北天坐标系ENU（X,Y,Z）
    print("探测器东北天坐标",enu_last)

    Rxyz_now = llh2ecef(LBH_plane)  # 目标经纬高坐标系转地心坐标系
    enu_now = ecef2enu(LBH_dd, Rxyz_now)  # 目标地心坐标系转东北天坐标系ENU（X,Y,Z）
    print("目标东北天坐标", enu_now)
    print("目标与探测器水平距离", math.sqrt(enu_now[0] * enu_now[0] + enu_now[1] * enu_now[1]))  # 单位-- m。求得相对于0，0，0的水平距离


    plane_vec = np.array(enu_now)   # 目标的相对于探测器1的东北天坐标
    XYZ_vector = np.matmul(R_yanzheng,plane_vec)  # 目标的东北天--》导航系
    # XYZ_vector1 = np.matmul(R_yanzheng_wai,plane_vec)
    print("目标相对与探测器坐标",XYZ_vector)   # 导航系下的
    # print("验证坐标系旋转",XYZ_vector1)
    XYZ2XZY = np.matmul(np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]), XYZ_vector)  # 导航系--》相机系
    print("验证XYZ坐标",XYZ2XZY)
    XYZ2XZY = XYZ2XZY / XYZ2XZY[2]   # 转齐次坐标
    print("验证XYZ坐标,转齐次坐标",XYZ2XZY)
    print("目标在图像中的真实uv",uv)
    uv_vector_fromK = np.matmul(K,XYZ2XZY)  # 相机系--》像素系


    uv_vector_fromratio = computrdUVPos(XYZ_vector, 4.569807720109516, 3.6565438793883596)   # 12.5, 10 如何得到？？ 解答，仿真软件的成像，用于标定
    print("验证XYZ转UV,uv_vector_fromK",uv_vector_fromK)
    print("验证XYZ转UV,uv_vector_fromratio",uv_vector_fromratio)



def oula_xuanzhuan(EulerAngles):
    Zyaw = np.array([[math.cos(EulerAngles[0]),math.sin(EulerAngles[0]),0],[-math.sin(EulerAngles[0]),math.cos(EulerAngles[0]),0],[0,0,1]]);
    Xpitch = np.array([[1,0,0],[0,math.cos(EulerAngles[1]),math.sin(EulerAngles[1])],[0,-math.sin(EulerAngles[1]),math.cos(EulerAngles[1])]]);
    Yroll = np.array([[math.cos(EulerAngles[2]),0,-math.sin(EulerAngles[2])],[0,1,0],[math.sin(EulerAngles[2]),0,math.cos(EulerAngles[2])]]);
    R =  np.matmul(np.matmul(Yroll, Xpitch), Zyaw)
    return R


def ralative_ENU(LBH_last,LBH_now):
    """
    已知前后两帧目标经纬高 求 当前帧基于上一帧目标的东北天坐标
    :param LBH_last: 上一帧经纬高
    :param LBH_now: 当前帧经纬高
    :return: 当前帧相对于上一帧的东北天坐标系下坐标
    """
    Rxyz_last = llh2ecef(LBH_last)
    enu_last = ecef2enu(LBH_last, Rxyz_last)
    print("上一帧探测器东北天坐标", enu_last)
    print("目标基于上一帧探测器东北天坐标系坐标：X Y Z")
    Rxyz_now = llh2ecef(LBH_now)
    enu_now = ecef2enu(LBH_last, Rxyz_now)
    print("当前帧探测器东北天坐标", enu_now)
    print("目标基于当前帧探测器东北天坐标系坐标：X-{} Y-{} Z-{}".format(enu_now[0],enu_now[1],enu_now[2]))

    return enu_now

def ENU2UV(K,y_p_r):
    """
    东北天->导航系->图像系->坐标
    oula_xuanzhuan_tmp,R_dd2image,K

    :param K: 内参矩阵，只与内参相关
    :param y_p_r: 探测器航向，俯仰，滚转
    :return:
    """
    oula_xuanzhuan_tmp = oula_xuanzhuan(y_p_r)
    R_dd2image = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    return  np.matmul(np.matmul(K,R_dd2image),oula_xuanzhuan_tmp)


def ENU2XYZ(y_p_r):
    """
    东北天->导航系->图像系
    oula_xuanzhuan_tmp,R_dd2image,K

    :param y_p_r: 探测器航向，俯仰，滚转
    :return:
    """
    oula_xuanzhuan_tmp = oula_xuanzhuan(y_p_r)
    R_dd2image = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    return  np.matmul(R_dd2image,oula_xuanzhuan_tmp)


def UV2XYZ(UV,K):
    """
    图像系->相机系
    oula_xuanzhuan_tmp,R_dd2image,K

    :param y_p_r: 探测器航向，俯仰，滚转
    :return:
    """
    K_inv = np.linalg.inv(K)
    UV = np.array([UV[0],UV[1],1])
    XYZ = np.matmul(K_inv,UV)



    return  XYZ


def linalg_solve(R,uv):
    b = np.array([uv[0],uv[1],1]).T
    A = np.mat('{},{},{};{},{},{};{},{},{}'.format(R[0][0],R[0][1],R[0][2],R[1][0],R[1][1],R[1][2],R[2][0],R[2][1],R[2][2]))
    XYZ = np.linalg.solve(A,b)
    print("解得此方程")
    print("X/Z,",XYZ[0])
    print("Y/Z,",XYZ[1])
    return  XYZ



def linalg_solve_22(R,uv):
    b = np.array([uv[0]-R[0][2],uv[1]-R[1][2]]).T
    A = np.mat('{},{};{},{}'.format(R[0][0],R[0][1],R[1][0],R[1][1]))
    XYZ = np.linalg.solve(A,b)
    print("解得此方程")
    print("X/Z,",XYZ[0])
    print("Y/Z,",XYZ[1])
    return  XYZ



def ENU2BLH():


    return 0



if __name__ == '__main__':


    """
    # R = np.array([[1, 2,1], [4.5, 6,1],[1,1,1]])
    # b = [3, 10.5]   
    # XY =  linalg_solve(R, b)
    """

    pi = 3.1415926;
    DEG2RAD = pi / 180.0;
    RAD2DEG = 180.0 / pi;

    #相机内参矩阵
    K = np.array([[8020, 0.0, 320.0], [0.0, 8020, 204.8], [0.0, 0.0, 1.0e+00]])
    print("K:",K)

    #前后帧目标经纬高
    target_blh = [[120.32159772,36.054907059,12.049120136],[120.32159772,36.054907059,12.049120136]]
    #前后帧探测器经纬高
    dd_blh = [[120.279849667912,36.058055505632,1014.605175501667],[120.279842741550,36.057998310434,1014.788692085072]]
    #前后帧探测器-航向-俯仰-滚转
    dd_y_p_r = [[95.45522268672802,-14.94768108101971,11.851206010414064],[95.29337752019588,-15.019529378031152,12.224938390969163]]
    #前后帧成像坐标，u：x轴；v:y轴
    uv = [[321,251],[317,252]]

    # 经纬高转东北天坐标
    # dd_blh[0]为上一帧探测器经纬高，记为东北天（相对系）原点； 
    # dd_blh[1]为当前帧探测器经纬高，计算当前帧探测器相对于上一帧探测器的东北天坐标
    enu_now = ralative_ENU(dd_blh[0], dd_blh[1])   
    
    
    #前后帧图片路径
    image1_path = r'/workspace/epipolar_geometry/000000.bmp'
    image2_path = r'/workspace/epipolar_geometry/000001.bmp'
    
    
    print("计算验证成像的正确性*************************************************************************************")
    for i in range(2):
        print("验证投影点********************************************************")
        this_target_blh = np.array(target_blh[i])
        this_dd_blh = np.array(dd_blh[i])
        this_dd_y_p_r = np.array(dd_y_p_r[i])
        this_uv = np.array(uv[i])
        print("目标经纬高:", this_target_blh)
        print("探测器经纬高:",this_dd_blh)
        print("探测器-航向-俯仰-滚转",this_dd_y_p_r)
        print("图像坐标",this_uv)
        vaild_point(this_dd_y_p_r, this_dd_blh, this_target_blh, this_uv)
        print("验证投影点********************************************************")
    
    
    print("计算验证投影直线算法*****************************************************************************************")
    # 上一帧目标经纬高
    l_plane = target_blh[0][0];
    b_plane = target_blh[0][1];
    h_plane = target_blh[0][2];
    LBH_plane = [l_plane,b_plane,h_plane]
    
    # 前一帧帧探测器-航向-俯仰-滚转
    yaw_last_frame = dd_y_p_r[0][0];
    pitch_last_frame = dd_y_p_r[0][1];
    roll_last_frame = dd_y_p_r[0][2];
    ypl_last = [yaw_last_frame * DEG2RAD, pitch_last_frame * DEG2RAD, roll_last_frame * DEG2RAD]  # 角度转弧度
    print("前一帧探测器机体坐标系(航向，俯仰，滚转)：", yaw_last_frame, pitch_last_frame, roll_last_frame)
    # 当前帧探测器-航向-俯仰-滚转
    yaw_now_frame = dd_y_p_r[1][0];
    pitch_now_frame = dd_y_p_r[1][1];
    roll_now_frame = dd_y_p_r[1][2];
    ypl_now = [yaw_now_frame * DEG2RAD, pitch_now_frame * DEG2RAD, roll_now_frame * DEG2RAD]
    print("当前帧探测器机体坐标系(航向，俯仰，滚转)：", yaw_now_frame, pitch_now_frame, roll_now_frame)
    l_last_frame = dd_blh[0][0];
    b_last_frame = dd_blh[0][1];
    h_last_frame = dd_blh[0][2];
    LBH_last = [l_last_frame, b_last_frame, h_last_frame]
    print("前一帧探测器机体坐标(经度，维度，高度)：", l_last_frame, b_last_frame, h_last_frame)
    # 当前帧探测器经纬高
    l_now_frame = dd_blh[1][0];
    b_now_frame = dd_blh[1][1];
    h_now_frame = dd_blh[1][2];
    LBH_now = [l_now_frame, b_now_frame, h_now_frame]
    print("当前帧探测器机体坐标(经度，维度，高度)：", l_now_frame, b_now_frame, h_now_frame)
    EulerAngles = [(yaw_now_frame - yaw_last_frame) * DEG2RAD, (pitch_now_frame - pitch_last_frame) * DEG2RAD, 0.0]  # 此处roll为什么是0.0？？
    print("欧拉角弧度", EulerAngles[0], EulerAngles[1], EulerAngles[2])
    # Zyaw = np.array([[math.cos(EulerAngles[0]),math.sin(EulerAngles[0]),0],[-math.sin(EulerAngles[0]),math.cos(EulerAngles[0]),0],[0,0,1]]);
    # Xpitch = np.array([[1,0,0],[0,math.cos(EulerAngles[1]),math.sin(EulerAngles[1])],[0,-math.sin(EulerAngles[1]),math.cos(EulerAngles[1])]]);
    # Yroll = np.array([[math.cos(EulerAngles[2]),0,-math.sin(EulerAngles[2])],[0,1,0],[math.sin(EulerAngles[2]),0,math.cos(EulerAngles[2])]]);
    # R_yanzheng = np.matmul(np.matmul(Yroll,Xpitch),Zyaw)
    # R = R_Yroll_Xpitch_Zyaw = np.matmul(np.matmul(Yroll, Xpitch), Zyaw)
    
    Rxyz_last = llh2ecef(LBH_last)   # 探测器1经纬高坐标系转地心坐标系
    enu_last = ecef2enu(LBH_last, Rxyz_last)    # 探测器1地心坐标系转东北天坐标系ENU（X,Y,Z）
    print("上一帧探测器1东北天坐标", enu_last)
    Rxyz_now = llh2ecef(LBH_now)    #探测器2经纬高坐标系转地心坐标系  
    t = enu_now = ecef2enu(LBH_last, Rxyz_now)  # 探测器2地心坐标系转东北天坐标系ENU（X,Y,Z）
    print("当前帧探测器2东北天坐标", enu_now)
    
    #开始求R,T
    #x2 = Rdd2图 * R欧2 * [（Rdd2图 * R欧1）逆 * x1  +  北天冬delta_t  ]
    oula_xuanzhuan1 = oula_xuanzhuan(ypl_last)  # 上一帧欧拉角 求R
    oula_xuanzhuan1_inv = np.linalg.inv(oula_xuanzhuan1)
    oula_xuanzhuan2 = oula_xuanzhuan(ypl_now)   #当前帧欧拉角 求R
    R_dd2image = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])  # 导航系--》相机系 R
    R_dd2image_inv = np.linalg.inv(np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]))
    print("上一帧欧拉旋转矩阵1",oula_xuanzhuan1)
    print("当前帧欧拉旋转矩阵2",oula_xuanzhuan2)
    print("R_dd2image",R_dd2image)
    
    R = np.matmul(np.matmul(np.matmul(R_dd2image,oula_xuanzhuan2),oula_xuanzhuan1_inv),R_dd2image_inv)
    t = np.matmul(np.matmul(R_dd2image,oula_xuanzhuan2),-np.array(t))
    
    
    
    print("求得R,t为*********************************************************")
    print("R:",R)
    print("t:",t)
    print("*********************************************************")
    
    print("F = K逆的转置 * t反对称 * R * K逆")
    
    F1 = np.linalg.inv(K).T
    F2 = vec_inv(t)
    F3 = R
    F4 = np.linalg.inv(K)
    F = np.matmul(np.matmul(np.matmul(F1,F2),F3),F4)
    print(F)
    print("*********************************************************")
    
    pts1 = [uv[0]]
    pts2 = [uv[1]]
    calculate_projected_line(F,pts1,pts2,image1_path,image2_path)
    
    







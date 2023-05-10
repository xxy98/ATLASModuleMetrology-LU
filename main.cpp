#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#define ct
#include <ctime>
#define ch
#include <chrono>
#define cmath
#include <cmath>
#include <Eigen/Dense>
#define eigen
#include "src/tz.cpp"
#include "date/date.h"
#include "date/tz.h"
#include "date/tz_private.h"
#define HHdate

//#define debug
//#define t1017
//#define t1118
//#define KU1215
//#define t0207
//#define t0207z
//#define t02151p
//#define t021555
//#define t0227
//#define t0317probe
//#define UU0328_1
//#define UU0328_2
//#define UU0329CrossF
//#define UU0329ChequeF
//#define UU0329FF
//#define UU0329Zring168
//#define UU0329ZTTL162
//#define UU0329ZTTL180
//#define UU0330_1
//#define UU0330_2
#define UU0330_3
//#define LU0502_1_Op
//#define LU0502_2_Op
//#define LU0502_3_Op
//#define LU0502_1_TP
//#define LU0502_2_TP
//#define LU0502_3_TP

#if defined(debug) || defined(t1017) || defined(t1118) || defined(t0207) || defined(t0207z) || defined(t02151p) || defined(t021555) || defined(t0227) || defined(t0317probe) || \
    defined(LU0502_1_Op) || defined(LU0502_2_Op) || defined(LU0502_3_Op) || defined(LU0502_1_TP) || defined(LU0502_2_TP) || defined(LU0502_3_TP)
    #define LU
#endif
#if defined(KU1215)
    #define KU
#endif
#if defined(UU0328_1) || defined(UU0328_2) || defined(UU0329CrossF) || defined(UU0329ChequeF) || defined(UU0329FF) || defined(UU0329Zring168) || defined(UU0329ZTTL162) || \
    defined(UU0329ZTTL180) || defined(UU0330_1) || defined(UU0330_2) || defined(UU0330_3)
    #define UU
#endif

#define PI 3.1415926

#include "structures.h"
#define ANA_STRUCTURES_H
#include "mathematics.h"
#define ANA_MATHEMATICS_H
#include "measurements.h"
#define ANA_MEASUREMENTS_H
#include "functions.h"
#define ANA_FUNCTIONS_H

int main() {
#ifdef LU
    std::cout << "Lunds Universitet configuration applied.\n";
#endif
#ifdef KU
    std::cout << "Københavns Universitet configuration applied.\n";
#endif
#ifdef UU
    std::cout << "Uppsala Universitet configuration applied.\n";
#endif

    int num;
    int runt;
    static std::string MeasProgVers = "V1";
    static std::string work_dir = "C:/Users/Xiangyu/Documents/LU/EPP/Ana/";
    static std::string file_dir = "../ATLASModuleMetrology/";

#ifdef debug
    static std::string file_name = "Xu-20221017-condensed.TXT";
#endif
#ifdef t1017
    static std::string file_name = "Xu-20221017-condensed.TXT";
#endif
#ifdef t1118
    static std::string file_name = "Xu-20221118-condensed.TXT";
#endif
#ifdef KU1215
    static std::string file_name = "KU-test-20221215-03.txt";
#endif
#ifdef t0207
    static std::string file_name = "Xu-20230206-condensed-nbi-r.txt";
#endif
#ifdef t0207z
    static std::string file_name = "Xu-20230206-condensed-z.txt";
#endif
#ifdef t02151p
    static std::string file_name = "Xu-20230215-condensed-z-1p.txt";
#endif
#ifdef t021555
    static std::string file_name = "Xu-20230215-condensed-z-55.txt";
#endif
#ifdef t0227
    static std::string file_name = "Xu-20230227-condensed-z-2.txt";
#endif
#ifdef t0317probe
    static std::string file_name = "Xu-20230317-condensed-probe.txt";
#endif
#ifdef UU0328_1
    static std::string file_name = "UU-20230328-full.mxi Listing.txt";
#endif
#ifdef UU0328_2
    static std::string file_name = "UU-20230328-full.mxi Listing-2.txt";
#endif
#ifdef UU0329CrossF
    static std::string file_name = "UU-20230329-+test3.txt";
#endif
#ifdef UU0329ChequeF
    static std::string file_name = "UU-20230329-chequetest3.txt";
#endif
#ifdef UU0329FF
    static std::string file_name = "UU-20230329-Ftest14-4-3-4-3.txt";
#endif
#ifdef UU0329Zring168
    static std::string file_name = "UU-20230329-z-ring168-3x3.txt";
#endif
#ifdef UU0329ZTTL162
    static std::string file_name = "UU-20230329-z-ttl162-3x3.txt";
#endif
#ifdef UU0329ZTTL180
    static std::string file_name = "UU-20230329-z-ttl180-3x3.txt";
#endif
#ifdef UU0330_1
    static std::string file_name = "UU-20230330-full-ttl180.mxi Listing-c.txt";
#endif
#ifdef UU0330_2
    static std::string file_name = "UU-20230330-full-ttl180.mxi Listing-c2.txt";
#endif
#ifdef UU0330_3
    static std::string file_name = "UU-20230330-full-ttl180.mxi Listing-c3.txt";
#endif
#ifdef LU0502_1_Op
    static std::string file_name = "Xu-202303502-01.1-condensed-optical.TXT";
#endif
#ifdef LU0502_2_Op
    static std::string file_name = "Xu-202303502-02-condensed-optical.TXT";
#endif
#ifdef LU0502_3_Op
    static std::string file_name = "Xu-202303502-03-condensed-optical.TXT";
#endif
#ifdef LU0502_1_TP
    static std::string file_name = "Xu-202303502-01-condensed-probe.TXT";
#endif
#ifdef LU0502_2_TP
    static std::string file_name = "Xu-202303502-02-condensed-probe.TXT";
#endif
#ifdef LU0502_3_TP
    static std::string file_name = "Xu-202303502-03-condensed-probe.TXT";
#endif

    test_path(work_dir + file_dir, file_name, false);
    num = get_line_num(work_dir + file_dir, file_name);

    std::ifstream fi(work_dir + file_dir + file_name, std::ios::in);
    std::string line;
    int count = 0;

    Measurement st[num];
    for (int i = 0; i < num; i++) {
        st[i].Set(0, EMPTY, 0, 0, 0, 0, -1, -1);
    }

    while (getline(fi, line)){
        count++;
#if defined(LU) || defined(UU)
        if (count > 3 && count < 4 + num){
            st[count-4] = Meas_Input(line);
            //st[count-4].Output();
#endif
#if defined(KU)
        if (count > 2 && count < 3 + num){
            st[count - 3] = Meas_Input(line);
            //st[count - 3].Output();
#endif
        }
    }
    fi.close();

    std::cout << "======== Here ends the file reading ========\n";

    static std::string out_folder = "Output/";
    time_t now = std::time(0);
    tm *now_l = std::localtime(&now);
    char now_chr[20];
    std::string now_str;
    strftime(now_chr, 64, "%Y-%m-%d-%H-%M-%S", now_l);
    now_str = now_chr;
    static std::string out_name = "Output-"+now_str+".txt";

    std::cout << "======== Here starts the file writing ========\n";

    int ECorB = 0; // 0 for End Cap "SE", 1 for Barrel "SB"
    int ModTyp = 1; // ML, MS, M0-M2, MA-MF
    std::string ModTypS = "";
    std::string ModRef = "00000000";
    std::string institute = "";
    std::string operatorp = "";
    std::string instrutyp = "";

    //std::ofstream fo(work_dir + file_dir + out_folder + out_name, std::ios::out);
    std::ofstream fo;
    fo.open(work_dir + file_dir + out_folder + out_name, std::ios::out);
    //fo.open(work_dir + file_dir + out_folder + out_name, std::ios::app);
    if (!fo) {std::cout<<"Error! Output file can't be created.\n";}

    double hxp1x, hxp1y, hxp2x, hxp2y, pb1x, pb1y, pb2x, pb2y;

    class Point PB1, PB2;
    std::vector<class Point> HyXY;
    //std::vector<class PointsWDistr> PB;
    std::vector<class Point> PB;
    std::vector<class Point> Si;
    std::vector<std::vector <class Point>> ABC_1_v_t, ABC_2_v_t, ABC_3_v_t, ABC_4_v_t;
    std::vector<std::vector <class Point>> HCC_1_v_t, HCC_2_v_t, HCC_3_v_t, HCC_4_v_t;
    std::vector<std::vector<std::vector <class Point>>> ABC, HCC;

#ifdef debug
    /*
    class Point p1, p2, p3, p4;
    class Line l1, l2;
    p1.SetPoint(0,0,-1);
    p2.SetPoint(1,1,0);
    p3.SetPoint(0,1,0);
    p4.SetPoint(2,-1,1);
    l1 = LinkPoints(p1, p2);
    l2 = LinkPoints(p3, p4);

    class Vector IP1, IP2;
    IP1 = LineDistVec(l1,l2,true);
    IP2 = LineDistVec(l1,l2,false);
    std::cout<<LineDist(l1,l2)<<std::endl;
    IP1.Out(true,true,true,true);
    std::cout<<std::endl;
    IP2.Out(true,true,true,true);
    */

    /*
    std::vector<class Point>pv1;
    pv1.push_back(p1);
    pv1.push_back(p2);
    pv1.push_back(p3);
    pv1.push_back(p4);
    class PointsWDistr pd;
    pd = PointsAve(pv1);
    std::cout <<
        " x(ave):" << pd.x <<
        " y(ave):" << pd.y <<
        " z(ave):" << pd.z <<
        " x(SD):" << pd.xSD <<
        " y(SD):" << pd.ySD <<
        " z(SD):" << pd.zSD <<
        std::endl;
    */
    //double ang = VecsAngle(l1.vec, l2.vec);
    //std::cout << ang/2/3.1415926*360;
    //std::cout << pow(PointLineDist(p3, l1),2)<<std::endl;
    //std::cout << LineDist(l1, l2) << std::endl;
    //class Point ic = TwoLineIntersection(l1, l2);
    //ic.Out();

    //class Vector v;
    //v = LineDistVec(l1, l2);
    //v.Out(true,true,true,true);

    class Point p[4];
    class Plane pl[3];
    p[0].SetPoint(0,0,0);
    p[1].SetPoint(1,1,0);
    p[2].SetPoint(1,0,1);
    p[3].SetPoint(0,1,1);
    pl[0].ThreePointSetPlane(p[0], p[1], p[2]);
    pl[1].ThreePointSetPlane(p[0], p[1], p[3]);
    pl[2].ThreePointSetPlane(p[0], p[2], p[3]);
    std::cout<<PointPlaneDist(p[3],pl[0])<<"\n";
    std::cout<<PointPlaneDist(p[2],pl[1])<<"\n";
    std::cout<<PointPlaneDist(p[1],pl[2]);

#endif

#ifdef t1017
    class Point p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14;
    class Line l15, l16, l17, l18, l19, l20, l25, l26;
    class Line l15a;
    class Vector v21, v22, v23;
    //class Vector v21a;
    class Point  p21, p22, p23;
    p3.SetPoint(st[2].xr, st[2].ya, st[2].z);
    p4.SetPoint(st[3].xr, st[3].ya, st[3].z);
    p5.SetPoint(st[4].xr, st[4].ya, st[4].z);
    p6.SetPoint(st[5].xr, st[5].ya, st[5].z);
    p7.SetPoint(st[6].xr, st[6].ya, st[6].z);
    p8.SetPoint(st[7].xr, st[7].ya, st[7].z);
    p9.SetPoint(st[8].xr, st[8].ya, st[8].z);
    p10.SetPoint(st[9].xr, st[9].ya, st[9].z);
    p11.SetPoint(st[10].xr, st[10].ya, st[10].z);
    p12.SetPoint(st[11].xr, st[11].ya, st[11].z);
    p13.SetPoint(st[12].xr, st[12].ya, st[12].z);
    p14.SetPoint(st[13].xr, st[13].ya, st[13].z);
    l15 = LinkPoints(p4, p3);
    l15a= LinkPoints(p3, p4);
    l16 = LinkPoints(p5, p6);
    l17 = LinkPoints(p7, p8);
    l18 = LinkPoints(p9, p10);
    l19 = LinkPoints(p11, p12);
    l20 = LinkPoints(p13, p14);
    v21 = LineDistVec(l15a, l16,false);
    //v21a = LineDistVec(l15a, l16);
    v22 = LineDistVec(l17, l18);
    v23 = LineDistVec(l19, l20);
    //l15.vec.Out(true,true,true);
    //l15a.vec.Out(true,true,true);
    //v21.Out(true,true,true);
    //v21a.Out(true,true,true);
    p21.SetPoint((v21.op.x+v21.dp.x)/2,(v21.op.y+v21.dp.y)/2,(v21.op.z+v21.dp.z)/2);
    p22.SetPoint((v22.op.x+v22.dp.x)/2,(v22.op.y+v22.dp.y)/2,(v22.op.z+v22.dp.z)/2);
    p23.SetPoint((v23.op.x+v23.dp.x)/2,(v23.op.y+v23.dp.y)/2,(v23.op.z+v23.dp.z)/2);
    /*
    p21.Out();
    std::cout<<std::endl;
    p22.Out();
    std::cout<<std::endl;
    p23.Out();
    */
    v21.Out(true,false,true,true);
#endif

#ifdef t1118
    class Point p[12];
    class Circle c[40];
    class Line l[28];
    class Point pi[10];
    class Vector v[10];
    class PointsWDistr pwd[3];

    for (int i=0;i<12;i++) {
        p[i].SetPoint(st[i].xr,st[i].ya,st[i].z);
    }
    int ci = 0;
    for (int i=0;i<10;i++) {
        for(int j=0;j<4;j++) {
            c[ci].SetCircle(st[25+i*7+j].xr,st[25+i*7+j].ya,st[25+i*7+j].z,st[25+i*7+j].size);
            ci++;
        }
        l[6+2*i] = LinkPoints(c[i*4],c[i*4+1]);
        l[7+2*i] = LinkPoints(c[i*4+2],c[i*4+3]);
        v[i] = LineDistVec(l[6+2*i],l[7+2*i]);
        //pi[i].SetPoint((v[i].op.x+v[i].dp.x)/2,(v[i].op.y+v[i].dp.y)/2,(v[i].op.z+v[i].dp.z)/2);
        pi[i] = MidVec(v[i]);
        pi[i].Out();
        std::cout<<std::endl;
    }
    std::vector<class Point>pv1;
    std::vector<class Point>pv2;
    std::vector<class Point>pv3;
    pv1.push_back(pi[0]);
    pv1.push_back(pi[1]);
    pv2.push_back(pi[2]);
    pv2.push_back(pi[3]);
    pv1.push_back(pi[4]);
    pv2.push_back(pi[5]);
    pv3.push_back(pi[6]);
    pv3.push_back(pi[7]);
    pv2.push_back(pi[8]);
    pv1.push_back(pi[9]);
    pwd[0] = PointsAve(pv1);
    pwd[0].count = pv1.size();
    pwd[1] = PointsAve(pv2);
    pwd[1].count = pv2.size();
    pwd[2] = PointsAve(pv3);
    pwd[2].count = pv3.size();
    pwd[0].Out(true, true);
    pwd[1].Out(true, true);
    pwd[2].Out(true, true);
#endif

#ifdef KU1215
    class Point p[32];
    class Line l[16];
    class Vector v[8];
    class Point pi[8];
    class Line edge[8];
    for (int i=0;i<32;i++) {
        p[i].SetPoint(st[i].xr,st[i].ya,st[i].z);
        p[i].Out();
        std::cout<<std::endl;
    }
    for (int i=0;i<16;i++) {
        l[i] = LinkPoints(p[2*i], p[2*i+1]);
    }
    v[0] = LineDistVec(l[0],l[7]);
    pi[0] = MidVec(v[0]);
    v[1] = LineDistVec(l[1],l[2]);
    pi[1] = MidVec(v[1]);
    v[2] = LineDistVec(l[3],l[4]);
    pi[2] = MidVec(v[2]);
    v[3] = LineDistVec(l[5],l[6]);
    pi[3] = MidVec(v[3]);
    v[4] = LineDistVec(l[8],l[15]);
    pi[4] = MidVec(v[4]);
    v[5] = LineDistVec(l[9],l[10]);
    pi[5] = MidVec(v[5]);
    v[6] = LineDistVec(l[11],l[12]);
    pi[6] = MidVec(v[6]);
    v[7] = LineDistVec(l[13],l[14]);
    pi[7] = MidVec(v[7]);
    edge[0] = LinkPoints(pi[0],pi[1]);
    edge[1] = LinkPoints(pi[1],pi[2]);
    edge[2] = LinkPoints(pi[2],pi[3]);
    edge[3] = LinkPoints(pi[3],pi[0]);
    edge[4] = LinkPoints(pi[4],pi[5]);
    edge[5] = LinkPoints(pi[5],pi[6]);
    edge[6] = LinkPoints(pi[6],pi[7]);
    edge[7] = LinkPoints(pi[7],pi[4]);

    std::cout<<"Angles: "<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<VecsAngle(edge[i].vec,edge[(i+1)%4].vec,false)/90-1<<" ";
    }
    std::cout<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<VecsAngle(edge[i+4].vec,edge[(i+1)%4+4].vec,false)/90-1<<" ";
    }
    std::cout<<"\nDistance: "<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<PointLineDist(pi[i],edge[(i+1)%4])<<" ";
    }
    std::cout<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<PointLineDist(pi[i+4],edge[(i+1)%4+4])<<" ";
    }
#endif

#ifdef t0207
    class Point p[32];
    class Line l[24];
    class Vector v[8];
    class Point pi[8];
    class Line edge[8];
    for (int i=0;i<16;i++) {
        p[i].SetPoint(st[i].xr,st[i].ya,st[i].z);
        p[i].Out();
        std::cout<<std::endl;
    }
    for (int i=0;i<8;i++) {
        l[i] = LinkPoints(p[2*i], p[2*i+1]);
    }
    for (int i=16;i<32;i++) {
        p[i].SetPoint(st[i+16].xr,st[i+16].ya,st[i+16].z);
        p[i].Out();
        std::cout<<std::endl;
    }
    for (int i=8;i<16;i++) {
        l[i] = LinkPoints(p[2*i], p[2*i+1]);
    }

    v[0] = LineDistVec(l[0],l[7]);    pi[0] = MidVec(v[0]);
    v[1] = LineDistVec(l[1],l[2]);    pi[1] = MidVec(v[1]);
    v[2] = LineDistVec(l[3],l[4]);    pi[2] = MidVec(v[2]);
    v[3] = LineDistVec(l[5],l[6]);    pi[3] = MidVec(v[3]);
    v[4] = LineDistVec(l[8],l[15]);   pi[4] = MidVec(v[4]);
    v[5] = LineDistVec(l[9],l[10]);   pi[5] = MidVec(v[5]);
    v[6] = LineDistVec(l[11],l[12]);  pi[6] = MidVec(v[6]);
    v[7] = LineDistVec(l[13],l[14]);  pi[7] = MidVec(v[7]);

    edge[0] = LinkPoints(pi[0],pi[1]);
    edge[1] = LinkPoints(pi[1],pi[2]);
    edge[2] = LinkPoints(pi[2],pi[3]);
    edge[3] = LinkPoints(pi[3],pi[0]);
    edge[4] = LinkPoints(pi[4],pi[5]);
    edge[5] = LinkPoints(pi[5],pi[6]);
    edge[6] = LinkPoints(pi[6],pi[7]);
    edge[7] = LinkPoints(pi[7],pi[4]);

    std::cout<<"Angles: "<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<VecsAngle(edge[i].vec,edge[(i+1)%4].vec,false)/90-1<<" ";
    }
    std::cout<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<VecsAngle(edge[i+4].vec,edge[(i+1)%4+4].vec,false)/90-1<<" ";
    }
    std::cout<<"\nDistance: "<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<PointLineDist(pi[i],edge[(i+1)%4])<<" ";
    }
    std::cout<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<PointLineDist(pi[i+4],edge[(i+1)%4+4])<<" ";
    }
#endif

#ifdef t0207z
    class Point p[48];
    class Line l[36];
    class Vector v[12];
    class Point pi[12];
    class Line edge[12];
    class Point c[3];
    for (int i=0;i<16;i++) {
        int dif = 0;
        p[i].SetPoint(st[i+dif].xr,st[i+dif].ya,st[i+dif].z);
        p[i].Out();
        std::cout<<std::endl;
    }
    for (int i=0;i<8;i++) {
        l[i] = LinkPoints(p[2*i], p[2*i+1]);
    }
    for (int i=16;i<32;i++) {
        int dif = 19;
        p[i].SetPoint(st[i+dif].xr,st[i+dif].ya,st[i+dif].z);
        p[i].Out();
        std::cout<<std::endl;
    }
    for (int i=8;i<16;i++) {
        l[i] = LinkPoints(p[2*i], p[2*i+1]);
    }
    for (int i=32;i<48;i++) {
        int dif = 35;
        p[i].SetPoint(st[i+dif].xr,st[i+dif].ya,st[i+dif].z);
        p[i].Out();
        std::cout<<std::endl;
    }
    for (int i=16;i<24;i++) {
        l[i] = LinkPoints(p[2*i], p[2*i+1]);
    }

    v[0] = LineDistVec(l[0],l[7]);    pi[0] = MidVec(v[0]);
    v[1] = LineDistVec(l[1],l[2]);    pi[1] = MidVec(v[1]);
    v[2] = LineDistVec(l[3],l[4]);    pi[2] = MidVec(v[2]);
    v[3] = LineDistVec(l[5],l[6]);    pi[3] = MidVec(v[3]);
    v[4] = LineDistVec(l[8],l[15]);   pi[4] = MidVec(v[4]);
    v[5] = LineDistVec(l[9],l[10]);   pi[5] = MidVec(v[5]);
    v[6] = LineDistVec(l[11],l[12]);  pi[6] = MidVec(v[6]);
    v[7] = LineDistVec(l[13],l[14]);  pi[7] = MidVec(v[7]);
    v[8] = LineDistVec(l[16],l[23]);  pi[8] = MidVec(v[8]);
    v[9] = LineDistVec(l[17],l[18]);  pi[9] = MidVec(v[9]);
    v[10] = LineDistVec(l[19],l[20]); pi[10] = MidVec(v[10]);
    v[11] = LineDistVec(l[21],l[22]); pi[11] = MidVec(v[11]);

    edge[0] = LinkPoints(pi[0],pi[1]);
    edge[1] = LinkPoints(pi[1],pi[2]);
    edge[2] = LinkPoints(pi[2],pi[3]);
    edge[3] = LinkPoints(pi[3],pi[0]);
    edge[4] = LinkPoints(pi[4],pi[5]);
    edge[5] = LinkPoints(pi[5],pi[6]);
    edge[6] = LinkPoints(pi[6],pi[7]);
    edge[7] = LinkPoints(pi[7],pi[4]);
    edge[8] = LinkPoints(pi[8],pi[9]);
    edge[9] = LinkPoints(pi[9],pi[10]);
    edge[10] = LinkPoints(pi[10],pi[11]);
    edge[11] = LinkPoints(pi[11],pi[8]);


    std::vector<class Point> pv1;
    for (int i=0;i<4;i++){
        pv1.push_back(pi[i]);
    }
    c[0] = PointsAve(pv1);
    std::vector<class Point> pv2;
    for (int i=4;i<8;i++){
        pv2.push_back(pi[i]);
    }
    c[1] = PointsAve(pv2);
    std::vector<class Point> pv3;
    for (int i=8;i<12;i++){
        pv3.push_back(pi[i]);
    }
    c[2] = PointsAve(pv3);

    std::cout<<"Angles: "<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<VecsAngle(edge[i].vec,edge[(i+1)%4].vec,false)/90-1<<" ";
    }
    std::cout<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<VecsAngle(edge[i+4].vec,edge[(i+1)%4+4].vec,false)/90-1<<" ";
    }
    std::cout<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<VecsAngle(edge[i+8].vec,edge[(i+1)%4+8].vec,false)/90-1<<" ";
    }

    std::cout<<"\nDistance: "<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<PointLineDist(pi[i],edge[(i+1)%4])<<" ";
    }
    std::cout<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<PointLineDist(pi[i+4],edge[(i+1)%4+4])<<" ";
    }
    std::cout<<std::endl;
    for (int i=0;i<4;i++) {
        std::cout<<PointLineDist(pi[i+8],edge[(i+1)%4+8])<<" ";
    }

    std::cout<<"\nZ: "<<std::endl;
        std::cout<<c[0].z<<" "<<c[1].z<<' '<<c[2].z<<std::endl;
#endif

#ifdef t02151p
    class Point p[16];
    class PointsWDistr c;
    for (int i=0;i<16;i++) {
        int dif = 0;
        p[i].SetPoint(st[i+dif].xr,st[i+dif].ya,st[i+dif].z);
        p[i].Out();
        std::cout<<std::endl;
    }

    std::vector<class Point> pv;
    for (int i=0;i<16;i++){
        pv.push_back(p[i]);
    }
    c = PointsAve(pv);

    std::cout<<"\nZ: "<<std::endl;
    std::cout<<std::fixed<< std::setprecision(5)<<c.z<<" ± "<<c.zSD<<" mm"<<std::endl;
#endif

#ifdef t021555
    class Point p1[25], p2[25];
    class PointsWDistr c1, c2, c2x[5], c2y[5];

    for (int i=0;i<25;i++) {
        int dif = 0;
        p1[i].SetPoint(st[i+dif].xr,st[i+dif].ya,st[i+dif].z);
        p1[i].Out();
        std::cout<<std::endl;
    }
    std::vector<class Point> pv1;
    for (int i=0;i<25;i++){
        pv1.push_back(p1[i]);
    }
    c1 = PointsAve(pv1);

    for (int i=0;i<25;i++) {
        int dif = 25;
        p2[i].SetPoint(st[i+dif].xr,st[i+dif].ya,st[i+dif].z);
        p2[i].Out();
        std::cout<<std::endl;
    }
    std::vector<class Point> pv2;
    for (int i=0;i<25;i++){
        pv2.push_back(p2[i]);
    }
    c2 = PointsAve(pv2);

    std::vector<class Point> pv2x[5];
    std::vector<class Point> pv2y[5];
    for (int i=0; i<5; i++) {
        for (int j=0;j<5; j++) {
            pv2x[i].push_back(p2[i*5+j]);
        }
    }
    pv2y[0].push_back(p2[0]);   pv2y[0].push_back(p2[9]);   pv2y[0].push_back(p2[10]);  pv2y[0].push_back(p2[19]);  pv2y[0].push_back(p2[20]);
    pv2y[1].push_back(p2[1]);   pv2y[1].push_back(p2[8]);   pv2y[1].push_back(p2[11]);  pv2y[0].push_back(p2[18]);  pv2y[0].push_back(p2[21]);
    pv2y[2].push_back(p2[2]);   pv2y[2].push_back(p2[7]);   pv2y[2].push_back(p2[12]);  pv2y[0].push_back(p2[17]);  pv2y[0].push_back(p2[22]);
    pv2y[3].push_back(p2[3]);   pv2y[3].push_back(p2[6]);   pv2y[3].push_back(p2[13]);  pv2y[0].push_back(p2[16]);  pv2y[0].push_back(p2[23]);
    pv2y[4].push_back(p2[4]);   pv2y[4].push_back(p2[5]);   pv2y[4].push_back(p2[14]);  pv2y[0].push_back(p2[15]);  pv2y[0].push_back(p2[24]);


    for (int i=0; i<5; i++) {
        c2x[i] = PointsAve(pv2x[i]);
        c2y[i] = PointsAve(pv2y[i]);
    }

    std::cout<<"\nZ in x: "<<std::endl;
    for (int i=0; i<5; i++) {
        std::cout<<std::fixed<< std::setprecision(5)<<c2x[i].z<<" ± "<<c2x[i].zSD<<" mm"<<std::endl;
    }

    std::cout<<"\nZ in y: "<<std::endl;
    for (int i=0; i<5; i++) {
        std::cout<<std::fixed<< std::setprecision(5)<<c2y[i].z<<" ± "<<c2y[i].zSD<<" mm"<<std::endl;
    }

    std::cout<<"\nZ: "<<std::endl;
    std::cout<<std::fixed<< std::setprecision(5)<<c1.z<<" ± "<<c1.zSD<<" mm"<<std::endl;
    std::cout<<std::fixed<< std::setprecision(5)<<c2.z<<" ± "<<c2.zSD<<" mm"<<std::endl;
#endif

#ifdef t0227
    class Point pr[16], prc[16], p4e[14], p4c[25], p5e[14], p5c[25];
    double dprc[16], dp4e[14], dp4c[25], dp5e[14], dp5c[25];
    double dprcave, dp4eave, dp4cave, dp5eave, dp5cave;
    double dprcavesd, dp4eavesd, dp4cavesd, dp5eavesd, dp5cavesd;
    class PointsWDistr cr, c4, c5, c4x[5], c4y[5], c5x[5],c5y[5];
    class Line lr[8], l4[8], l5[8];
    class Point pir[4], p4r[4], p5r[4];
    class Plane plr, pl4, pl5;
    class Plane plrLLS, plrcLLS, pl4LLS, pl5LLS;
    class Vector vr1, vr2, vrp, v41, v42, v4p, v51, v52, v5p;
    std::vector<class Point> pvr, pvrc, pv4c, pv5c, pv4cx[5], pv4cy[5], pv5cx[5], pv5cy[5];

    fo<<"Test Run t0227.\n";

    for (int i=0;i<16;i++) {
        int ph = 0;
        pr[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
        //pr[i].Out();
        //pr[i].Out(fo);
        std::cout<<std::endl;
    }
    for (int i=0;i<16;i++) {
        int ph = 34;
        prc[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
        //prc[i].Out();
        //prc[i].Out(fo);
        std::cout<<std::endl;
    }
    for (int i=0;i<14;i++) {
        int ph = 50;
        p4e[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
        //p4e[i].Out();
        //p4e[i].Out(fo);
        std::cout<<std::endl;
    }
    for (int i=0;i<25;i++) {
        int ph = 80;
        p4c[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
        //p4c[i].Out();
        //p4c[i].Out(fo);
        std::cout<<std::endl;
    }
    for (int i=0;i<14;i++) {
        int ph = 105;
        p5e[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
        //p5e[i].Out();
        //p5e[i].Out(fo);
        std::cout<<std::endl;
    }
    for (int i=0;i<25;i++) {
        int ph = 135;
        p5c[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
        //p5c[i].Out();
        //p5c[i].Out(fo);
        std::cout<<std::endl;
    }

    for (int i=0;i<8;i++) {lr[i] = LinkPoints(pr[2*i], pr[2*i+1]);}
    for (int i=0;i<2;i++) {
        l4[0+i*4] = LinkPoints(p4e[0+7*i], p4e[1+7*i]);
        l4[1+i*4] = LinkPoints(p4e[1+7*i], p4e[2+7*i]);
        l4[2+i*4] = LinkPoints(p4e[3+7*i], p4e[4+7*i]);
        l4[3+i*4] = LinkPoints(p4e[5+7*i], p4e[6+7*i]);
    }
    for (int i=0;i<2;i++) {
        l5[0+i*4] = LinkPoints(p5e[0+7*i], p5e[1+7*i]);
        l5[1+i*4] = LinkPoints(p5e[1+7*i], p5e[2+7*i]);
        l5[2+i*4] = LinkPoints(p5e[3+7*i], p5e[4+7*i]);
        l5[3+i*4] = LinkPoints(p5e[5+7*i], p5e[6+7*i]);
    }

    for (int i=0;i<4;i++) {pir[i] = TwoLineIntersection(lr[i], lr[(i+1)%4]);}
    for (int i=0;i<4;i++) {p4r[i] = TwoLineIntersection(l4[i], l4[(i+1)%4]);}
    for (int i=0;i<4;i++) {p5r[i] = TwoLineIntersection(l5[i], l5[(i+1)%4]);}

    plr.ThreePointSetPlane(pir[0],pir[1],pir[2]);
    //std::cout<<"\nThe 4th corner is "<<PointPlaneDist(pir[3],plr)<<" mm away from the reference plane.\n";
    fo<<"\nThe 4th corner is "<<PointPlaneDist(pir[3],plr)<<" mm away from the reference plane.\n";
    pl4.ThreePointSetPlane(p4r[0],p4r[1],p4r[2]);
    //std::cout<<"\nThe 4th corner is "<<PointPlaneDist(p4r[3],pl4)<<" mm away from the #4 plane.\n";
    fo<<"\nThe 4th corner is "<<PointPlaneDist(p4r[3],pl4)<<" mm away from the #4 plane.\n";
    pl5.ThreePointSetPlane(p5r[0],p5r[1],p5r[2]);
    //std::cout<<"\nThe 4th corner is "<<PointPlaneDist(p5r[3],pl5)<<" mm away from the #5 plane.\n";
    fo<<"\nThe 4th corner is "<<PointPlaneDist(p5r[3],pl5)<<" mm away from the #5 plane.\n";

    for (int i=0;i<16;i++) {pvr.push_back(pr[i]);}
    for (int i=0;i<16;i++) {pvrc.push_back(prc[i]);}
    for (int i=0;i<25;i++) {pv4c.push_back(p4c[i]);}
    for (int i=0;i<25;i++) {pv5c.push_back(p5c[i]);}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {pv4cx[i].push_back(p4c[i*5+j]);}}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {pv4cy[i].push_back(p4c[i+j*5]);}}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {pv5cx[i].push_back(p5c[i*5+j]);}}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {pv5cy[i].push_back(p5c[i+j*5]);}}
    cr = PointsAve(pvrc);
    c4 = PointsAve(pv4c);
    c5 = PointsAve(pv5c);
    for (int i=0;i<5;i++) {c4x[i] = PointsAve(pv4cx[i]);}
    for (int i=0;i<5;i++) {c4y[i] = PointsAve(pv4cy[i]);}
    for (int i=0;i<5;i++) {c5x[i] = PointsAve(pv5cx[i]);}
    for (int i=0;i<5;i++) {c5y[i] = PointsAve(pv5cy[i]);}

    plrLLS.MultiPointSetPlaneSVD(pvr);
    plrcLLS.MultiPointSetPlaneSVD(pvrc);
    pl4LLS.MultiPointSetPlaneSVD(pv4c);
    pl5LLS.MultiPointSetPlaneSVD(pv5c);

    dprcave = 0;
    for (int i=0;i<16;i++) {
        dprc[i] = PointPlaneDist(prc[i], plr);
        dprcave = dprcave + dprc[i];
    }
    dprcave = dprcave / 16;
    for (int i=0;i<16;i++) {dprcavesd = dprcavesd + pow(fabs(dprcave - dprc[i]),2);}
    dprcavesd = sqrt(dprcavesd / 16);

    dp4eave = 0;
    for (int i=0;i<14;i++) {
        dp4e[i] = PointPlaneDist(p4e[i], plr);
        dp4eave = dp4eave + dp4e[i];
    }
    dp4eave = dp4eave / 14;
    for (int i=0;i<14;i++) {dp4eavesd = dp4eavesd + pow(fabs(dp4eave - dp4e[i]),2);}
    dp4eavesd = sqrt(dp4eavesd / 14);

    dp4cave = 0;
    for (int i=0;i<25;i++) {
        dp4c[i] = PointPlaneDist(p4c[i], plr);
        dp4cave = dp4cave + dp4c[i];
    }
    dp4cave = dp4cave / 25;
    for (int i=0;i<14;i++) {dp4cavesd = dp4cavesd + pow(fabs(dp4cave - dp4c[i]),2);}
    dp4cavesd = sqrt(dp4cavesd / 25);

    dp5eave = 0;
    for (int i=0;i<14;i++) {
        dp5e[i] = PointPlaneDist(p5e[i], plr);
        dp5eave = dp5eave + dp5e[i];
    }
    dp5eave = dp5eave / 14;
    for (int i=0;i<14;i++) {dp5eavesd = dp5eavesd + pow(fabs(dp5eave - dp5e[i]),2);}
    dp5eavesd = sqrt(dp5eavesd / 14);

    dp5cave = 0;
    for (int i=0;i<25;i++) {
        dp5c[i] = PointPlaneDist(p5c[i], plr);
        dp5cave = dp5cave + dp5c[i];
    }
    dp5cave = dp5cave / 25;
    for (int i=0;i<25;i++) {dp5cavesd = dp5cavesd + pow(fabs(dp5cave - dp5c[i]),2);}
    dp5cavesd = sqrt(dp5cavesd / 25);

    fo<<"\nZ in x for #4: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c4x[i],plr)<<" ± "<<c4x[i].zSD<<" mm"<<std::endl;}
    fo<<"\nZ in y for #4: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c4y[i],plr)<<" ± "<<c4y[i].zSD<<" mm"<<std::endl;}
    fo<<"\nZ in x for #5: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c5x[i],plr)<<" ± "<<c5x[i].zSD<<" mm"<<std::endl;}
    fo<<"\nZ in y for #5: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c5y[i],plr)<<" ± "<<c5y[i].zSD<<" mm"<<std::endl;}

    /*
    std::cout<<"\nZ in x for #4: "<<std::endl;
    for (int i=0; i<5; i++) {std::cout<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c4x[i],plr)<<" ± "<<c4x[i].zSD<<" mm"<<std::endl;}
    std::cout<<"\nZ in y for #4: "<<std::endl;
    for (int i=0; i<5; i++) {std::cout<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c4y[i],plr)<<" ± "<<c4y[i].zSD<<" mm"<<std::endl;}
    std::cout<<"\nZ in x for #5: "<<std::endl;
    for (int i=0; i<5; i++) {std::cout<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c5x[i],plr)<<" ± "<<c5x[i].zSD<<" mm"<<std::endl;}
    std::cout<<"\nZ in y for #5: "<<std::endl;
    for (int i=0; i<5; i++) {std::cout<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c5y[i],plr)<<" ± "<<c5y[i].zSD<<" mm"<<std::endl;}
     */

    fo<<"\nZ on reference plane (Ave. Point): "<<std::endl;
    fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(cr,plr)<<" ± "<<cr.zSD<<" mm"<<std::endl;
    fo<<"\nZ on reference plane (Ave. Dist.): "<<std::endl;
    fo<<std::fixed<< std::setprecision(5)<<dprcave<<" ± "<<dprcavesd<<" mm"<<std::endl;
    fo<<"\nZ on #4 (Ave. Point): "<<std::endl;
    fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c4,plr)<<" ± "<<c4.zSD<<" mm"<<std::endl;
    fo<<"\nZ on #4 (Ave. Dist.): "<<std::endl;
    fo<<std::fixed<< std::setprecision(5)<<dp4cave<<" ± "<<dp4cavesd<<" mm"<<std::endl;
    fo<<"\nZ on #5 (Ave. Point): "<<std::endl;
    fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c5,plr)<<" ± "<<c5.zSD<<" mm"<<std::endl;
    fo<<"\nZ on #4 (Ave. Dist.): "<<std::endl;
    fo<<std::fixed<< std::setprecision(5)<<dp5cave<<" ± "<<dp5cavesd<<" mm"<<std::endl;
    fo<<std::endl;

    plrLLS.op.Out(fo); plrLLS.normvec.Out(fo); fo<<" and the standard deviation is "<<MultiPointsPlaneSD(pvr,plrLLS)<<" mm.\n";
    plrcLLS.op.Out(fo);plrcLLS.normvec.Out(fo);fo<<" and the standard deviation is "<<MultiPointsPlaneSD(pvrc,plrcLLS)<<" mm.\n";
    pl4LLS.op.Out(fo); pl4LLS.normvec.Out(fo); fo<<" and the standard deviation is "<<MultiPointsPlaneSD(pv4c,pl4LLS)<<" mm.\n";
    pl5LLS.op.Out(fo); pl5LLS.normvec.Out(fo); fo<<" and the standard deviation is "<<MultiPointsPlaneSD(pv5c,pl5LLS)<<" mm.\n";
    fo<<"\n In reference to plane,  #4 is "<<pl4LLS.op.z-plrLLS.op.z<<" mm and #5 is "<<pl5LLS.op.z-plrLLS.op.z<<" mm.\n";
    fo<<"\n In reference to center, #4 is "<<pl4LLS.op.z-plrcLLS.op.z<<" mm and #5 is "<<pl5LLS.op.z-plrcLLS.op.z<<" mm.\n";
#endif

#ifdef t0317probe
    class Point pr[25], p4[25], p5[25];
    double dp4[25], dp5[25];
    double dprcave, dp4eave, dp4cave, dp5eave, dp5cave;
    double dprcavesd, dp4eavesd, dp4cavesd, dp5eavesd, dp5cavesd;
    class PointsWDistr cr, c4, c5, crx[5], cry[5], c4x[5], c4y[5], c5x[5], c5y[5];
    class Line lr[8], l4[8], l5[8];
    class Point pir[4], p4r[4], p5r[4];
    class Plane plr, pl4, pl5;
    class Vector vr, v4, v5;
    std::vector<class Point> pvr, pv4, pv5, vrx[5], vry[5], v4x[5], v4y[5], v5x[5], v5y[5];

    fo<<"Test Run t0317probe.\n";

    for (int i=0;i<25;i++) {
        int ph = 0;
        pr[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
    }
    for (int i=0;i<25;i++) {
        int ph = 25;
        p4[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
    }
    for (int i=0;i<25;i++) {
        int ph = 50;
        p5[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);
    }
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {vrx[i].push_back(pr[i*5+j]);}}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {vry[i].push_back(pr[i+j*5]);}}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {v4x[i].push_back(p4[i*5+j]);}}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {v4y[i].push_back(p4[i+j*5]);}}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {v5x[i].push_back(p5[i*5+j]);}}
    for (int i=0;i<5;i++) {for (int j=0;j<5;j++) {v5y[i].push_back(p5[i+j*5]);}}
    for (int i=0;i<5;i++) {crx[i] = PointsAve(vrx[i]);}
    for (int i=0;i<5;i++) {cry[i] = PointsAve(vry[i]);}
    for (int i=0;i<5;i++) {c4x[i] = PointsAve(v4x[i]);}
    for (int i=0;i<5;i++) {c4y[i] = PointsAve(v4y[i]);}
    for (int i=0;i<5;i++) {c5x[i] = PointsAve(v5x[i]);}
    for (int i=0;i<5;i++) {c5y[i] = PointsAve(v5y[i]);}

    for (int i=0;i<25;i++) {pvr.push_back(pr[i]);}
    for (int i=0;i<25;i++) {pv4.push_back(p4[i]);}
    for (int i=0;i<25;i++) {pv5.push_back(p5[i]);}

    plr.MultiPointSetPlaneSVD(pvr);
    pl4.MultiPointSetPlaneSVD(pv4);
    pl5.MultiPointSetPlaneSVD(pv5);

    fo<<"\nZ in x for #30: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(crx[i],plr)<<" ± "<<crx[i].zSD<<" mm"<<std::endl;}
    fo<<"\nZ in y for #30: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(cry[i],plr)<<" ± "<<cry[i].zSD<<" mm"<<std::endl;}
    fo<<"\nZ in x for #4: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c4x[i],plr)<<" ± "<<c4x[i].zSD<<" mm"<<std::endl;}
    fo<<"\nZ in y for #4: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c4y[i],plr)<<" ± "<<c4y[i].zSD<<" mm"<<std::endl;}
    fo<<"\nZ in x for #5: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c5x[i],plr)<<" ± "<<c5x[i].zSD<<" mm"<<std::endl;}
    fo<<"\nZ in y for #5: "<<std::endl;
    for (int i=0; i<5; i++) {fo<<std::fixed<< std::setprecision(5)<<PointPlaneDist(c5y[i],plr)<<" ± "<<c5y[i].zSD<<" mm"<<std::endl;}
    fo<<std::endl;

    plr.op.Out(fo);
    plr.normvec.Out(fo);
    fo<<" and the standard deviation is "<<MultiPointsPlaneSD(pvr,plr)<<" mm.\n";
    pl4.op.Out(fo);
    pl4.normvec.Out(fo);
    fo<<" and the standard deviation is "<<MultiPointsPlaneSD(pv4,pl4)<<" mm.\n";
    pl5.op.Out(fo);
    pl5.normvec.Out(fo);
    fo<<" and the standard deviation is "<<MultiPointsPlaneSD(pv5,pl5)<<" mm.\n";
    fo<<"#4 is "<<pl4.op.z-plr.op.z<<" mm and #5 is "<<pl5.op.z-plr.op.z<<" mm.\n";

    /*
    plr.op.Out();
    plr.normvec.Out();
    std::cout<<" and the standard deviation is "<<MultiPointsPlaneSD(pvr,plr)<<" mm.\n";
    pl4.op.Out();
    pl4.normvec.Out();
    std::cout<<" and the standard deviation is "<<MultiPointsPlaneSD(pv4,pl4)<<" mm.\n";
    pl5.op.Out();
    pl5.normvec.Out();
    std::cout<<" and the standard deviation is "<<MultiPointsPlaneSD(pv5,pl5)<<" mm.\n";
    std::cout<<"#4 is "<<pl4.op.z-plr.op.z<<" mm and #5 is "<<pl5.op.z-plr.op.z<<" mm.\n";
     */
#endif
#ifdef UU0328_1
#endif
#ifdef UU0328_2
#endif
#ifdef UU0329CrossF
#endif
#ifdef UU0329ChequeF
#endif
#ifdef UU0329FF
#endif
#ifdef UU0329Zring168
#endif
#ifdef UU0329ZTTL162
#endif
#ifdef UU0329ZTTL180
#endif
#if defined(UU0330_1) || defined(UU0330_2) || defined(UU0330_3)
    class Point CS[16], H_R1H0_P1[8], H_R1H0_P2[8],
                        H_R1H1_P1[8], H_R1H1_P2[8],
                        PB_P1[4], PB_P2[4], SiRef[23],
                        H_R1H0_0[4], H_R1H0_1[4], H_R1H0_2[4], HCC_R1H0_4[4], ABC_R1H0[10][4],
                        H_R1H1_0[4], H_R1H1_1[4], H_R1H1_2[4], HCC_R1H1_5[4], ABC_R1H1[11][4],
                        PB_1[4], PB_2, PB_3[4], PB_4[4], PB_5[4],
                        H_R1H0_P1_ave, H_R1H0_P2_ave, H_R1H1_P1_ave, H_R1H1_P2_ave, PB_P1_ave, PB_P2_ave;
    std::vector<class Point> H_R1H0_0_v, H_R1H0_1_v, H_R1H0_2_v, HCC_R1H0_4_v, ABC_R1H0_v[10],
                             H_R1H1_0_v, H_R1H1_1_v, H_R1H1_2_v, HCC_R1H1_5_v, ABC_R1H1_v[11],
                             PB_1_v, PB_3_v, PB_4_v, PB_5_v;
    class PointsWDistr H_R1H0_0_ave, H_R1H0_1_ave, H_R1H0_2_ave, HCC_R1H0_4_ave, ABC_R1H0_ave[10],
                       H_R1H1_0_ave, H_R1H1_1_ave, H_R1H1_2_ave, HCC_R1H1_5_ave, ABC_R1H1_ave[11],
                       PB_1_ave, PB_3_ave, PB_4_ave, PB_5_ave;
    class Point corner[4];
    class Line CSl[8], Edge[4];
    class Plane SiRefP;
    std::vector<class Point> SiRefV;
    double SiRefpSD;

    for (int i=0;i<16;i++) {int ph = 0; CS[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<8 ;i++) {int ph = 34; H_R1H0_P1[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<8 ;i++) {int ph = 53; H_R1H0_P2[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<8 ;i++) {int ph = 72; H_R1H1_P1[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<8 ;i++) {int ph = 91; H_R1H1_P2[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 110; PB_P1[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 117; PB_P2[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<23;i++) {int ph = 124; SiRef[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 147; H_R1H0_0[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 151; H_R1H0_1[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 155; H_R1H0_2[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int num=0;num<10;num++) {for (int i=0;i<4;i++) {int ph = 159+num*4; ABC_R1H0[num][i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}}
    for (int i=0;i<4 ;i++) {int ph = 199; HCC_R1H0_4[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 203; H_R1H1_0[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 207; H_R1H1_1[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 211; H_R1H1_2[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int num=0;num<11;num++) {for (int i=0;i<4;i++) {int ph = 215+num*4; ABC_R1H1[num][i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}}
    for (int i=0;i<4 ;i++) {int ph = 259; HCC_R1H1_5[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 263; PB_1[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    PB_2.SetPoint(st[267].xr,st[267].ya,st[267].z);
    for (int i=0;i<4 ;i++) {int ph = 268; PB_3[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 272; PB_4[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4 ;i++) {int ph = 275; PB_5[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}

    for (int i=0;i<8;i++) {CSl[i] = LinkPoints(CS[2*i], CS[2*i+1]);}
    for (int i=0;i<4;i++) {corner[i] = TwoLineIntersection(CSl[2*i], CSl[2*i+1]);}
    Edge[0] = LinkPoints(corner[0], corner[1]);
    Edge[1] = LinkPoints(corner[1], corner[2]);
    Edge[2] = LinkPoints(corner[2], corner[3]);
    Edge[3] = LinkPoints(corner[3], corner[0]);

    H_R1H0_P1_ave = CrossFiducial(H_R1H0_P1);
    H_R1H0_P2_ave = CrossFiducial(H_R1H0_P2);
    H_R1H1_P1_ave = CrossFiducial(H_R1H1_P1);
    H_R1H1_P2_ave = CrossFiducial(H_R1H1_P2);
    PB_P1_ave = SquareFiducial(PB_P1);
    PB_P2_ave = SquareFiducial(PB_P2);

    for (const auto & i : SiRef) {SiRefV.push_back(i);}
    SiRefP.MultiPointSetPlaneSVD(SiRefV);
    SiRefpSD = MultiPointsPlaneSD(SiRefV,SiRefP);

    //for (const auto & i : H_R1H0_0) {H_R1H0_0_v.push_back(i);} H_R1H0_0_ave = PointsAve(H_R1H0_0_v);
    for (int i=0;i<4;i++) {H_R1H0_0_v.push_back(H_R1H0_0[i]);} H_R1H0_0_ave = PointsAve(H_R1H0_0_v);
    for (int i=0;i<4;i++) {H_R1H0_1_v.push_back(H_R1H0_1[i]);} H_R1H0_1_ave = PointsAve(H_R1H0_1_v);
    for (int i=0;i<4;i++) {H_R1H0_2_v.push_back(H_R1H0_2[i]);} H_R1H0_2_ave = PointsAve(H_R1H0_2_v);
    for (int num=0;num<10;num++) {for (int i=0;i<4;i++) {ABC_R1H0_v[num].push_back(ABC_R1H0[num][i]);} ABC_R1H0_ave[num] = PointsAve(ABC_R1H0_v[num]);}
    for (int i=0;i<4;i++) {HCC_R1H0_4_v.push_back(HCC_R1H0_4[i]);} HCC_R1H0_4_ave = PointsAve(HCC_R1H0_4_v);
    for (int i=0;i<4;i++) {H_R1H1_0_v.push_back(H_R1H1_0[i]);} H_R1H1_0_ave = PointsAve(H_R1H1_0_v);
    for (int i=0;i<4;i++) {H_R1H1_1_v.push_back(H_R1H1_1[i]);} H_R1H1_1_ave = PointsAve(H_R1H1_1_v);
    for (int i=0;i<4;i++) {H_R1H1_2_v.push_back(H_R1H1_2[i]);} H_R1H1_2_ave = PointsAve(H_R1H1_2_v);
    for (int num=0;num<11;num++) {for (int i=0;i<4;i++) {ABC_R1H1_v[num].push_back(ABC_R1H1[num][i]);} ABC_R1H1_ave[num] = PointsAve(ABC_R1H1_v[num]);}
    for (int i=0;i<4;i++) {HCC_R1H1_5_v.push_back(HCC_R1H1_5[i]);} HCC_R1H1_5_ave = PointsAve(HCC_R1H1_5_v);
    for (int i=0;i<4;i++) {PB_1_v.push_back(PB_1[i]);} PB_1_ave = PointsAve(PB_1_v);
    for (int i=0;i<4;i++) {PB_3_v.push_back(PB_3[i]);} PB_3_ave = PointsAve(PB_3_v);
    for (int i=0;i<4;i++) {PB_4_v.push_back(PB_4[i]);} PB_4_ave = PointsAve(PB_4_v);
    for (int i=0;i<4;i++) {PB_5_v.push_back(PB_5[i]);} PB_5_ave = PointsAve(PB_5_v);

    for (int num=0;num<10;num++) {ABC_1_v_t.push_back(ABC_R1H0_v[num]);}
    for (int num=0;num<11;num++) {ABC_2_v_t.push_back(ABC_R1H1_v[num]);}
    HCC_1_v_t.push_back(H_R1H0_0_v);
    HCC_1_v_t.push_back(H_R1H0_1_v);
    HCC_1_v_t.push_back(H_R1H0_2_v);
    HCC_1_v_t.push_back(HCC_R1H0_4_v);
    HCC_2_v_t.push_back(H_R1H1_0_v);
    HCC_2_v_t.push_back(H_R1H1_1_v);
    HCC_2_v_t.push_back(H_R1H1_2_v);
    HCC_2_v_t.push_back(HCC_R1H1_5_v);
    HyXY.push_back(H_R1H0_P1_ave);
    HyXY.push_back(H_R1H0_P2_ave);
    HyXY.push_back(H_R1H1_P1_ave);
    HyXY.push_back(H_R1H1_P2_ave);
    PB1 = PB_P1_ave;
    PB2 = PB_P2_ave;
    Si = SiRefV;
    ABC.push_back(ABC_1_v_t);
    ABC.push_back(ABC_2_v_t);
    HCC.push_back(HCC_1_v_t);
    HCC.push_back(HCC_2_v_t);
    PB.push_back(PB_1_ave);
    PB.push_back(PB_2);
    PB.push_back(PB_3_ave);
    PB.push_back(PB_4_ave);
    PB.push_back(PB_5_ave);
    file_output_head(fo, 0, "M1", "N/A", "Uppsala Universitet", "Xiangyu XU", "OGP SmartScope Flash 200", 3, "v.23.05.02");
    file_output_body(fo, 1, HyXY, PB_P1_ave, PB_P2_ave, SiRefV, HCC, ABC, PB); //Where's H_RaHb_0 & HCC_RaHb_4?

    fo<<"\nH_R1H0_0 is "<<PointPlaneDist(H_R1H0_0_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"H_R1H0_1 is "<<PointPlaneDist(H_R1H0_1_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"H_R1H0_2 is "<<PointPlaneDist(H_R1H0_2_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"HCC_R1H0_4 is "<<PointPlaneDist(HCC_R1H0_4_ave,SiRefP)<<" mm above the reference plane.\n";
    for (int i=0;i<10;i++) {fo<<"ABC_R1H0_"<<i<<" is "<<PointPlaneDist(ABC_R1H0_ave[i],SiRefP)<<" mm above the reference plane.\n";}
    fo<<"H_R1H1_0 is "<<PointPlaneDist(H_R1H1_0_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"H_R1H1_1 is "<<PointPlaneDist(H_R1H1_1_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"H_R1H1_2 is "<<PointPlaneDist(H_R1H1_2_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"HCC_R1H1_5 is "<<PointPlaneDist(HCC_R1H1_5_ave,SiRefP)<<" mm above the reference plane.\n";
    for (int i=0;i<11;i++) {fo<<"ABC_R1H1_"<<i<<" is "<<PointPlaneDist(ABC_R1H1_ave[i],SiRefP)<<" mm above the reference plane.\n";}
    fo<<"PB_1 is "<<PointPlaneDist(PB_1_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"PB_2 is "<<PointPlaneDist(PB_2,SiRefP)<<" mm above the reference plane.\n";
    fo<<"PB_3 is "<<PointPlaneDist(PB_3_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"PB_4 is "<<PointPlaneDist(PB_4_ave,SiRefP)<<" mm above the reference plane.\n";
    fo<<"PB_5 is "<<PointPlaneDist(PB_5_ave,SiRefP)<<" mm above the reference plane.\n";

    fo<<"\n H_R1H0_P1/2, H_R1H1_P1/2, PB_P1/2\n";
    H_R1H0_P1_ave.Out(fo);
    fo<<std::endl;
    H_R1H0_P2_ave.Out(fo);
    fo<<std::endl;
    H_R1H1_P1_ave.Out(fo);
    fo<<std::endl;
    H_R1H1_P2_ave.Out(fo);
    fo<<std::endl;
    PB_P1_ave.Out(fo);
    fo<<std::endl;
    PB_P2_ave.Out(fo);

    fo<<"\n\nH_R1H0_0/1/2/4, ABC_R1H0_0-10, H_R1H1_0/1/2/5, ABC_R1H1_0-11, PB_1-5\n";
    fo<<PointPlaneDist(H_R1H0_0_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(H_R1H0_1_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(H_R1H0_2_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(HCC_R1H0_4_ave,SiRefP)<<"\n";
    for (int i=0;i<10;i++) {fo<<PointPlaneDist(ABC_R1H0_ave[i],SiRefP)<<"\n";}
    fo<<PointPlaneDist(H_R1H1_0_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(H_R1H1_1_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(H_R1H1_2_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(HCC_R1H1_5_ave,SiRefP)<<"\n";
    for (int i=0;i<11;i++) {fo<<PointPlaneDist(ABC_R1H1_ave[i],SiRefP)<<"\n";}
    fo<<PointPlaneDist(PB_1_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(PB_2,SiRefP)<<"\n";
    fo<<PointPlaneDist(PB_3_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(PB_4_ave,SiRefP)<<"\n";
    fo<<PointPlaneDist(PB_5_ave,SiRefP)<<"\n";

#endif

#if defined(LU0502_1_Op) || defined(LU0502_1_TP)
    class Point pl1p[12], pl2p[8];
    class Plane pl1pp, pl2pp;
    double pl1pSD, pl2pSD;
    std::vector<class Point> pl1pv, pl2pv;
    for (int i=0;i<12;i++) {int ph = 0; pl1p[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<8;i++) {int ph = 12; pl2p[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (const auto & i : pl1p) {pl1pv.push_back(i);}
    pl1pp.MultiPointSetPlaneSVD(pl1pv);
    pl1pSD = MultiPointsPlaneSD(pl1pv,pl1pp);
    for (const auto & i : pl2p) {pl2pv.push_back(i);}
    pl2pp.MultiPointSetPlaneSVD(pl2pv);
    pl2pSD = MultiPointsPlaneSD(pl2pv,pl2pp);
#ifdef LU0502_1_Op
    fo<<"Item 1, Optical.\n";
#endif
#ifdef LU0502_1_TP
    fo<<"Item 1, Touch Probe.\n";
#endif
    pl1pp.op.Out(fo);
    pl1pp.normvec.Out(fo);
    fo<<" and the standard deviation is "<<pl1pSD<<" mm.\n";
    pl2pp.op.Out(fo);
    pl2pp.normvec.Out(fo);
    fo<<" and the standard deviation is "<<pl2pSD<<" mm.\n";
    fo<<"The height difference is "<<pl2pp.op.z-pl1pp.op.z<<" mm.\n";
#endif
#if defined(LU0502_2_Op) || defined(LU0502_2_TP)
    class Point pl1p[8], pl2p[8];
    class Plane pl1pp, pl2pp;
    double pl1pSD, pl2pSD;
    std::vector<class Point> pl1pv, pl2pv;
    for (int i=0;i<8;i++) {int ph = 0; pl1p[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<8;i++) {int ph = 8; pl2p[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (const auto & i : pl1p) {pl1pv.push_back(i);}
    pl1pp.MultiPointSetPlaneSVD(pl1pv);
    pl1pSD = MultiPointsPlaneSD(pl1pv,pl1pp);
    for (const auto & i : pl2p) {pl2pv.push_back(i);}
    pl2pp.MultiPointSetPlaneSVD(pl2pv);
    pl2pSD = MultiPointsPlaneSD(pl2pv,pl2pp);
#ifdef LU0502_2_Op
    fo<<"Item 2, Optical.\n";
#endif
#ifdef LU0502_2_TP
    fo<<"Item 2, Touch Probe.\n";
#endif
    pl1pp.op.Out(fo);
    pl1pp.normvec.Out(fo);
    fo<<" and the standard deviation is "<<pl1pSD<<" mm.\n";
    pl2pp.op.Out(fo);
    pl2pp.normvec.Out(fo);
    fo<<" and the standard deviation is "<<pl2pSD<<" mm.\n";
    fo<<"The height difference is "<<pl2pp.op.z-pl1pp.op.z<<" mm.\n";
#endif
#if defined(LU0502_3_Op) || defined(LU0502_3_TP)
    class Point pl1p[12], pl2p[4], pl3p[8];
    class Plane pl1pp, pl2pp, pl3pp;
    double pl1pSD, pl2pSD, pl3pSD;
    class PointsWDistr pl2p_ave;
    std::vector<class Point> pl1pv, pl2pv, pl3pv;
    for (int i=0;i<12;i++) {int ph = 0; pl1p[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<4;i++) {int ph = 12; pl2p[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (int i=0;i<8;i++) {int ph = 16; pl3p[i].SetPoint(st[i+ph].xr,st[i+ph].ya,st[i+ph].z);}
    for (const auto & i : pl1p) {pl1pv.push_back(i);}
    pl1pp.MultiPointSetPlaneSVD(pl1pv);
    pl1pSD = MultiPointsPlaneSD(pl1pv,pl1pp);
    for (const auto & i : pl2p) {pl2pv.push_back(i);}
    pl2p_ave = PointsAve(pl2pv);
    //pl2pp.MultiPointSetPlaneSVD(pl2pv);
    //pl2pSD = MultiPointsPlaneSD(pl2pv,pl2pp);
    for (const auto & i : pl3p) {pl3pv.push_back(i);}
    pl3pp.MultiPointSetPlaneSVD(pl3pv);
    pl3pSD = MultiPointsPlaneSD(pl3pv,pl3pp);
#ifdef LU0502_3_Op
    fo<<"Item 3, Optical.\n";
#endif
#ifdef LU0502_3_TP
    fo<<"Item 3, Touch Probe.\n";
#endif
    pl1pp.op.Out(fo);
    pl1pp.normvec.Out(fo);
    fo<<" and the standard deviation is "<<pl1pSD<<" mm.\n";
    //pl2pp.op.Out(fo);
    //pl2pp.normvec.Out(fo);
    //fo<<" and the standard deviation is "<<pl2pSD<<" mm.\n";
    pl3pp.op.Out(fo);
    pl3pp.normvec.Out(fo);
    fo<<" and the standard deviation is "<<pl3pSD<<" mm.\n";
    fo<<"The height difference is "<<PointPlaneDist(pl2p_ave, pl1pp)<<" mm.\n";
    fo<<"The height difference is "<<pl3pp.op.z-pl1pp.op.z<<" mm.\n";
#endif

#if defined(debug)
    file_output_head(fo, 0, "Not-a-module", "N/A", "This University", "John Smith", "OGP SmartScope", 1, "v.");
    file_output_body(fo, 1, HyXY, PB1, PB2, Si, ABC, PB);
#endif

    fo.close();
    std::cout << "Output has been written to " << work_dir + file_dir + out_folder + out_name << " successfully.\n";
    //test_path(work_dir + file_dir + out_folder, out_name, false);
    return 0;
}
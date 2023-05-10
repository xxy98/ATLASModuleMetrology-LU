//
// Created by Xiangyu on 2022/11/07.
//

#ifndef ANA_FUNCTIONS_H
#define ANA_FUNCTIONS_H

#endif //ANA_FUNCTIONS_H

#ifndef ANA_MATHEMATICS_H
#define ANA_MATHEMATICS_H
#include "mathematics.h"
#endif

#ifndef ch
#include <chrono>
#endif
#ifndef ct
#include <ctime>
#endif

#ifndef HHdate
#include "date/date.h"
#include "date/.z.h" //These third-party heads are used for proper timezone outputting. Otherwise you will se ".z)" in put_time() gives some random rubbish based on your lang. Half a day spent for that.
#include "date/tz_private.h"
#include "src/.z.cpp"
#endif
#ifndef stdio
#include "stdio.h"
#endif

template <class T>
int ArrayLen(T& array) {
    return (sizeof(array) / sizeof(array[0]));
}

void test_path(std::string path, std::string name, bool output) {
    std::string file_path = path + name;
    std::cout << file_path << std::endl;
    std::ifstream f(file_path, std::ios::in);
    if (!f) {
        f.close();
        std::cout << "Error opening source file." << std::endl;
        return;
    }
    std::string line;
    //std::cout << "Successful read of " << file_path << " ." << std::endl;
    while (getline(f, line)){
        if (output) {
#ifdef debug
            std::cout << line << std::endl;
#endif
        }
    }
    f.close();
}

int get_line_num(std::string path, std::string name) {
    std::string file_path = path + name;
    std::cout << file_path << std::endl;
    std::ifstream f(file_path, std::ios::in);
    if (!f) {
        f.close();
        std::cout << "Error opening source file." << std::endl;
        return -1;
    }
    std::string line;
    int num = 0;
    while (getline(f, line)){
        num++;
    }
    f.close();
#ifdef LU
    num = num - 3;
#endif
#ifdef UU
    num = num - 4;
#endif
#ifdef KU
    num = num - 2;
#endif
    std::cout << "This file has " << num << " steps.\n";
    return num;
}

//Chao De
constexpr size_t HASH_STRING_PIECE(const char *string_piece,size_t hashNum=0){
    return *string_piece?HASH_STRING_PIECE(string_piece+1,(hashNum*131)+*string_piece):hashNum;
}

//Chao De
constexpr size_t operator "" _HASH(const char *string_pice,size_t){
    return HASH_STRING_PIECE(string_pice);
}

Measurement Meas_Input(std::string line){
    Measurement line_m;

    int st_i;
    int ref_i[2];
    Feature ft_e;
    double xr_d, ya_d, z_d, sz_d;
    std::string st_s, ft_s, xr_s, ya_s, z_s, sz_s, ref_s;
    int count = 0;
#ifdef LU
    static std::string blank = "          ";
    st_s  = line.substr(0,4);
    ft_s  = line.substr(6,16);
    xr_s  = line.substr(22,10);
    ya_s  = line.substr(36,10);
    z_s   = line.substr(48,10);
    sz_s  = line.substr(60,10);
    ref_s = line.substr(70,10);
#endif
#ifdef UU
    static std::string blank = "         ";
    st_s  = line.substr(0,4);
    ft_s  = line.substr(6,16);
    xr_s  = line.substr(22,9);
    ya_s  = line.substr(35,9);
    z_s   = line.substr(47,9);
    sz_s  = line.substr(59,9);
    ref_s = line.substr(70,10);
#endif
#ifdef KU
    static std::string blank = "          ";
    st_s  = line.substr(0,4);
    ft_s  = line.substr(5,16);
    xr_s  = line.substr(20,11);
    ya_s  = line.substr(34,11);
    z_s   = line.substr(46,11);
    sz_s  = line.substr(59,10);
    ref_s = line.substr(69,10);
#endif
    st_i = atoi(st_s.c_str());
    switch (HASH_STRING_PIECE(ft_s.c_str())) {
        case "Plane           "_HASH:
            ft_e = Plane;
            break;
        case "Point           "_HASH:
            ft_e = Point;
            break;
        case "Line            "_HASH:
            ft_e = Line;
            break;
        case "Intersect       "_HASH:
            ft_e = Intersect;
            break;
        case "Circle          "_HASH:
            ft_e = Circle;
            break;
        case "Datum Origin    "_HASH:
            ft_e = DatumOrigin;
            break;
        case "Datum Plane     "_HASH:
            ft_e = DatumPlane;
            break;
    }
    xr_d = (xr_s == blank ? 0 : std::stod(xr_s));
    ya_d = (ya_s == blank ? 0 : std::stod(ya_s));
    z_d  = (z_s  == blank ? 0 : std::stod(z_s));
    sz_d = (sz_s == blank ? 0 : std::stod(sz_s));
    int ref_len = ref_s.length();
    int num = -1;
    bool space = true;
    std::string ref_ss;
    ref_i[0] = 0;
    ref_i[1] = 0;
    for (int i = 0; i < ref_len; i++) {
        char ref_c = ref_s.at(i);
        if (ref_c==' '){
            if (!space){
                ref_i[num] = atoi(ref_ss.c_str());
                ref_ss = "";
            }
            space = true;
        }
        else{
            if (space) {
                num++;
            }
            ref_ss += ref_c;
            space = false;
        }
    }
    line_m.Set(st_i, ft_e, xr_d, ya_d, z_d, sz_d, ref_i[0], ref_i[1]);
    return line_m;
}

std::string ISOtime() {
    auto now = std::chrono::system_clock::now();
    auto ms = duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
    auto us = duration_cast<std::chrono::microseconds>(now.time_since_epoch()) % 1000;
    auto timer = std::chrono::system_clock::to_time_t(now);
    std::tm bt = *std::localtime(&timer);

    auto t = make_zoned(date::current_zone(), std::chrono::system_clock::now());

    std::ostringstream oss;
//    oss << std::put_time(&bt, "%FT%T");               //This gives NOTHING.
    oss << std::put_time(&bt, "%Y-%m-%dT%H:%M:%S");
//    oss << format("%FT%T",t);                         //This somehow will go to millisecond, microsecond and ... to nanosecond! I DON'T KNOW WHY!
//    oss << format("%FT%H:%M:%S",t);
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count() << std::setw(3) << us.count();
    oss << format("%z",t);
    return oss.str();
}

std::string format_double(double num) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << std::setfill('0') << std::setw(9) << std::internal << std::showpos << num;
    return oss.str();
}

void file_output_head(std::ofstream &fo, int ECorB = 0, std::string ModTypS = "", std::string ModRef = "", std::string institute = "", std::string operatorp = "", std::string instrutyp = "", int runt = 0, std::string MeasProgVers = "") {
    fo << "#---Header\n";
    fo << "EC or Barrel:                " << (ECorB ? "SB" : "SE") << "\n";
    fo << "Module Type:                 " << ModTypS << "\n";
    fo << "Module ref. Number:          " << ModRef << "\n";
    fo << "Date:                        " << ISOtime() << "\n";
    fo << "Institute:                   " << institute << "\n";
    fo << "Operator:                    " << operatorp << "\n";
    fo << "Instrument type:             " << instrutyp << "\n";
    fo << "Run Number:                  " << runt << "\n";
    fo << "Measurement program version: " << MeasProgVers << "\n";
}

void file_output_body(std::ofstream &fo, int moduletyp, std::vector<class Point> HyXY, class Point PBP1, class Point PBP2,std::vector<class Point> sensor, std::vector<std::vector<std::vector <class Point>>> gluehcc, std::vector<std::vector<std::vector <class Point>>> glueabc, std::vector<class Point> gluepb) {
    std::string hy1, hy2, hy3, hy4;
    int numofabc0, numofabc1, numofabc2, numofabc3, numofpb = 5;

    switch (moduletyp) {
        case 0:
            hy1 = "R0H0";
            hy2 = "R0H1";
            numofabc0 = 8;
            numofabc1 = 9;
            break;
        case 1:
            hy1 = "R1H0";
            hy2 = "R1H1";
            numofabc0 = 10;
            numofabc1 = 11;
            break;
        case 2:
            hy1 = "R2H0";
            numofabc0 = 12;
            break;
        case 3:
            hy1 = "R3H0";
            hy2 = "R3H1";
            hy3 = "R3H2";
            hy4 = "R3H3";
            numofabc0 = 7;
            numofabc1 = 7;
            numofabc2 = 7;
            numofabc3 = 7;
            break;
        case 4:
            hy1 = "R4H0";
            hy2 = "R4H1";
            numofabc0 = 8;
            numofabc1 = 8;
            break;
        case 5:
            hy1 = "R5H0";
            hy2 = "R5H1";
            numofabc0 = 9;
            numofabc1 = 9;
            break;
        default:
            std::cout << "Wrong module type!\n";
    }

    fo << "#---Positions\n#Location         X [mm]        Y [mm]\n";
    fo << "H_"<<hy1<<"_P1         " << format_double(HyXY[0].x) << "     " << format_double(HyXY[0].y) << "\n";
    fo << "H_"<<hy1<<"_P2         " << format_double(HyXY[1].x) << "     " << format_double(HyXY[1].y) << "\n";
    if (moduletyp != 2) {
        fo << "H_"<<hy2<<"_P1         " << format_double(HyXY[2].x) << "     " << format_double(HyXY[2].y) << "\n";
        fo << "H_"<<hy2<<"_P2         " << format_double(HyXY[3].x) << "     " << format_double(HyXY[3].y) << "\n";
    }
    if (moduletyp == 3) {
        fo << "H_"<<hy3<<"_P1         " << format_double(HyXY[4].x) << "     " << format_double(HyXY[4].y) << "\n";
        fo << "H_"<<hy3<<"_P2         " << format_double(HyXY[5].x) << "     " << format_double(HyXY[5].y) << "\n";
        fo << "H_"<<hy4<<"_P1         " << format_double(HyXY[6].x) << "     " << format_double(HyXY[6].y) << "\n";
        fo << "H_"<<hy4<<"_P2         " << format_double(HyXY[7].x) << "     " << format_double(HyXY[7].y) << "\n";
    }
    fo << "PB_P1             " << format_double(PBP1.x) << "     " << format_double(PBP1.y) << "\n";
    fo << "PB_P2             " << format_double(PBP2.x) << "     " << format_double(PBP2.y) << "\n";

    fo << "#---Glue heights:\n# Location    Type      X [mm]       Y [mm]       Z [mm]\n";
    for (int i=0;i<sensor.size();i++) {fo << "Sensor        1         " << format_double(sensor[i].x) << "    " << format_double(sensor[i].y) << "    " << format_double(sensor[i].z) << "\n";}
    for(int i=0;i<gluehcc[0][0].size();i++) {fo << "H_" << hy1 << "_0      2         " << format_double(gluehcc[0][0][i].x) << "    " << format_double(gluehcc[0][0][i].y) << "    "  << format_double(gluehcc[0][0][i].z) << "\n";}
    for(int i=0;i<gluehcc[0][1].size();i++) {fo << "H_" << hy1 << "_1      2         " << format_double(gluehcc[0][1][i].x) << "    " << format_double(gluehcc[0][1][i].y) << "    "  << format_double(gluehcc[0][1][i].z) << "\n";}
    for(int i=0;i<gluehcc[0][2].size();i++) {fo << "H_" << hy1 << "_2      2         " << format_double(gluehcc[0][2][i].x) << "    " << format_double(gluehcc[0][2][i].y) << "    "  << format_double(gluehcc[0][2][i].z) << "\n";}
    switch (moduletyp) {
        case 0:
            for(int i=0;i<gluehcc[0][3].size();i++) {fo << "HCC" << hy1 << "_2     2         " << format_double(gluehcc[0][3][i].x) << "    " << format_double(gluehcc[0][3][i].y) << "    "  << format_double(gluehcc[0][3][i].z) << "\n";}
            break;
        case 1:
            for(int i=0;i<gluehcc[0][3].size();i++) {fo << "HCC" << hy1 << "_4     2         " << format_double(gluehcc[0][3][i].x) << "    " << format_double(gluehcc[0][3][i].y) << "    "  << format_double(gluehcc[0][3][i].z) << "\n";}
            break;
        case 2:
            for(int i=0;i<gluehcc[0][3].size();i++) {fo << "HCC" << hy1 << "_2     2         " << format_double(gluehcc[0][3][i].x) << "    " << format_double(gluehcc[0][3][i].y) << "    "  << format_double(gluehcc[0][3][i].z) << "\n";}
            for(int i=0;i<gluehcc[0][4].size();i++) {fo << "HCC" << hy1 << "_3     2         " << format_double(gluehcc[0][4][i].x) << "    " << format_double(gluehcc[0][4][i].y) << "    "  << format_double(gluehcc[0][4][i].z) << "\n";}
            break;
        default:
            std::cout << "Wrong module type!\n";
    }
    for(int i=0;i<numofabc0;i++) {for(int j=0;j<glueabc[0][i].size();j++) {fo << "ABC_" << hy1 << "_" << std::left << std::setw(2) << i << "   2         " << format_double(glueabc[0][i][j].x) << "    " << format_double(glueabc[0][i][j].y) << "    "  << format_double(glueabc[0][i][j].z) << "\n";}}
    if (moduletyp != 2) {
        for(int i=0;i<gluehcc[1][0].size();i++) {fo << "H_" << hy2 << "_0      2         " << format_double(gluehcc[1][0][i].x) << "    " << format_double(gluehcc[1][0][i].y) << "    "  << format_double(gluehcc[1][0][i].z) << "\n";}
        for(int i=0;i<gluehcc[1][1].size();i++) {fo << "H_" << hy2 << "_1      2         " << format_double(gluehcc[1][1][i].x) << "    " << format_double(gluehcc[1][1][i].y) << "    "  << format_double(gluehcc[1][1][i].z) << "\n";}
        for(int i=0;i<gluehcc[1][2].size();i++) {fo << "H_" << hy2 << "_2      2         " << format_double(gluehcc[1][2][i].x) << "    " << format_double(gluehcc[1][2][i].y) << "    "  << format_double(gluehcc[1][2][i].z) << "\n";}
        switch (moduletyp) {
            case 0:
                for(int i=0;i<gluehcc[1][3].size();i++) {fo << "HCC" << hy2 << "_3     2         " << format_double(gluehcc[1][3][i].x) << "    " << format_double(gluehcc[1][3][i].y) << "    "  << format_double(gluehcc[1][3][i].z) << "\n";}
                break;
            case 1:
                for(int i=0;i<gluehcc[1][3].size();i++) {fo << "HCC" << hy2 << "_5     2         " << format_double(gluehcc[1][3][i].x) << "    " << format_double(gluehcc[1][3][i].y) << "    "  << format_double(gluehcc[1][3][i].z) << "\n";}
                break;
            case 3:
                for(int i=0;i<gluehcc[1][3].size();i++) {fo << "HCC" << hy2 << "_4     2         " << format_double(gluehcc[1][3][i].x) << "    " << format_double(gluehcc[1][3][i].y) << "    "  << format_double(gluehcc[1][3][i].z) << "\n";}
                for(int i=0;i<gluehcc[1][4].size();i++) {fo << "HCC" << hy2 << "_5     2         " << format_double(gluehcc[1][4][i].x) << "    " << format_double(gluehcc[1][4][i].y) << "    "  << format_double(gluehcc[1][4][i].z) << "\n";}
                break;
            case 4:
                for(int i=0;i<gluehcc[1][3].size();i++) {fo << "HCC" << hy2 << "_2     2         " << format_double(gluehcc[1][3][i].x) << "    " << format_double(gluehcc[1][3][i].y) << "    "  << format_double(gluehcc[1][3][i].z) << "\n";}
                for(int i=0;i<gluehcc[1][4].size();i++) {fo << "HCC" << hy2 << "_3     2         " << format_double(gluehcc[1][4][i].x) << "    " << format_double(gluehcc[1][4][i].y) << "    "  << format_double(gluehcc[1][4][i].z) << "\n";}
                break;
            case 5:
                for(int i=0;i<gluehcc[1][3].size();i++) {fo << "HCC" << hy2 << "_4     2         " << format_double(gluehcc[1][3][i].x) << "    " << format_double(gluehcc[1][3][i].y) << "    "  << format_double(gluehcc[1][3][i].z) << "\n";}
                for(int i=0;i<gluehcc[1][4].size();i++) {fo << "HCC" << hy2 << "_5     2         " << format_double(gluehcc[1][4][i].x) << "    " << format_double(gluehcc[1][4][i].y) << "    "  << format_double(gluehcc[1][4][i].z) << "\n";}
                break;
            default:
                std::cout << "Wrong module type!\n";
        }
        for(int i=0;i<numofabc1;i++) {for(int j=0;j<glueabc[1][i].size();j++) {fo << "ABC_" << hy2 << "_" << std::left << std::setw(2) << i << "   2         " << format_double(glueabc[1][i][j].x) << "    " << format_double(glueabc[1][i][j].y) << "    "  << format_double(glueabc[1][i][j].z) << "\n";}}
    }
    if (moduletyp == 3) {
        for(int i=0;i<gluehcc[2][0].size();i++) {fo << "H_" << hy1 << "_0      2         " << format_double(gluehcc[2][0][i].x) << "    " << format_double(gluehcc[2][0][i].y) << "    "  << format_double(gluehcc[2][0][i].z) << "\n";}
        for(int i=0;i<gluehcc[2][1].size();i++) {fo << "H_" << hy1 << "_1      2         " << format_double(gluehcc[2][1][i].x) << "    " << format_double(gluehcc[2][1][i].y) << "    "  << format_double(gluehcc[2][1][i].z) << "\n";}
        for(int i=0;i<gluehcc[2][2].size();i++) {fo << "H_" << hy1 << "_2      2         " << format_double(gluehcc[2][2][i].x) << "    " << format_double(gluehcc[2][2][i].y) << "    "  << format_double(gluehcc[2][2][i].z) << "\n";}
        for(int i=0;i<gluehcc[2][3].size();i++) {fo << "HCC" << hy3 << "_6     2         " << format_double(gluehcc[2][3][i].x) << "    " << format_double(gluehcc[2][3][i].y) << "    "  << format_double(gluehcc[2][3][i].z) << "\n";}
        for(int i=0;i<gluehcc[2][4].size();i++) {fo << "HCC" << hy3 << "_7     2         " << format_double(gluehcc[2][4][i].x) << "    " << format_double(gluehcc[2][4][i].y) << "    "  << format_double(gluehcc[2][4][i].z) << "\n";}
        for(int i=0;i<numofabc2;i++) {for(int j=0;j<glueabc[2][i].size();j++) {fo << "ABC_" << hy3 << "_" << std::left << std::setw(2) << i << "   2         " << format_double(glueabc[2][i][j].x) << "    " << format_double(glueabc[2][i][j].y) << "    "  << format_double(glueabc[2][i][j].z) << "\n";}}
        for(int i=0;i<gluehcc[3][0].size();i++) {fo << "H_" << hy1 << "_0      2         " << format_double(gluehcc[3][0][i].x) << "    " << format_double(gluehcc[3][0][i].y) << "    "  << format_double(gluehcc[3][0][i].z) << "\n";}
        for(int i=0;i<gluehcc[3][1].size();i++) {fo << "H_" << hy1 << "_1      2         " << format_double(gluehcc[3][1][i].x) << "    " << format_double(gluehcc[3][1][i].y) << "    "  << format_double(gluehcc[3][1][i].z) << "\n";}
        for(int i=0;i<gluehcc[3][2].size();i++) {fo << "H_" << hy1 << "_2      2         " << format_double(gluehcc[3][2][i].x) << "    " << format_double(gluehcc[3][2][i].y) << "    "  << format_double(gluehcc[3][2][i].z) << "\n";}
        for(int i=0;i<numofabc3;i++) {for(int j=0;j<glueabc[3][i].size();j++) {fo << "ABC_" << hy4 << "_" << std::left << std::setw(2) << i << "   2         " << format_double(glueabc[3][i][j].x) << "    " << format_double(glueabc[3][i][j].y) << "    "  << format_double(glueabc[3][i][j].z) << "\n";}}
    }

    for(int i=0;i<numofpb;i++) {fo << "PB_" << i << "          2         " << format_double(gluepb[i].x) << "    " << format_double(gluepb[i].y) << "    "  << format_double(gluepb[i].z) << "\n";}

//======== Other Heights here ========

    fo << "#---Other heights:\n# Location     Type     X [mm]       Y [mm]       Z [mm]\n";
}
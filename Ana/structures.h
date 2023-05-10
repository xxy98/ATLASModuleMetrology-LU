//
// Created by Xiangyu on 2022/11/07.
//

#ifndef ANA_STRUCTURES_H
#define ANA_STRUCTURES_H

#endif //ANA_STRUCTURES_H

enum Feature
{
    EMPTY,      //0
    Plane,      //1
    DatumPlane, //2
    Point,      //3
    Line,       //4
    Intersect,  //5
    DatumOrigin,//6
    Circle      //7
};

class Measurement{
public:
    int step;
    //std::string feat;
    Feature feat;
    double xr;
    double ya;
    double z;
    double size;
    int ref[2];
    //void set(int st, std::string ft, double xr, double ya, double z, double size, int rf1, int rf2);
    void Set(int st, Feature ft, double xr, double ya, double z, double size, int rf1, int rf2);
    void Output();
};

void Measurement::Set(int st, Feature ft, double x, double y, double zz, double sz, int rf1, int rf2) {
    step = st;
    feat = ft;
    xr = x;
    ya = y;
    z = zz;
    size = sz;
    ref[0] = rf1;
    ref[1] = rf2;
}

void Measurement::Output() {
    std::string ft_s;
    switch (feat) {
        case EMPTY:
            ft_s = "EMPTY";
            break;
        case Plane:
            ft_s = "Plane";
            break;
        case DatumPlane:
            ft_s = "Datum Plane";
            break;
        case Point:
            ft_s = "Point";
            break;
        case Line:
            ft_s = "Line";
            break;
        case Intersect:
            ft_s = "Intersect";
            break;
        case DatumOrigin:
            ft_s = "Datum Origin";
            break;
        case Circle:
            ft_s = "Circle";
            break;
    }
    std::cout << "Step: " << step
              << " | " << ft_s << " |";
    if (feat != DatumPlane){
        std::cout << " X: " << xr
                  << " mm Y: " << ya << " mm";
    }
    std::cout << " Z: " << z << " mm";
    if (feat != DatumPlane && feat != Point && feat != DatumOrigin){
        std::cout << " | Size: " << size << " mm";
    }
    if (feat == DatumPlane || feat == DatumOrigin){
        std::cout << " | Ref: " << ref[0];
    }
    if (feat == Line || feat == Intersect){
        std::cout << " | Ref: " << ref[0] << " " << ref[1];
    }
    std::cout << std::endl;
}


//
// Created by Xiangyu on 2022/11/07.
//

#ifndef ANA_MATHEMATICS_H
#define ANA_MATHEMATICS_H

#endif //ANA_MATHEMATICS_H

#ifndef cmath
#include <cmath>
#define cmath
#endif //cmath
#ifndef eigen
#include <Eigen/Dense>
#define eigen
#endif //eigen

//Declaration

double PointDist(class Point, class Point);
class PointsWDistr PointsAve(std::vector<class Point> ps);
class Vector SetVector(class Point, class Point);
class Vector CrossProduct(class Vector, class Vector);
double DotProduct(class Vector, class Vector);
class Point MovePoint(class Point, class Vector);
class Vector MoveVec(class Point, class Vector);
class Vector RevVec(class Vector);
double VecsAngle(class Vector, class Vector, bool);
class Point MidVec(class Vector);
class Line LinkPoints(class Point, class Point);
double PointLineDist(class Point, class Line);
class Point TwoLineIntersection(class Line, class Line);
class Point LinearIntersection(class Line, class Line);
bool IsTwoLineOnPlane(class Line, class Line);
double PointPlaneDist(class Point, class Plane);
double MultiPointsPlaneSD(std::vector<class Point>, class::Plane);
double LineDist (class Line, class Line);
class Vector LineDistVec (class Line, class Line, bool);

//Points begin here.

class Point {
public:
    double x;
    double y;
    double z;
    void SetPoint(double xi, double yi, double zi);
    void Out(int p);
    void Out(std::ostream &f, int p);
    void SetOrigin(class Point op);
};

void Point::SetPoint(double xi, double yi, double zi) {
    x = xi;
    y = yi;
    z = zi;
}

void Point::Out(int p=5) {
    std::cout << std::setprecision(p) << "(" << x << ", " << y << ", " << z << ")";
}

void Point::Out(std::ostream &f, int p=5) {
    f << std::setprecision(5) << "(" << x << ", " << y << ", " << z << ")";
}

void Point::SetOrigin(class Point op) {
    x = x - op.x;
    y = y - op.y;
    z = z - op.z;
}

double PointsDist(class Point p1, class Point p2) {
    return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
}

class PointsWDistr : public Point {
public:
    double xSD;
    double ySD;
    double zSD;
    int count;
    void Out(bool numds = true, bool propds = false);
    void Out(std::ofstream &f, bool numds = true, bool propds = false);
    void Refresh();
};

void PointsWDistr::Out(bool numds, bool propds) {
    std::cout << "There're " << count << " sample(s).\n";
    std::cout << std::fixed << std::setprecision(5) << "AVE: (" << std::showpos << std::setw(11) << x << ", "
                                                                   << std::showpos << std::setw(11) << y << ", "
                                                                   << std::showpos << std::setw(11) << z << ")\n";
    if(numds){std::cout << std::fixed << std::setprecision(5) << "SD:  (" << std::noshowpos << std::setw(11) << xSD << ", "
                                                                             << std::noshowpos << std::setw(11) << ySD << ", "
                                                                             << std::noshowpos << std::setw(11) << zSD << ")\n";}
    if(propds){std::cout << std::scientific << std::setprecision(5) << "RSD: (" << std::noshowpos << std::setw(11) << fabs(xSD/x) << ", "
                                                                                   << std::noshowpos << std::setw(11) << fabs(ySD/y) << ", "
                                                                                   << std::noshowpos << std::setw(11) << fabs(zSD/z) << ")\n";}
}

void PointsWDistr::Out(std::ofstream &f, bool numds, bool propds) {
    f << "There're " << count << " sample(s).\n";
    f << std::fixed << std::setprecision(5) << "AVE: (" << std::showpos << std::setw(11) << x << ", "
                                                           << std::showpos << std::setw(11) << y << ", "
                                                            << std::showpos << std::setw(11) << z << ")\n";
    if(numds){f << std::fixed << std::setprecision(5) << "SD:  (" << std::noshowpos << std::setw(11) << xSD << ", "
                                                                     << std::noshowpos << std::setw(11) << ySD << ", "
                                                                     << std::noshowpos << std::setw(11) << zSD << ")\n";}
    if(propds){f << std::scientific << std::setprecision(5) << "RSD: (" << std::noshowpos << std::setw(11) << fabs(xSD/x) << ", "
                                                                           << std::noshowpos << std::setw(11) << fabs(ySD/y) << ", "
                                                                           << std::noshowpos << std::setw(11) << fabs(zSD/z) << ")\n";}
}

class PointsWDistr PointsAve(std::vector<class Point>ps) {
    class PointsWDistr pr;
    int num = ps.size();
    double xa = 0, ya = 0, za = 0, xd = 0, yd = 0, zd = 0;
    pr.count = num;
    for(int i = 0; i<num; i++) {
        xa += ps[i].x / num;
        ya += ps[i].y / num;
        za += ps[i].z / num;
    }
    for(int i = 0; i<num; i++) {
        xd += pow(ps[i].x - xa,2) / num;
        yd += pow(ps[i].y - ya,2) / num;
        zd += pow(ps[i].z - za,2) / num;
    }
    xd = sqrt(xd);
    yd = sqrt(yd);
    zd = sqrt(zd);
    pr.SetPoint(xa,ya,za);
    pr.xSD = xd;
    pr.ySD = yd;
    pr.zSD = zd;
    return pr;
}

class Circle : public Point {
public:
    double r;
    void SetCircle(double xi, double yi, double zi, double ri);
    void Out();
    void Out(std::ostream &f);
};

void Circle::SetCircle(double xi, double yi, double zi, double ri) {
    x = xi;
    y = yi;
    z = zi;
    r = ri;
}

void Circle::Out() {
    std::cout << "(" << x << ", " << y << ", " << z << "; " << r << ")";
}

void Circle::Out(std::ostream &f) {
    f << "(" << x << ", " << y << ", " << z << "; " << r << ")";
}

//Vectors begin here.

class Vector {
public:
    class Point op;
    class Point dp;
    double vx;
    double vy;
    double vz;
    double Length();
    void SetDp();
    void SetV();
    void Mag(double m);
    bool NormCheck();
    void Norm();
    void Out(int p, bool opo, bool diro, bool dpo, bool leno);
    void Out(std::ostream &f, int p, bool opo, bool diro, bool dpo, bool leno);
};

double Vector::Length() {
    return sqrt(vx*vx + vy*vy + vz*vz);
}

void Vector::SetDp() {
    dp.x = op.x + vx;
    dp.y = op.y + vy;
    dp.z = op.z + vz;
}

void Vector::SetV() {
    vx = dp.x - op.x;
    vy = dp.y - op.y;
    vz = dp.z - op.z;
}

void Vector::Mag(double m){
    vx = vx * m;
    vy = vy * m;
    vz = vz * m;
    SetDp();
}

bool Vector::NormCheck() {
    return ((Length()==1)?true:false);
}

void Vector::Norm() {
    Mag(1 / Length());
}

void Vector::Out(int p=5, bool opo = true, bool diro = true, bool dpo = false, bool leno = true) {
    if (opo) {
        std::cout << "| Ini. p.: ";
        op.Out(p);
        std::cout << " ";
    }
    if (diro) {
        std::cout << "| Dir: (" << vx << ", " << vy << ", " << vz <<") ";
    }
    if (dpo) {
        std::cout << "| Ter. p.: ";
        dp.Out(p);
        std::cout << " ";
    }
    if (leno) {
        std::cout << "| Len.:" << Length() << " ";
    }
    std::cout << "|\n";
}

void Vector::Out(std::ostream &f, int p=5, bool opo = true, bool diro = true, bool dpo = false, bool leno = true) {
    if (opo) {
        f << "| Ini. p.: ";
        op.Out(f, p);
        f << " ";
    }
    if (diro) {
        f << "| Dir: (" << vx << ", " << vy << ", " << vz <<") ";
    }
    if (dpo) {
        f << "| Ter. p.: ";
        dp.Out(f, p);
        f << " ";
    }
    if (leno) {
        f << "| Len.:" << Length() << " ";
    }
    f << "|\n";
}

class Vector SetVector(class Point p1, class Point p2) {
    class Vector v;
    v.op = p1;
    v.dp = p2;
    v.vx = p2.x - p1.x;
    v.vy = p2.y - p1.y;
    v.vz = p2.z - p1.z;
    return v;
}

class Vector CrossProduct(class Vector v1, class Vector v2) {
    class Vector vec;
    class Point opi;
    opi.SetPoint(0,0,0);
    vec.op= opi;
    vec.vx = v1.vy*v2.vz - v1.vz*v2.vy;
    vec.vy = v1.vz*v2.vx - v1.vx*v2.vz;
    vec.vz = v1.vx*v2.vy - v1.vy*v2.vx;
    vec.SetDp();
    return vec;
}

double DotProduct(class Vector v1, class Vector v2) {
    return (v1.vx*v2.vx + v1.vy*v2.vy + v1.vz*v2.vz);
}

class Point MovePoint(class Point op, class Vector vec) {
    class Point dp;
    dp.x = op.x + vec.vx;
    dp.y = op.y + vec.vy;
    dp.z = op.z + vec.vz;
    return dp;
}

class Vector MoveVec(class Point op, class Vector vec) {
    class Vector v;
    v.op = op;
    v.vx = vec.vx;
    v.vy = vec.vy;
    v.vz = vec.vz;
    v.SetDp();
    return v;
}

class Vector RevVec(class Vector vec) {
    class Vector rev;
    rev.op = vec.dp;
    rev.dp = vec.op;
    rev.SetV();
    return rev;
}

double VecsAngle(class Vector v1, class Vector v2, bool Rad=true) {
    if(Rad){
        return acos((DotProduct(v1, v2))/((v1.Length())*(v2.Length())));
    } else {
        return acos((DotProduct(v1, v2))/((v1.Length())*(v2.Length())))/PI*180;
    }

}

class Point MidVec(class Vector vec){
    class Point p;
    p.SetPoint((vec.op.x+vec.dp.x)/2,(vec.op.y+vec.dp.y)/2,(vec.op.z+vec.dp.z)/2);
    return p;
}

//Lines begin here.

class Line {
public:
    class Point op1;
    class Point op2;
    class Vector vec;
    void SetWithTwoPoints();
    void SetWithVec();
};

void Line::SetWithTwoPoints() {
    vec = SetVector(op1, op2);
}

void Line::SetWithVec() {
    op2.x = op1.x + vec.vx;
    op2.y = op1.y + vec.vy;
    op2.z = op1.z + vec.vz;
}

class Line LinkPoints(class Point p1, class Point p2){
    class Line l;
    l.vec = SetVector(p1, p2);
    l.op1 = p1;
    l.op2 = p2;
    return l;
}

double PointLineDist(class Point p, class Line l) {
    class Vector vht = SetVector(l.op1, p);
    class Vector vtp = CrossProduct(vht, l.vec);
    return vtp.Length()/l.vec.Length();
}

class Point TwoLineIntersection(class Line l1, class Line l2) {
    class Point Is;
    class Vector v1 = l1.vec;
    class Vector v2 = l2.vec;
    v1.Norm();
    v2.Norm();
    double per = PointLineDist(l2.op1, l1);
    double htl = PointsDist(l2.op1, l1.op1);
    double sdl = sqrt(htl*htl - per*per);
    double ipl = per / tan(VecsAngle(l1.vec, l2.vec));
    class Point ppp, ppn, pnp, pnn, pp, pn;
    v1.Mag(sdl+ipl);
    ppp = MovePoint(l1.op1, v1);
    double pppd = PointLineDist(ppp,l2);
    v1.Norm();
    v1.Mag(sdl-ipl);
    ppn = MovePoint(l1.op1, v1);
    double ppnd = PointLineDist(ppn,l2);
    v1.Norm();
    v1.Mag(-1*sdl+ipl);
    pnp = MovePoint(l1.op1, v1);
    double pnpd = PointLineDist(pnp,l2);
    v1.Norm();
    v1.Mag(-1*sdl-ipl);
    pnn = MovePoint(l1.op1, v1);
    double pnnd = PointLineDist(pnn,l2);
    v1.Norm();
    double ppd, pnd, d;
    /*
    ppp.Out();std::cout<<pppd<<std::endl;
    ppn.Out();std::cout<<ppnd<<std::endl;
    pnp.Out();std::cout<<pnpd<<std::endl;
    pnn.Out();std::cout<<pnnd<<std::endl;
    */
    if (pppd<ppnd){
        pp = ppp;
        ppd = pppd;
    } else {
        pp = ppn;
        ppd = ppnd;
    }
    if (pnpd<pnnd){
        pn = pnp;
        pnd = pnpd;
    } else {
        pn = pnn;
        pnd = pnnd;
    }
    if (ppd<pnd){
        Is = pp;
        d = ppd;
    } else {
        Is = pn;
        d = pnd;
    }
    return Is;
}

class Point LinearIntersection(class Line line1, class Line line2) {
    double x1, x2, l, l1, l2;
    double d11, d12, d21, d22, d1, d2;
    class Point p11, p12, p21, p22, p1, p2;
    class Vector v = line1.vec;
    v.Norm();
    x1 = PointLineDist(line1.op1, line2);
    x2 = PointLineDist(line1.op2, line2);
    l = line1.vec.Length();
    l1 = l * (x1 / (x1 - x2));
    l2 = l * (x1 / (x1 + x2));
    v.Mag(l1);
    p11 = MovePoint(line1.op1, v);
    d11 = PointLineDist(p11, line2);
    v.Norm();
    v.Mag(l1*-1);
    p12 = MovePoint(line1.op1,v);
    d12 = PointLineDist(p12, line2);
    v.Norm();
    v.Mag(l2);
    p21 = MovePoint(line1.op1,v);
    d21 = PointLineDist(p21, line2);
    v.Norm();
    v.Mag(l2*-1);
    p22 = MovePoint(line1.op1,v);
    d22 = PointLineDist(p22, line2);
    p1 = d11<d12?p11:p12;
    d1 = d11<d12?d11:d12;
    p2 = d21<d22?p21:p22;
    d2 = d21<d22?d21:d22;
    return d1<d2?p1:p2;
}

bool IsTwoLineOnPlane(class Line l1, class Line l2) {
    return (fabs(((l2.op1.x - l1.op1.x) * l1.vec.vy * l2.vec.vz)
                 + ((l2.op1.y - l1.op1.y) * l1.vec.vz * l2.vec.vx)
                 + ((l2.op1.z - l1.op1.z) * l1.vec.vx * l2.vec.vy)
                 - ((l2.op1.x - l1.op1.x) * l1.vec.vz * l2.vec.vy)
                 - ((l2.op1.y - l1.op1.y) * l1.vec.vx * l2.vec.vz)
                 - ((l2.op1.z - l1.op1.z) * l1.vec.vy * l2.vec.vx)) == 0);
}

class Plane {
public:
    class Point op;
    class Vector normvec;
    void SetPlane(class Point p, class Vector v);
    void ThreePointSetPlane(class Point p1, class Point p2, class Point p3);
    void MultiPointSetPlaneSVD(std::vector<class Point>ps);
    bool IsOnPlane(class Point tp);
};

void Plane::SetPlane(class Point p, class Vector v) {
    op = p;
    normvec = v;
}

void Plane::ThreePointSetPlane(class Point p1, class Point p2, class Point p3) {
    class Vector pv1, pv2, pvr;
    pv1 = SetVector(p1, p2);
    pv2 = SetVector(p1, p3);
    pvr = CrossProduct(pv1, pv2);
    pvr.Norm();
    op = p1;
    normvec = pvr;
}

void Plane::MultiPointSetPlaneSVD(std::vector<class Point> ps) {
    normvec.op.SetPoint(0,0,0);
    const int num = ps.size();
    Eigen::MatrixXd PointCloud(num, 3);
    for (int i=0;i<num;i++) {
        PointCloud(i,0) = ps[i].x;
        PointCloud(i,1) = ps[i].y;
        PointCloud(i,2) = ps[i].z;
    }
#ifdef debug
    std::cout << PointCloud << std::endl;
#endif
    Eigen::RowVector3d centroid = PointCloud.colwise().mean();
    Eigen::MatrixXd demean = PointCloud;
    demean.rowwise() -= centroid;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(demean, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd V = svd.matrixV();
    //Eigen::MatrixXd U = svd.matrixU();
    //Eigen::MatrixXd S = U.inverse() * demean * V.transpose().inverse();

    Eigen::RowVector3d normal;
    normal << V(0,2), V(1,2), V(2,2);

    double d = -normal * centroid.transpose();
    normvec.dp.SetPoint(V(0,2),V(1,2),V(2,2));
    normvec.SetV();
    op.SetPoint(0,0,-1*d/V(2,2));
    //op.SetPoint(-1*d/3/V(0,2),-1*d/3/V(1,2),-1*d/3/V(2,2));
}

bool Plane::IsOnPlane(class Point tp) {
    class Vector vec = SetVector(op, tp);
    return ((fabs(DotProduct(vec, normvec))==0)?true:false);
}

double PointPlaneDist(class Point p, class Plane pl) {
    class Vector v;
    v = SetVector(pl.op, p);
    return fabs(DotProduct(v,pl.normvec)/pl.normvec.Length());
}

double MultiPointsPlaneSD(std::vector<class Point> ps, class::Plane p) {
    const int num = ps.size();
    double diffsum;
    for (int i=0;i<num;i++) {diffsum += pow(PointPlaneDist(ps[i], p),2);}
    return std::sqrt(diffsum);
}

double LineDist (class Line l1, class Line l2){
    class Vector DistVec = CrossProduct(l1.vec, l2.vec);
    //DistVec.Mag(1/l1.vec.Length()/l2.vec.Length());
    DistVec.Norm();
    class Vector ArbiVec = SetVector(l1.op1, l2.op1);
    return fabs(DotProduct(DistVec, ArbiVec));
}

class Vector LineDistVec (class Line l1, class Line l2, bool fl = true){
    class Vector Dist;
    class Vector DistVec = CrossProduct(l1.vec, l2.vec);
    if (LineDist(l1,l2)==0) {
        DistVec.Mag(0);
        if (fl) {
            DistVec.op = TwoLineIntersection(l1,l2);
        } else {
            DistVec.op = LinearIntersection(l1,l2);
        }

        DistVec.SetDp();
        return DistVec;
    }
    DistVec.Norm();
    DistVec.Mag(LineDist(l1,l2));
    //DistVec.Mag(LineDisd(l1,l2)/DistVec.Length());
    class Line l1Par, l1ParP, l1ParN;
    class Point IP1, IP2;
    double distP, distN;
    int mag = 0;
    l1ParP.op1 = MovePoint(l1.op1, DistVec);
    l1ParP.op2 = MovePoint(l1.op2, DistVec);
    l1ParP.SetWithTwoPoints();
    distP = LineDist(l1ParP,l2);
    DistVec.Mag(-1);
    l1ParN.op1 = MovePoint(l1.op1, DistVec);
    l1ParN.op2 = MovePoint(l1.op2, DistVec);
    l1ParN.SetWithTwoPoints();
    distN = LineDist(l1ParN,l2);
    DistVec.Mag(-1);
    if (distP < distN) {
        l1Par = l1ParP;
        mag = 1;
    } else {
        l1Par = l1ParN;
        mag = -1;
    }
    //std::cout<<LineDist(l1ParP,l2)<<" "<<LineDist(l1,l2)<<" "<<LineDist(l1ParN,l2)<<std::endl;
    //std::cout<<LineDist(l1Par,l2)<<","<<LineDist(l1,l2)<<","<<LineDist(l1Par,l2)/ LineDist(l1,l2)<<std::endl;
    DistVec.Mag(mag);
    if(fl) {
        IP2 = LinearIntersection(l1Par, l2);
    } else {
        IP2 = TwoLineIntersection(l1Par, l2);
    }
    IP1 = MovePoint(IP2, RevVec(DistVec));
    Dist.op = IP1;
    Dist.dp = IP2;
    Dist.SetV();
    return Dist;
}
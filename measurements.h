//
// Created by Xiangyu on 2023/04/11.
//

#ifndef ANA_MEASUREMENTS_H
#define ANA_MEASUREMENTS_H

#endif //ANA_MEASUREMENTS_H

#ifndef ANA_MATHEMATICS_H
#include "mathematics.h"
#define ANA_MATHEMATICS_H
#endif

class Point CrossFiducial(class Point[8], int[]);
class Point SquareFiducial(class Point[4], int[]);

int defaultCrossFiducialArr[8] = {0,1,2,3,4,5,6,7};
class Point CrossFiducial(class Point p[8], int seq[8] = defaultCrossFiducialArr) {
    class Line edge[4], cross[2];
    class Point interedge[4], center;
    edge[0] = LinkPoints(p[seq[0]], p[seq[1]]);
    edge[1] = LinkPoints(p[seq[2]], p[seq[3]]);
    edge[2] = LinkPoints(p[seq[4]], p[seq[5]]);
    edge[3] = LinkPoints(p[seq[6]], p[seq[7]]);
    interedge[0] = TwoLineIntersection(edge[3], edge[0]);
    interedge[1] = TwoLineIntersection(edge[0], edge[1]);
    interedge[2] = TwoLineIntersection(edge[1], edge[2]);
    interedge[3] = TwoLineIntersection(edge[2], edge[3]);
    cross[0] = LinkPoints(interedge[0], interedge[2]);
    cross[1] = LinkPoints(interedge[1], interedge[3]);
    center = TwoLineIntersection(cross[0], cross[1]);
    return center;
}

int defaultSquareFiducialArr[4] = {0,1,2,3,};
class Point SquareFiducial(class Point p[4], int seq[4] = defaultSquareFiducialArr) {
    class Line edge[2];
    class Point center;
    edge[0] = LinkPoints(p[seq[0]], p[seq[1]]);
    edge[1] = LinkPoints(p[seq[2]], p[seq[3]]);
    center = TwoLineIntersection(edge[0], edge[1]);
    return center;
}

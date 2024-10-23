#ifndef PHYSTRUCTS_H
#define PHYSTRUCTS_H

#include "VectorND.h"

template<typename Float,int DIM>
struct Particle {
    VectorND<Float,DIM> pos;
    VectorND<Float,DIM> vel;
    VectorND<Float,DIM> posnew;
    VectorND<Float,DIM> velnew;
    int collision;
    float lasts;
    float thiss;
};

#endif

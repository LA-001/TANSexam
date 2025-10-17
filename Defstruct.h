#ifndef DEFSTRUCT_H
#define DEFSTRUCT_H

struct Hit {
    double r;
    double phi;
    double z;
    int etichetta;
};

struct Vrt {
    double x0, y0, z0;
    int moltiplicita;
};

struct Ric {
    double z0,zRic;
    int molti;
};

#endif
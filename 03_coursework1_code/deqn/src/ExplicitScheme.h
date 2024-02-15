#ifndef EXPLICIT_SCHEME_H_
#define EXPLICIT_SCHEME_H_

#include "Scheme.h"
#include "InputFile.h"

class ExplicitScheme : public Scheme {
    private:
        Mesh* mesh;

        void updateBoundaries();
        void advance(); //replaces reset
        void diffuse(double dt);
        void reflectBoundaries(int cell, int boundary_id);

        //void talkToSurroundingCells(int cell);
        void getNeighbourCellData(int cell, int boundary_id);
    public:
        ExplicitScheme(const InputFile* input, Mesh* m);

        void doAdvance(const double dt);

        void init();
};
#endif

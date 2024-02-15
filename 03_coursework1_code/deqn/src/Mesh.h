#ifndef DIFFUSION_MESH_H_
#define DIFFUSION_MESH_H_

#include "InputFile.h"

class Mesh {
    private:
        const InputFile* input;

        double*** uX; //pointer to array of frames, each frame array of cells
        int numFrames; //need for array length
        int currentFrame;
        int* cellSize;

        double** posX; //in each cell
        double** posY;
        double* posGlobalX; //over entire field
        double* posGlobalY;

        double* min_coords;
        double* max_coords;

        int NDIM;

        int* n; 
        int* min;
        int* max;

        double* dx;

        /*
         * A mesh has four neighbours, and they are 
         * accessed in the following order:
         * - top
         * - right
         * - bottom
         * - left
         */

        void allocate();
        bool allocated;
    public:
        Mesh(const InputFile* input);
                
        double** getU0();
        double** getU1();
        double** getUX(int frame); //my function

        double* getDx();
        int* getNx();
        int* getMin();
        int* getMax();
        double* getMinCoord(); //added
        double* getMaxCoord(); //added
        int getDim();

        double** getPosInCellX();
        double** getPosInCellY();
        double* getPosGlobalX();
        double* getPosGlobalY();

        int* getNeighbours();

        double getTotalTemperature();

        //my added functions
        void advance();
        int getCurrentFrame();
        int* getCellSize();

        //index translation functions, will be called as little as possible as not very fast, (get original index)
        int getOI(int cell, int index);
        int getOIi(int cell, int index);
        int getOIj(int cell, int index);

};
#endif

#include "Mesh.h"

#include <cstdlib>
#include <iostream>
#include <omp.h>

#include "DIVISIONS.h"

#define POLY2(i, j, imin, jmin, ni) (((i) - (imin)) + ((j)-(jmin)) * (ni))


Mesh::Mesh(const InputFile* input):
    input(input)
{
    allocated = false;

    NDIM = 2;

    n = new int[NDIM];
    min = new int[NDIM];
    max = new int[NDIM];
    dx = new double[NDIM];

    int nx = input->getInt("nx", 0);
    int ny = input->getInt("ny", 0);

    //hacky way to fix divisions by forcing nx and ny to be multiple
    nx -= (nx % DIVISIONS);
    ny -= (ny % DIVISIONS);

    min_coords = new double[NDIM];
    max_coords = new double[NDIM];

    min_coords[0] = input->getDouble("xmin", 0.0);
    max_coords[0] = input->getDouble("xmax", 1.0);
    min_coords[1] = input->getDouble("ymin", 0.0);
    max_coords[1] = input->getDouble("ymax", 1.0);

    cellSize = new int[NDIM]; //added

    // setup first dimension.
    n[0] = nx;
    min[0] = 1;
    max[0] = nx;

    dx[0] = ((double) max_coords[0]-min_coords[0])/nx;

    cellSize[0] = (nx/DIVISIONS) + 2; //halo

    // setup second dimension.
    n[1] = ny;
    min[1] = 1;
    max[1] = ny;

    dx[1] = ((double) max_coords[1]-min_coords[1])/ny;
    
    cellSize[1] = (ny/DIVISIONS) + 2; //halo

    float endTime = input->getDouble("end_time", 0.0);
    float stepTime = input->getDouble("initial_dt", 0.0);
    numFrames = int(endTime / stepTime) + 1;
    currentFrame = 0;

    allocate();
}

void Mesh::allocate()
{
    allocated = true;

    int nx = n[0];
    int ny = n[1];

    /* Allocate and initialise coordinate arrays */
    posX = new double*[NUM_CELLS];
    posY = new double*[NUM_CELLS];

    double xmin = min_coords[0];
    double ymin = min_coords[1];

    #pragma omp parallel num_threads(NUM_CELLS)
    {      
        int cell = omp_get_thread_num();

        int cellX = cell % DIVISIONS;
        int cellY = cell / DIVISIONS;

        posX[cell] = new double[cellSize[0]];
        for (int j = 0; j < cellSize[0]; j++){
            posX[cell][j] = xmin + cellX*dx[0]*(cellSize[0]-2) + dx[0]*(j-1);
        }
        posY[cell] = new double[cellSize[1]];
        for (int i = 0; i < cellSize[1]; i++){
            posY[cell][i] = ymin + cellY*dx[1]*(cellSize[1]-2) + dx[1]*(i-1);
        }
    }

    //renamed cellX,cellY to not confuse with physics cells
    posGlobalX = new double[nx+2];
    posGlobalY = new double[ny+2];

    for (int j = 0; j < nx + 2; j++){
        posGlobalX[j] = xmin + dx[0]*(j-1);
    }
    for (int i = 0; i < ny + 2; i++){
        posGlobalY[i] = ymin + dx[1]*(i-1);
    }

    //allocate frames
    int uXsize = std::min(numFrames, omp_get_max_threads());
    uX = new double**[uXsize];
    for (int i = 0; i < uXsize; i++){
        /* Allocate cell pointers */
        uX[i] = new double*[NUM_CELLS];
        for (int j = 0; j < NUM_CELLS; j++){
            uX[i][j] = new double[cellSize[0] * cellSize[1]];
        }
    }
}

double** Mesh::getU0()
{
    if (currentFrame >= numFrames){
        std::cout << "Current frame too high! (getU0, current frame = " << currentFrame << ", num frames = " << numFrames << ")" << std::endl;
        throw;
    }
    return uX[currentFrame % omp_get_max_threads()];
}

double** Mesh::getU1()
{
    if (currentFrame + 1 >= numFrames){
        std::cout << "Current frame too high! (getU1, current frame + 1 = " << currentFrame + 1 << ", num frames = " << numFrames << ")" << std::endl;
        throw;
    }
    return uX[(currentFrame + 1) % omp_get_max_threads()];
}

double** Mesh::getUX(int frame)
{
    if (frame >= numFrames || frame < 0){
        std::cout << "Out of bounds of UX: " << frame << ". Total frames = " << numFrames << std::endl;
        throw;
    }
    return uX[frame % omp_get_max_threads()];
}

int Mesh::getCurrentFrame()
{
    return currentFrame;
}

double* Mesh::getDx()
{
    return dx;
}

int* Mesh::getMin()
{
    return min;
}

int* Mesh::getMax()
{
    return max;
}

int Mesh::getDim()
{
    return NDIM;
}

int* Mesh::getNx()
{
    return n;
}

double* Mesh::getMinCoord(){
    return min_coords;
}
double* Mesh::getMaxCoord(){
    return max_coords;
}

double** Mesh::getPosInCellX()
{
    return posX;
}

double** Mesh::getPosInCellY()
{
    return posY;
}

double* Mesh::getPosGlobalX()
{
    return posGlobalX;
}

double* Mesh::getPosGlobalY()
{
    return posGlobalY;
}


double Mesh::getTotalTemperature()
{
    if(allocated) {

        int cellSizeX = getCellSize()[0];
        int cellSizeY = getCellSize()[1];

        double** u0 = getU0();
        
        double temperature = 0.0;
        #pragma omp parallel num_threads(NUM_CELLS)
        {
            int cellNum = omp_get_thread_num();  
            double* u0cell = u0[cellNum];

            double localTemperature = 0.0;

            for (int i = 1; i < cellSizeY - 1; i++){
                for (int j = 1; j < cellSizeX - 1; j++){
                    localTemperature += u0cell[i * cellSizeX + j];
                }
            }

            #pragma omp atomic
            temperature+=localTemperature;
        }

        return temperature;
    } else {
        return 0.0;
    }
}

//my frame increase function
void Mesh::advance(){
    currentFrame++;
}

int Mesh::getOI(int cell, int index)
{
    int j = index % getCellSize()[0]; //x
    int i = index / getCellSize()[0]; //y
    int cellJ = cell % DIVISIONS;
    int cellI = cell / DIVISIONS;

    j += cellJ * getCellSize()[0];
    i += cellI * getCellSize()[1];

    return i * getNx()[0] + j;
}

int Mesh::getOIi(int cell, int index)
{
    int i = index / getCellSize()[0];
    int cellI = cell / DIVISIONS;
    i += cellI * getCellSize()[1];

    return i;
}

int Mesh::getOIj(int cell, int index)
{
    int j = index % getCellSize()[0]; //x
    int cellJ = cell % DIVISIONS;
    j += cellJ * getCellSize()[0];

    return j;
}

int* Mesh::getCellSize(){
    return cellSize;
}
#include "ExplicitScheme.h"

#include <iostream>
#include <omp.h>

#include "DIVISIONS.h"

#define POLY2(i, j, imin, jmin, ni) (((i) - (imin)) + (((j)-(jmin)) * (ni)))

static double totalSeconds = 0;

ExplicitScheme::ExplicitScheme(const InputFile* input, Mesh* m) :
    mesh(m)
{
    int nx = mesh->getNx()[0];
    int ny = mesh->getNx()[1];
}

void ExplicitScheme::doAdvance(const double dt)
{
    diffuse(dt);

    updateBoundaries();

    advance();
}

void ExplicitScheme::updateBoundaries()
{
    bool top, right, bottom, left;

    //cant multi thread this as accesses other cells
    for (int cell = 0; cell < NUM_CELLS; cell++){
        top = cell < DIVISIONS;
        right = cell % DIVISIONS == DIVISIONS - 1;
        left = cell % DIVISIONS == 0;
        bottom = cell >= DIVISIONS*(DIVISIONS - 1);

        //share neighbour cell data to the halo so that in cell calculations are accurate
        if (!top)
            getNeighbourCellData(cell, 0);

        if (!right)
            getNeighbourCellData(cell, 1);

        if (!bottom)
            getNeighbourCellData(cell, 2);

        if (!left)
            getNeighbourCellData(cell, 3);

        //reflect boundaries
        if (top)
            reflectBoundaries(cell, 0);

        if (right)
            reflectBoundaries(cell, 1);

        if (bottom)
            reflectBoundaries(cell, 2);

        if (left)
            reflectBoundaries(cell, 3);
    }
}

void ExplicitScheme::init()
{
    updateBoundaries();
}

void ExplicitScheme::advance()
{
    mesh->advance();
}

void ExplicitScheme::diffuse(double dt)
{
    double** u0 = mesh->getU0();
    double** u1 = mesh->getU1();

    double dx = mesh->getDx()[0];
    double dy = mesh->getDx()[1];

    double rx = dt/(dx*dx);
    double ry = dt/(dy*dy);

    int cellSizeX = mesh->getCellSize()[0];
    int cellSizeY = mesh->getCellSize()[1];

    #pragma omp parallel num_threads(NUM_CELLS)
    {
        int cellNum = omp_get_thread_num();  
        double* u0cell = u0[cellNum];
    
        for (int n = cellSizeX + 1; n < cellSizeX*(cellSizeY-1); n++){
            if (n % cellSizeX == 0 || n % cellSizeX == cellSizeX - 1)
                continue;

            u1[cellNum][n] = (1.0-2.0*rx-2.0*ry)*u0cell[n] + rx*u0cell[n - 1] + rx*u0cell[n + 1]
                + ry*u0cell[n - cellSizeX] + ry*u0cell[n + cellSizeX];
        }
    }
}


void ExplicitScheme::getNeighbourCellData(int cell, int boundary_id)
{
    double** u1 = mesh->getU1();
    double* u1cell = mesh->getU1()[cell]; 

    int xSize = mesh->getCellSize()[0];    
    int ySize = mesh->getCellSize()[1];

    int neighbourCell;
    int offsetFix;
    switch(boundary_id){
        //top
        case 0:
            neighbourCell = cell - DIVISIONS;
            offsetFix = xSize*(ySize-2);
            for (int j = 0; j < xSize; j++){
                u1cell[j] = u1[neighbourCell][j + offsetFix]; //top row = second bottom row of above cell
            }
            break;
        //right
        case 1:
            neighbourCell = cell + 1;
            offsetFix = -(xSize-2);
            for (int i = xSize - 1; i < xSize*ySize; i+=xSize){
                u1cell[i] = u1[neighbourCell][i + offsetFix]; //right column = second left column of cell to the right
            }
            break;
        //bottom
        case 2:
            neighbourCell = cell + DIVISIONS;
            offsetFix = -(xSize*(ySize - 2));
            for (int j = xSize*(ySize - 1); j < xSize*ySize; j++){
                u1cell[j] = u1[neighbourCell][j + offsetFix]; //bottom row = second top row of cell below
            }
            break;
        //left
        case 3:
            neighbourCell = cell - 1;
            offsetFix = (xSize - 2);
            for (int i = 0; i < xSize*ySize; i+=xSize){
                u1cell[i] = u1[neighbourCell][i + offsetFix]; //left column = second right column of cell to the left
            }
            break;
        default: std::cerr << "Error in reflectBoundaries(): unknown boundary id (" << boundary_id << ")" << std::endl;
    }


}



void ExplicitScheme::reflectBoundaries(int cell, int boundary_id)
{
    double* u1cell = mesh->getU1()[cell]; 

    int xSize = mesh->getCellSize()[0];  
    int ySize = mesh->getCellSize()[1];  

    switch(boundary_id) {
        case 0: 
            /* top */
            {
                for (int j = 0; j < xSize; j++){
                    u1cell[j] = u1cell[j + xSize]; //top boundary = row below
                }
            } break;
        case 1:
            /* right */
            {
                for (int i = xSize - 1; i < xSize*ySize; i+=xSize){
                    u1cell[i] = u1cell[i - 1]; //right boundary = column to the left
                }
            } break;
        case 2: 
            /* bottom */
            {
                for (int j = xSize*(ySize - 1); j < xSize*ySize; j++){
                    u1cell[j] = u1cell[j - xSize]; //bottom boundary = row above
                }
            } break;
        case 3: 
            /* left */
            {
                for (int i = 0; i < xSize*ySize; i+=xSize){
                    u1cell[i] = u1cell[i + 1]; //left boundary = column to the right
                }
            } break;
        default: std::cerr << "Error in reflectBoundaries(): unknown boundary id (" << boundary_id << ")" << std::endl;
    }
}

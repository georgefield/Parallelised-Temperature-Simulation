#include "Diffusion.h"

#include "ExplicitScheme.h"

#include <iostream>
#include <cstdlib>
#include <time.h>
#include <omp.h>

#include "DIVISIONS.h"

Diffusion::Diffusion(const InputFile* input, Mesh* m) :
    mesh(m) 
{

    std::string scheme_str = input->getString("scheme", "explicit");

    if(scheme_str.compare("explicit") == 0) {
        scheme = new ExplicitScheme(input, mesh);
    } else {
        std::cerr << "Error: unknown scheme \"" << scheme_str << "\"" << std::endl;
        exit(1);
    }

    subregion = input->getDoubleList("subregion", std::vector<double>());

    if (subregion.size() != 0 && subregion.size() != 4) {
        std::cerr << "Error:  region must have 4 entries (xmin, ymin, xmax, ymax)" << std::endl;
        exit(1);
    }

    init();
}

Diffusion::~Diffusion()
{
    delete scheme;
}

void Diffusion::init()
{
    double** u0 = mesh->getU0();

    double** posX = mesh->getPosInCellX();
    double** posY = mesh->getPosInCellY();

    int cellSizeX = mesh->getCellSize()[0];
    int cellSizeY = mesh->getCellSize()[1];

    if(!subregion.empty()) {

        #pragma omp parallel num_threads(NUM_CELLS)
        {      
            int cell = omp_get_thread_num();
            bool inCell = false;

            for (int i = 0; i < cellSizeY; i++){
                for (int j = 0; j < cellSizeX; j++){
                    int cellIndex = i*cellSizeX + j;

                    if (posX[cell][j] > subregion[0] 
                    && posX[cell][j] <= subregion[2] 
                    && posY[cell][i] > subregion[1] 
                    && posY[cell][i] <= subregion[3]){
                        u0[cell][cellIndex] = 10.0;
                    }
                    else{
                        u0[cell][cellIndex] = 0.0;
                    }
                }
            }
        }
    } else {

        #pragma omp parallel num_threads(NUM_CELLS)
        {          
            int cell = omp_get_thread_num();  
  
            for (int i = 0; i < cellSizeY; i++){
                for (int j = 0; j < cellSizeX; j++){
                    u0[cell][i * cellSizeX + j] = 0.0;
                }
            }
        }
    }


    scheme->init();
}

void Diffusion::doCycle(const double dt)
{
    scheme->doAdvance(dt);
}

#include "VtkWriter.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <omp.h>

#include "DIVISIONS.h"

VtkWriter::VtkWriter(std::string basename, Mesh* mesh) :
    dump_basename(basename),
    vtk_header("# vtk DataFile Version 3.0\nvtk output\nASCII\n"),
    mesh(mesh)
{
    std::ofstream file;
    std::stringstream fname;

    fname << dump_basename << ".visit";

    std::string file_name = fname.str();

    file.open(file_name.c_str());

    file << "!NBLOCKS "
        << 1 << std::endl;
}

void VtkWriter::writeVisit(int stepMax)
{
    // Master process writes out the .visit file to coordinate the .vtk files
    std::ofstream file;
    std::stringstream fname;

    fname << dump_basename << ".visit";

    std::string file_name = fname.str();

    // Open file in append mode
    file.open(file_name.c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::app);

    for (int step = 0; step <= stepMax; step++){
        file << dump_basename
            << "." 
            << step 
            << "." 
            << 1 
            << ".vtk" << std::endl;
    }
}

void VtkWriter::writeVtk(int step, double time)
{
    std::ofstream file;
    std::stringstream fname;

    fname << dump_basename 
        << "." 
        << step 
        << "."
        << 1
        << ".vtk";

    std::string file_name = fname.str();

    file.open(file_name.c_str());

    file.setf(std::ios::fixed, std::ios::floatfield);
    file.precision(8);

    file << vtk_header;

    file << "DATASET RECTILINEAR_GRID" << std::endl;
    file << "FIELD FieldData 2" << std::endl;
    file << "TIME 1 1 double" << std::endl;
    file << time << std::endl;
    file << "CYCLE 1 1 int" << std::endl;
    file << step << std::endl;
    if (mesh->getDim() == 2) {
        file << "DIMENSIONS " << mesh->getNx()[0]+1
            << " " << mesh->getNx()[1]+1
            << " 1" << std::endl;
    } else if (mesh->getDim() == 3) {
        file << "DIMENSIONS " << mesh->getNx()[0]+1
            << " " << mesh->getNx()[1]+1
            << " " << mesh->getNx()[2]+1 << std::endl;
    }

    file << "X_COORDINATES " << mesh->getNx()[0]+1 << " float" << std::endl;
    for(int i = 1; i <= mesh->getNx()[0]+1; i++) {
        file << mesh->getPosGlobalX()[i] << " ";
    }

    file << std::endl;

    file << "Y_COORDINATES " << mesh->getNx()[1]+1 << " float" << std::endl;
    for(int j = 1; j <= mesh->getNx()[1]+1; j++) {
        file << mesh->getPosGlobalY()[j] << " ";
    }

    file << std::endl;

    file << "Z_COORDINATES 1 float" << std::endl;
    file << "0.0000" << std::endl;

    file << "CELL_DATA " << mesh->getNx()[0] * mesh->getNx()[1] << std::endl;

    file << "FIELD FieldData 1" << std::endl;

    file << "u 1 " << mesh->getNx()[0]*mesh->getNx()[1] << " double" << std::endl;

    int cellSizeX = mesh->getCellSize()[0];
    int cellSizeY = mesh->getCellSize()[1];

    double** uX = mesh->getUX(step);

    //loop math to output in order despite cells
    for (int cellY = 0; cellY < DIVISIONS; cellY++){
        for (int i = 1; i < cellSizeY-1; i++){
            for (int cellX = 0; cellX < DIVISIONS; cellX++){
                for (int j = 1; j < cellSizeX-1; j++){
                    file << uX[cellY * DIVISIONS + cellX][j + i*(cellSizeX)] << " ";
                }    
            }
            file << std::endl;
        }
    }

    file.close();
}

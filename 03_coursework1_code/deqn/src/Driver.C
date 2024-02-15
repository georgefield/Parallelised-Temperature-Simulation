#include "Driver.h"
#include <omp.h>
#include <time.h>
#include <iostream>

#include "DIVISIONS.h"


    Driver::Driver(const InputFile* input, const std::string& pname)
: problem_name(pname)
{

    std::cout << "+++++++++++++++++++++" << std::endl;
    std::cout << "  Running deqn v0.1  " << std::endl;
#ifdef DEBUG
    std::cout << "- input file: " << problem_name << std::endl;
#endif

    dt_max = input->getDouble("dt_max",  0.2);
    dt = input->getDouble("initial_dt", 0.2);

    t_start = input->getDouble("start_time", 0.0);
    t_end = input->getDouble("end_time", 10.0);

    vis_frequency = input->getInt("vis_frequency",-1);
    summary_frequency = input->getInt("summary_frequency", 1);

#ifdef DEBUG
    std::cout << "- dt_max: " << dt_max << std::endl;
    std::cout << "- initial_dt: " << dt << std::endl;
    std::cout << "- start_time: " << t_start << std::endl;
    std::cout << "- end_time: " << t_end << std::endl;
    std::cout << "- vis_frequency: " << vis_frequency << std::endl;
    std::cout << "- summary_frequency: " << summary_frequency << std::endl;
#endif
    std::cout << "+++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    mesh = new Mesh(input);
    diffusion = new Diffusion(input, mesh);
    writer = new VtkWriter(pname, mesh);

    /* Initial mesh dump */
    if(vis_frequency != -1)
        writer->writeVtk(0, 0.0);
}

Driver::~Driver() {
    delete mesh;
    delete diffusion;
    delete writer;
}

void Driver::writeVtks(int start, int count)
{
    double t_current;
    int step;

    #pragma omp parallel for private(t_current, step)
    for(int i = start; i < start + count; i++) {

        t_current = t_start + (i * dt);
        step = (t_current/dt) + 1 + 0.01;

        if(step % vis_frequency == 0 && vis_frequency != -1)
            writer->writeVtk(step, t_current);
    }

}

void Driver::run() {

    int step = 0;
    int count = 0;
    double t_current;
    for(t_current = t_start; t_current + (dt/2.0) < t_end; t_current += dt) { //+(dt/2.0) to stop floating point errors causing extra loop

        step = (t_current/dt) + 1 + 0.01; //0.01 to stop floating point errors causing cast down to wrong step number
        std::cout << "+ step: " << step << ", dt:   " << dt << std::endl;

        diffusion->doCycle(dt);

        if(step % summary_frequency == 0 && summary_frequency != -1) {
            double temperature = mesh->getTotalTemperature();
            std::cout << "+\tcurrent total temperature: " << temperature << std::endl;
        }

        count++;

        //save memory attempt (failed)
        /*if (count % omp_get_max_threads() == 0){
            double writeStartTime = omp_get_wtime();
            
            writeVtks(count - omp_get_max_threads(), omp_get_max_threads());

            writeTime += omp_get_wtime() - writeStartTime;
        }*/
    }


    writeVtks(0, count);
    writer->writeVisit(count);

    std::cout << std::endl;
    std::cout << "+++++++++++++++++++++" << std::endl;
    std::cout << "   Run completete.   " << std::endl;
    std::cout << "+++++++++++++++++++++" << std::endl;
}

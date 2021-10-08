#include "INMOST_ICE.h"

using namespace std;
using namespace INMOST;

int main(int argc, char* argv[])
{
#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // parse input line
    std::string current_exec_name = argv[0]; 
    std::vector<std::string> all_args;

    std::string config_json_path;

    if (argc > 1) 
    {
        all_args.assign(argv + 1, argv + argc);
    }

    if (all_args.size() > 1)
    {
        INMOST_ICE_ERR("should be only *.json file on input");
    }
    else
    {
        config_json_path = all_args.back();
    }

    // parse all inputs
    ConfigParams params(config_json_path);

    // log input
    
    if (rank == 0)
    {
        params.GetMeshParams().Log();
        params.GetModelParams().Log();
        params.GetAdvectionParams().Log();
        params.GetOutputParams().Log();
        params.GetMomentumParams().Log();
        params.GetIceForcingParams().Log();
        params.GetWaterForcingParams().Log();
        params.GetAirForcingParams().Log();
        //params.GetDynamicsTestParams().Log();
    }
    BARRIER
    

    // ice mesh initialization
    IceMesh arctic_mesh(params.GetMeshParams(),
                        params.GetOutputParams(),
                        params.GetModelParams());
    

    //Inititalization of Dynamics Test Class
    //DynamicsTest dynamics_test(arctic_mesh,
    //                           params.GetDynamicsTestParams(),
    //                           params.GetModelParams());

    // Assign initial scalars: h, a, m
    //dynamics_test.AssignInitialScalars();

    // Assign initial vectors: u_ice
    //dynamics_test.AssignInitialVectors();

    // Update water velocity: u_water 
    //dynamics_test.UpdateWaterVelocity();

    


    // Inititlize linear system solver for advection
    Solver linear_solver("inner_ilu2"); 
    linear_solver.SetParameter("absolute_tolerance", "1e-9");

    // Inititalization of mass and concentration advection class
    AdvectionSolver mass_transport(
            arctic_mesh, 
            params.GetAdvectionParams(),
            params.GetModelParams(),
            params.GetMeshParams(),
            ModelVariableNameToNotation["ice mass"],
            ModelVariableNameToNotation["ice velocity"],
            linear_solver,
            false);

    AdvectionSolver concentration_transport(
            arctic_mesh, 
            params.GetAdvectionParams(),
            params.GetModelParams(),
            params.GetMeshParams(),
            ModelVariableNameToNotation["ice concentration"],
            ModelVariableNameToNotation["ice velocity"],
            linear_solver,
            false);

    // Assembling LHS for transport solvers
    mass_transport.AssembleLHS();
    concentration_transport.AssembleLHS();

    // Inititalization of momentum balance solver
    MomentumSolver momentum_solver(
            arctic_mesh,
            params.GetMomentumParams(),
            params.GetModelParams(),
            params.GetMeshParams(),
            linear_solver,
            true);

    // time stepping parameters
    double current_time_hours = 0.0;
    size_t total_step_num = size_t(params.GetModelParams().GetTotalTimeHours()/params.GetModelParams().GetTimeStepHours());

    // log number of steps
    if (arctic_mesh.GetMesh()->GetProcessorRank() == 0)
    {
        cout << "Total time steps: " << total_step_num << endl;
    }

    // initialize ice forcing
    IceForcing ice_forcing(arctic_mesh,
                           params.GetIceForcingParams(),
                           params.GetModelParams()
                           );
    
    // initialize water forcing
    WaterForcing water_forcing(arctic_mesh,
                               params.GetWaterForcingParams(),
                               params.GetModelParams()
                               );
    
    // initialize air forcing
    AirForcing air_forcing(arctic_mesh,
                           params.GetAirForcingParams(),
                           params.GetModelParams()
                          );
    
    
    // log number of model time steps during ice forcing time gap
    if (arctic_mesh.GetMesh()->GetProcessorRank() == 0)
    {
        cout << "Model time steps in one ice forcing time gap: " << ice_forcing.GetNumOcurrences() << endl;
    }

    // log number of model time steps during water forcing time gap
    if (arctic_mesh.GetMesh()->GetProcessorRank() == 0)
    {
        cout << "Model time steps in one water forcing time gap: " << water_forcing.GetNumOcurrences() << endl;
    }

    // log number of model time steps during air forcing time gap
    if (arctic_mesh.GetMesh()->GetProcessorRank() == 0)
    {
        cout << "Model time steps in one air forcing time gap: " << air_forcing.GetNumOcurrences() << endl;
    }

    // interpolate initial scalars
    ice_forcing.UpdateScalars(0);
    water_forcing.UpdateScalars(0);
    air_forcing.UpdateScalars(0);

    // interpolate initial vectors
    ice_forcing.UpdateVectors(0);
    water_forcing.UpdateVectors(0);
    air_forcing.UpdateVectors(0);

    //// zero step without transport
    //dynamics_test.UpdateAirVelocity(current_time_hours);
    //momentum_solver.Evaluate();
    //
    //stringstream ss;
    //ss << setfill('0') << setw(5) << 0;
    //arctic_mesh.PrintPVTU("Step"+ ss.str() + ".pvtu");
    //current_time_hours += params.GetModelParams().GetTimeStepHours();
//
    //// calculate screenshots counter
    //size_t screenshots_counter = max(total_step_num/params.GetOutputParams().GetNumberOfScreenshots(), size_t(1));
    //
    //// time stepping
    //for (size_t step_num = 1; step_num < total_step_num; ++step_num)
    //{
    //    // increase current time
    //    current_time_hours += params.GetModelParams().GetTimeStepHours();
//
    //    // update wind speed
    //    dynamics_test.UpdateAirVelocity(current_time_hours);
    //    BARRIER
//
    //    // concenttantion transport
    //    concentration_transport.AssembleRHS();
    //    concentration_transport.Evaluate();
    //    BARRIER
//
    //    // mass transport
    //    mass_transport.AssembleRHS();
    //    mass_transport.Evaluate();
    //    BARRIER
//
    //    // momentum solver step
    //    momentum_solver.Evaluate();
//
    //    // print mesh
    //    if (step_num%screenshots_counter == 0)
    //    {
    //        stringstream sss;
    //        sss << setfill('0') << setw(5) << step_num;
    //        arctic_mesh.PrintPVTU("Step"+ sss.str() + ".pvtu");
    //    }
    //}
    arctic_mesh.PrintPVTU("Test.pvtu");
}
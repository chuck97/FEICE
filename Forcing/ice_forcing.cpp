#include "ice_forcing.h"

using namespace std;
using namespace INMOST;

IceForcing::IceForcing(IceMesh& im,
                       IceForcingParams& ifp,
                       ModelParams& mp):
                       ice_mesh(im),
                       ice_forcing_params(ifp),
                       model_params(mp)
{
    double forcing_time_gap = ice_forcing_params.GetTimeDiff();
    double model_time_step = model_params.GetTimeStepHours();
    if ((forcing_time_gap - model_time_step*(int)(forcing_time_gap/model_time_step)) > 1e-3)
    {
        INMOST_ICE_ERR("Ice forcing time gap should be consistent with model time step");
    }
    else
    {
        num_occur = (int)round(forcing_time_gap/model_time_step);
    }
    BARRIER
};

int IceForcing::GetNumOcurrences()
{
    return num_occur;
}

void IceForcing::UpdateScalars(int serial_number)
{
    NetcdfCoordNames coords = ice_forcing_params.GetNetcdfCoordNames();

    for (auto& [key, val]: ice_forcing_params.GetModelVarToFileVar())
    {
        if (key == ModelVariableNotation::u_ice)
        {
            continue;
        }
        else
        {
            if (key == ModelVariableNotation::a)
            {
                INMOST::Tag tag_to_interpolate = ice_mesh.GetData().NodeData[key];
                string nc_variable_name = val.nc_name[0];
                TopazScalarInterpolation(ice_mesh,
                                         tag_to_interpolate,
                                         nc_variable_name,
                                         ice_forcing_params.GetFilePath(),
                                         coords.x_name,
                                         coords.y_name,
                                         serial_number,
                                         val.scale_factor_name,
                                         val.invalid_value_name,
                                         val.offset_name,
                                         0.0,
                                         0.0,
                                         1.1,
                                         false);
            }
            else
            {
                INMOST::Tag tag_to_interpolate = ice_mesh.GetData().NodeData[key];
                string nc_variable_name = val.nc_name[0];
                TopazScalarInterpolation(ice_mesh,
                                         tag_to_interpolate,
                                         nc_variable_name,
                                         ice_forcing_params.GetFilePath(),
                                         coords.x_name,
                                         coords.y_name,
                                         serial_number,
                                         val.scale_factor_name,
                                         val.invalid_value_name,
                                         val.offset_name,
                                         0.0,
                                         0.0,
                                         numeric_limits<double>::max(),
                                         false);
            }
        }
        BARRIER
    }
};

void IceForcing::UpdateVectors(int serial_number)
{
    NetcdfCoordNames coords = ice_forcing_params.GetNetcdfCoordNames();
    for (auto& [key, val]: ice_forcing_params.GetModelVarToFileVar())
    {
        if (key != ModelVariableNotation::u_ice)
        {
            continue;
        }
        else
        {
            INMOST::Tag tag_to_interpolate = ice_mesh.GetData().NodeData[key];
            string nc_variable_name1 = val.nc_name[0];
            string nc_variable_name2 = val.nc_name[1];
            TopazVectorInterpolation(ice_mesh,
                                     tag_to_interpolate,
                                     nc_variable_name1,
                                     nc_variable_name2,
                                     ice_forcing_params.GetFilePath(),
                                     coords.x_name,
                                     coords.y_name,
                                     coords.lon_name,
                                     coords.lat_name,
                                     serial_number,
                                     val.scale_factor_name,
                                     val.invalid_value_name,
                                     val.offset_name,
                                     0.0,
                                     0.0,
                                     numeric_limits<double>::max(),
                                     false);        
        }
        BARRIER
    }
};
#include "air_forcing.h"

using namespace std;
using namespace INMOST;

AirForcing::AirForcing(IceMesh& im,
                       AirForcingParams& afp,
                       ModelParams& mp):
                       ice_mesh(im),
                       air_forcing_params(afp),
                       model_params(mp)
{
    double forcing_time_gap = air_forcing_params.GetTimeDiff();
    double model_time_step = model_params.GetTimeStepHours();
    if ((forcing_time_gap - model_time_step*(int)(forcing_time_gap/model_time_step)) > 1e-3)
    {
        INMOST_ICE_ERR("Air forcing time gap should be consistent with model time step");
    }
    else
    {
        num_occur = (int)round(forcing_time_gap/model_time_step);
    }
    BARRIER
};

int AirForcing::GetNumOcurrences()
{
    return num_occur;
}

void AirForcing::UpdateScalars(int serial_number)
{
    NetcdfCoordNames coords = air_forcing_params.GetNetcdfCoordNames();

    for (auto& [key, val]: air_forcing_params.GetModelVarToFileVar())
    {
        if (key == ModelVariableNotation::u_air)
        {
            continue;
        }
        else if (key == ModelVariableNotation::t_air)
        { 
            INMOST::Tag tag_to_interpolate = ice_mesh.GetData().NodeData[key];
            string nc_variable_name = val.nc_name[0];
            CamsScalarInterpolation(ice_mesh,
                                    tag_to_interpolate,
                                    nc_variable_name,
                                    air_forcing_params.GetFilePath(),
                                    coords.x_name,
                                    coords.y_name,
                                    serial_number,
                                    val.scale_factor_name,
                                    val.invalid_value_name,
                                    val.offset_name,
                                    273.0,
                                    273.0,
                                    373.0,
                                    false);
        }
        else
        {
            INMOST::Tag tag_to_interpolate = ice_mesh.GetData().NodeData[key];
            string nc_variable_name = val.nc_name[0];
            CamsScalarInterpolation(ice_mesh,
                                    tag_to_interpolate,
                                    nc_variable_name,
                                    air_forcing_params.GetFilePath(),
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
        BARRIER
    }
};

void AirForcing::UpdateVectors(int serial_number)
{
    NetcdfCoordNames coords = air_forcing_params.GetNetcdfCoordNames();
    
    for (auto& [key, val]: air_forcing_params.GetModelVarToFileVar())
    {
        if (key != ModelVariableNotation::u_air)
        {
            continue;
        }
        else
        {
            INMOST::Tag tag_to_interpolate = ice_mesh.GetData().NodeData[key];
            string nc_variable_name1 = val.nc_name[0];
            string nc_variable_name2 = val.nc_name[1];
            CamsVectorInterpolation(ice_mesh,
                                    tag_to_interpolate,
                                    nc_variable_name1,
                                    nc_variable_name2,
                                    air_forcing_params.GetFilePath(),
                                    coords.lon_name,
                                    coords.lat_name,
                                    serial_number,
                                    val.scale_factor_name,
                                    val.invalid_value_name,
                                    val.offset_name,
                                    0.0,
                                    0.0,
                                    100.0,
                                    false);        
        }
        BARRIER
    }
};
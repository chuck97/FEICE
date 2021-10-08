#include "water_forcing.h"

using namespace std;
using namespace INMOST;

WaterForcing::WaterForcing(IceMesh& im,
                           WaterForcingParams& wfp,
                           ModelParams& mp):
                           ice_mesh(im),
                           water_forcing_params(wfp),
                           model_params(mp)
{
    double forcing_time_gap = water_forcing_params.GetTimeDiff();
    double model_time_step = model_params.GetTimeStepHours();
    if ((forcing_time_gap - model_time_step*(int)(forcing_time_gap/model_time_step)) > 1e-3)
    {
        INMOST_ICE_ERR("Water forcing time gap should be consistent with model time step");
    }
    else
    {
        num_occur = (int)round(forcing_time_gap/model_time_step);
    }
    BARRIER
};

int WaterForcing::GetNumOcurrences()
{
    return num_occur;
}

void WaterForcing::UpdateScalars(int serial_number)
{
    NetcdfCoordNames coords = water_forcing_params.GetNetcdfCoordNames();

    for (auto& [key, val]: water_forcing_params.GetModelVarToFileVar())
    {
        if (key == ModelVariableNotation::u_water)
        {
            continue;
        }
        else if (key == ModelVariableNotation::hw)
        {
            INMOST::Tag tag_to_interpolate = ice_mesh.GetData().NodeData[key];
            string nc_variable_name = val.nc_name[0];
            TopazScalarInterpolation(ice_mesh,
                                     tag_to_interpolate,
                                     nc_variable_name,
                                     water_forcing_params.GetFilePath(),
                                     coords.x_name,
                                     coords.y_name,
                                     serial_number,
                                     val.scale_factor_name,
                                     val.invalid_value_name,
                                     val.offset_name,
                                     0.0,
                                     0.0,
                                     3.0,
                                     false);
        }
        else
        {
            INMOST::Tag tag_to_interpolate = ice_mesh.GetData().NodeData[key];
            string nc_variable_name = val.nc_name[0];
            TopazScalarInterpolation(ice_mesh,
                                     tag_to_interpolate,
                                     nc_variable_name,
                                     water_forcing_params.GetFilePath(),
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

void WaterForcing::UpdateVectors(int serial_number)
{
    NetcdfCoordNames coords = water_forcing_params.GetNetcdfCoordNames();

    for (auto& [key, val]: water_forcing_params.GetModelVarToFileVar())
    {
        if (key != ModelVariableNotation::u_water)
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
                                     water_forcing_params.GetFilePath(),
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
                                     3.0,
                                     false);        
        }
        BARRIER
    }
};
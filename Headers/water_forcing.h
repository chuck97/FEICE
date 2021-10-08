#pragma once
#include "external.h"
#include "enum_classes.h"
#include "config.h"
#include "mesh.h"
#include "forcing_interpolation.h"

class WaterForcing
{
public:
    WaterForcing(IceMesh& im,
                 WaterForcingParams& wfp,
                 ModelParams& mp);
    
    void UpdateScalars(int serial_number);
    void UpdateVectors(int serial_number);
    int GetNumOcurrences();

private:
    int num_occur; 
    int serial_number = 0;
    IceMesh& ice_mesh;
    WaterForcingParams& water_forcing_params;
    ModelParams& model_params;
};
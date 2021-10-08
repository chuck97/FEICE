#pragma once
#include "external.h"
#include "enum_classes.h"
#include "config.h"
#include "mesh.h"
#include "forcing_interpolation.h"

class AirForcing
{
public:
    AirForcing(IceMesh& im,
               AirForcingParams& afp,
               ModelParams& mp);
    
    void UpdateScalars(int serial_number);
    void UpdateVectors(int serial_number);
    int GetNumOcurrences();

private:
    int num_occur; 
    IceMesh& ice_mesh;
    AirForcingParams& air_forcing_params;
    ModelParams& model_params;
};
#pragma once
#include "external.h"
#include "enum_classes.h"
#include "config.h"
#include "mesh.h"
#include "forcing_interpolation.h"

class IceForcing
{
public:
    IceForcing(IceMesh& im,
               IceForcingParams& ifp,
               ModelParams& mp);
    
    void UpdateScalars(int serial_number);
    void UpdateVectors(int serial_number);
    int GetNumOcurrences();

private:
    int num_occur; 
    IceMesh& ice_mesh;
    IceForcingParams& ice_forcing_params;
    ModelParams& model_params;
};
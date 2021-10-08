#pragma once 
#include "external.h"
#include <proj_api.h>

std::vector<double> from_geo_2_topaz(double lon, double lat);

std::vector<double> from_topaz_2_geo(double x, double y);


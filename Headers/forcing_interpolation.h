#pragma once
#include "external.h"
#include "string_trim.h"
#include "enum_classes.h"
#include "mesh.h"
#include "interpolation2d.h"
#include "coords_rotation.h"
#include "coords_transformation.h"
#include <netcdf.h>

void TopazScalarInterpolation(IceMesh& ice_mesh,
                              INMOST::Tag tag_to_interpolate,
                              const std::string& nc_variable_name,
                              const std::string& filename,
                              const std::string& xname,
                              const std::string& yname,
                              int serial_number,
                              const std::string& scale_factor_name,
                              const std::string& invalid_value_name,
                              const std::string& offset_name,
                              double invalid_value_fill,
                              double no_extrapolation_fill,
                              double max_abs_value,
                              bool is_depth);

void CamsScalarInterpolation(IceMesh& ice_mesh,
                             INMOST::Tag tag_to_interpolate,
                             const std::string& nc_variable_name,
                             const std::string& filename,
                             const std::string& xname,
                             const std::string& yname,
                             int serial_number,
                             const std::string& scale_factor_name,
                             const std::string& invalid_value_name,
                             const std::string& offset_name,
                             double invalid_value_fill,
                             double no_extrapolation_fill,
                             double max_abs_value,
                             bool is_depth);

void TopazVectorInterpolation(IceMesh& ice_mesh,
                              INMOST::Tag tag_to_interpolate,
                              const std::string& nc_vec_variable_name1,
                              const std::string& nc_vec_variable_name2,
                              const std::string& filename,
                              const std::string& xname,
                              const std::string& yname,
                              const std::string& lonname,
                              const std::string& latname,
                              int serial_number,
                              const std::string& scale_factor_name,
                              const std::string& invalid_value_name,
                              const std::string& offset_name,
                              double invalid_value_fill,
                              double no_extrapolation_fill,
                              double max_abs_value,
                              bool is_depth
                              );

void CamsVectorInterpolation(IceMesh& ice_mesh,
                             INMOST::Tag tag_to_interpolate,
                             const std::string& nc_vec_variable_name1,
                             const std::string& nc_vec_variable_name2,
                             const std::string& filename,
                             const std::string& lonname,
                             const std::string& latname,
                             int serial_number,
                             const std::string& scale_factor_name,
                             const std::string& invalid_value_name,
                             const std::string& offset_name,
                             double invalid_value_fill,
                             double no_extrapolation_fill,
                             double max_abs_value,
                             bool is_depth
                             );
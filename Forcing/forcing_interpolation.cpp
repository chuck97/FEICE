#include "forcing_interpolation.h"

using namespace std;
using namespace INMOST;

vector<double> FindExtremalCoords(IceMesh& ice_mesh,
                                  NodeCoordsNotation coords_not)
{
    //Find extremal coords

    double min_topaz_x;
    double min_topaz_y;
    double max_topaz_x;
    double max_topaz_y;

    double tmp_min_x = ice_mesh.GetMesh()->BeginNode()->RealArray(ice_mesh.GetCoords()[coords_not])[0];
    double tmp_max_x = ice_mesh.GetMesh()->BeginNode()->RealArray(ice_mesh.GetCoords()[coords_not])[0];
    double tmp_min_y = ice_mesh.GetMesh()->BeginNode()->RealArray(ice_mesh.GetCoords()[coords_not])[1];
    double tmp_max_y = ice_mesh.GetMesh()->BeginNode()->RealArray(ice_mesh.GetCoords()[coords_not])[1];

    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                           nodeit != ice_mesh.GetMesh()->EndNode();
                           ++nodeit) 
	{
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double current_x = nodeit->RealArray(ice_mesh.GetCoords()[coords_not])[0];
            double current_y = nodeit->RealArray(ice_mesh.GetCoords()[coords_not])[1];
            if (current_x <= tmp_min_x)
            {
                tmp_min_x = current_x;
            }

            if (current_x >= tmp_max_x)
            {
                tmp_max_x = current_x;
            }

            if (current_y <= tmp_min_y)
            {
                tmp_min_y = current_y;
            }

            if (current_y >= tmp_max_y)
            {
                tmp_max_y = current_y;
            }
        }
    }
    min_topaz_x = tmp_min_x;
    min_topaz_y = tmp_min_y;
    max_topaz_x = tmp_max_x;
    max_topaz_y = tmp_max_y;

    BARRIER
    return {min_topaz_x, min_topaz_y, max_topaz_x, max_topaz_y};
};

vector<double> FindExtremalCamsCoords(IceMesh& ice_mesh,
                                      NodeCoordsNotation geo_not)
{
    //Find extremal coords

    double min_cams_x;
    double min_cams_y;
    double max_cams_x;
    double max_cams_y;

    double tmp_min_x = ice_mesh.GetMesh()->BeginNode()->RealArray(ice_mesh.GetCoords()[geo_not])[0]*(180.0/M_PI);
    double tmp_max_x = ice_mesh.GetMesh()->BeginNode()->RealArray(ice_mesh.GetCoords()[geo_not])[0]*(180.0/M_PI);
    double tmp_min_y = ice_mesh.GetMesh()->BeginNode()->RealArray(ice_mesh.GetCoords()[geo_not])[1]*(180.0/M_PI);
    double tmp_max_y = ice_mesh.GetMesh()->BeginNode()->RealArray(ice_mesh.GetCoords()[geo_not])[1]*(180.0/M_PI);

    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                           nodeit != ice_mesh.GetMesh()->EndNode();
                           ++nodeit) 
	{
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double current_x = nodeit->RealArray(ice_mesh.GetCoords()[geo_not])[0]*(180.0/M_PI);
            double current_y = nodeit->RealArray(ice_mesh.GetCoords()[geo_not])[1]*(180.0/M_PI);
            if (current_x <= tmp_min_x)
            {
                tmp_min_x = current_x;
            }

            if (current_x >= tmp_max_x)
            {
                tmp_max_x = current_x;
            }

            if (current_y <= tmp_min_y)
            {
                tmp_min_y = current_y;
            }

            if (current_y >= tmp_max_y)
            {
                tmp_max_y = current_y;
            }
        }
    }
    min_cams_x = tmp_min_x;
    min_cams_y = tmp_min_y;
    max_cams_x = tmp_max_x;
    max_cams_y = tmp_max_y;

    BARRIER
    return {min_cams_x, min_cams_y, max_cams_x, max_cams_y};
};

pair<size_t, double*> GetNetcdfCoords(int fileid,
                                      const string& coord_name)
{
    // Get coords from netcdf
    int retval;

    int coord_id;

    size_t coord_size;

    if((retval = nc_inq_dimid(fileid, coord_name.c_str(), &coord_id)))
        INMOST_ICE_ERR(to_string(retval));
    
    if((retval = nc_inq_dimlen(fileid, coord_id, &coord_size)))
        INMOST_ICE_ERR(to_string(retval));
    
    double* coords = new double[coord_size];

    int coords_id;

    if ((retval = nc_inq_varid(fileid, coord_name.c_str(), &coords_id)))
        INMOST_ICE_ERR(to_string(retval));

    if ((retval = nc_get_var_double(fileid, coords_id, coords)))
        INMOST_ICE_ERR(to_string(retval));
    
    BARRIER
    return {coord_size, coords};
};

vector<size_t> FindXStartCount(double* x_coords,
                               size_t x_size,
                               double dx,
                               double min_topaz_x,
                               double max_topaz_x)
{
    size_t x_start;
    size_t x_count;

    double max_x_begin = std::max(x_coords[0], min_topaz_x);
    double min_x_end = std::min(x_coords[x_size -1], max_topaz_x);

    if (max_x_begin > min_x_end)
    {
        x_count = 0;
        x_start = 0;
    }
    else if (max_x_begin == x_coords[0])
    {
        x_start = 0;
        x_count = std::min(size_t((max_topaz_x - x_coords[0])/dx) + 1, x_size);
    }
    else if (max_x_begin == min_topaz_x)
    {
        x_start = size_t((min_topaz_x - x_coords[0])/dx);
        x_count = std::min((size_t)((max_topaz_x - (x_coords[0] + dx*x_start))/dx) + 1, (x_size - x_start));
    }
    else
    {
        INMOST_ICE_ERR("bad x_start x_count procedure");
    }
    return {x_start, x_count};
};

vector<size_t> FindYStartCount(double* y_coords,
                               size_t y_size,
                               double dy,
                               double min_topaz_y,
                               double max_topaz_y)
{
    size_t y_start;
    size_t y_count;

    double max_y_begin = std::max(y_coords[0], min_topaz_y);
    double min_y_end = std::min(y_coords[y_size -1], max_topaz_y);

    if (max_y_begin > min_y_end)
    {
        y_count = 0;
        y_start = 0;
    }
    else if (max_y_begin == y_coords[0])
    {
        y_start = 0;
        y_count = std::min(size_t((max_topaz_y - y_coords[0])/dy) + 1, y_size);
    }
    else if (max_y_begin == min_topaz_y)
    {
        y_start = size_t((min_topaz_y - y_coords[0])/dy);
        y_count = std::min((size_t)((max_topaz_y - (y_coords[0] + dy*y_start))/dy) + 1, (y_size - y_start));
    }
    else
    {
        INMOST_ICE_ERR("bad y_start y_count procedure");
    }
    return{y_start, y_count};
};

double** GetScalarData(int fileid,
                       int serial_number,
                       const string& nc_variable_name,
                       const string& scale_factor_name,
                       const string& invalid_value_name,
                       const string& offset_name,
                       bool is_depth,
                       double invalid_value_fill,
                       size_t x_start,
                       size_t x_count,
                       size_t y_start,
                       size_t y_count
                      )
{
    int retval;

    // Find data id
    int data_variable_id;

    if ((retval = nc_inq_varid(fileid, nc_variable_name.c_str(), &data_variable_id)))
    {
        INMOST_ICE_ERR(to_string(retval));
    }

    //Get scale factor 
    double scale_factor;
    if (scale_factor_name.size() != 0)
    {
        if ((retval = nc_get_att(fileid, data_variable_id, scale_factor_name.c_str(), &scale_factor)))
        {
            INMOST_ICE_ERR(to_string(retval));
        }
    }
    else
    {
        scale_factor = 1.0;
    }
    
    //Get invalid value
    nc_type atttype;
    
    int invalid_value_int;
    short invalid_value_short;
    double invalid_value_double;

    if (invalid_value_name.size() != 0)
    {
        if ((retval = nc_inq_att(fileid, data_variable_id, invalid_value_name.c_str(), &atttype, NULL)))
            INMOST_ICE_ERR(to_string(retval));

        if (atttype == NC_SHORT)
        {
            if ((retval = nc_get_att(fileid, data_variable_id, invalid_value_name.c_str(), &invalid_value_short)))
                INMOST_ICE_ERR(to_string(retval));
        }
        else if ((atttype == NC_INT) or (atttype == NC_LONG))
        {
            if ((retval = nc_get_att(fileid, data_variable_id, invalid_value_name.c_str(), &invalid_value_int)))
                INMOST_ICE_ERR(to_string(retval));
        }
        else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
        {
            if ((retval = nc_get_att(fileid, data_variable_id, invalid_value_name.c_str(), &invalid_value_double)))
            {
                INMOST_ICE_ERR(to_string(retval));
            }
        }
        else
        {
            INMOST_ICE_ERR("att type is not int or double, can't read");
        }
    }
    else
    {
        invalid_value_int = std::nan("1");
        invalid_value_double = std::nan("1");
        invalid_value_short = std::nan("1");
    }

    //Get offset value
    double offset_value;
    if (offset_name.size() != 0)
    {
        if ((retval = nc_get_att(fileid, data_variable_id, offset_name.c_str(), &offset_value)))
        {
            INMOST_ICE_ERR(to_string(retval));
        }
    }
    else
    {
        offset_value = 0.0;
    }
    

    // Assemble start and count array
    size_t* start;
    size_t* count;

    if (is_depth)
    {
        start = new size_t[4];
        count = new size_t[4];
        start[0] = serial_number;
        start[1] = 0;
        start[2] = y_start;
        start[3] = x_start;
        count[0] = 1;
        count[1] = 1;
        count[2] = y_count;
        count[3] = x_count;
    }
    else
    {
        start = new size_t[3];
        count = new size_t[3];
        start[0] = serial_number;
        start[1] = y_start;
        start[2] = x_start;
        count[0] = 1;
        count[1] = y_count;
        count[2] = x_count;
    }

    // Get data in parallel

    int* data_local_int;
    short* data_local_short;
    double* data_local_double;

    if((retval = nc_inq_var(fileid, data_variable_id, NULL, &atttype, NULL, NULL, NULL)))
        INMOST_ICE_ERR(to_string(retval));

    
    if ((atttype == NC_SHORT))
    {
        data_local_short = new short[y_count*x_count];
        if ((retval = nc_get_vara_short(fileid, data_variable_id, start, count, data_local_short)))
            INMOST_ICE_ERR(to_string(retval));
    }
    else if ((atttype == NC_INT) or (atttype == NC_LONG))
    {
        data_local_int = new int[y_count*x_count];
        if ((retval = nc_get_vara_int(fileid, data_variable_id, start, count, data_local_int)))
            INMOST_ICE_ERR(to_string(retval));
    }
    else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
    {
        data_local_double = new double[y_count*x_count];
        if ((retval = nc_get_vara_double(fileid, data_variable_id, start, count, data_local_double)))
            INMOST_ICE_ERR(to_string(retval));
    }
    else
    {
        INMOST_ICE_ERR("var type is not int or double, can't read");
    }

    // convert data to double
    if (atttype == NC_SHORT)
    {
        data_local_double = new double[y_count*x_count];
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val;
                current_val = (data_local_short[j*x_count + i] == invalid_value_short) ? invalid_value_fill 
                                 : (double)(data_local_short[j*x_count + i])*scale_factor + offset_value;
                data_local_double[j*x_count + i] = current_val;
            }
        }
        delete[] data_local_short;
    }           
    else if ((atttype == NC_INT) or (atttype == NC_LONG))
    {
        data_local_double = new double[y_count*x_count];
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val;
                current_val = (data_local_int[j*x_count + i] == invalid_value_int) ? invalid_value_fill 
                                 : (double)(data_local_int[j*x_count + i])*scale_factor + offset_value;
                data_local_double[j*x_count + i] = current_val;
            }
        }
        delete[] data_local_int;
    }
    else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val;
                current_val = (data_local_double[j*x_count + i] == invalid_value_double) ? invalid_value_fill 
                                     : (data_local_double[j*x_count + i])*scale_factor + offset_value;
                data_local_double[j*x_count + i] = current_val;
            }
        }
    }
    else
    {
        INMOST_ICE_ERR("unknown value type");
    }
    BARRIER

    // convert 1d array to 2d
    double** data_local_double2d = new double*[y_count];
    
    for (size_t j = 0; j < y_count; ++j)
    {
        data_local_double2d[j] = new double[x_count];
        for(size_t i = 0; i < x_count; ++i)
        {
            data_local_double2d[j][i] = data_local_double[j*x_count + i];
        }
    }
    delete[] data_local_double;

    return  data_local_double2d;
};

void InterpolateScalarsTopaz(IceMesh& ice_mesh,
                             INMOST::Tag tag_to_interpolate,
                             double* x_coords,
                             double* y_coords,
                             double** data_local_double,
                             double no_extrapolation_fill,
                             double invalid_value_fill,
                             double max_abs_value,
                             size_t x_start,
                             size_t x_count,
                             size_t y_start,
                             size_t y_count,
                             double dx,
                             double dy)
{
    // get square variables     
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
        nodeit != ice_mesh.GetMesh()->EndNode();
        ++nodeit) 
    {
        double local_x_start = x_coords[x_start];
        double local_y_start = y_coords[y_start];
        double local_x_end = x_coords[x_start + x_count - 1];
        double local_y_end = y_coords[y_start + y_count - 1];

        if(nodeit->GetStatus() != Element::Ghost)
        {
            double x = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::topaz_stereographic])[0];
            double y = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::topaz_stereographic])[1];

            if ((x < local_x_start) or
                (x > local_x_end)   or
                (y < local_y_start) or
                (y > local_y_end))
            {
                nodeit->Real(tag_to_interpolate) = no_extrapolation_fill;
                continue;
            }
            else
            {
                size_t x_prev_pos = (size_t)((x - local_x_start)/dx);
                size_t y_prev_pos = (size_t)((y - local_y_start)/dy);
                
                double data_ld = data_local_double[y_prev_pos][x_prev_pos];
                double data_lu = data_local_double[y_prev_pos+1][x_prev_pos];
                double data_rd = data_local_double[y_prev_pos][x_prev_pos+1];
                double data_ru = data_local_double[y_prev_pos+1][x_prev_pos+1];
                
                double xl = x_coords[x_start + x_prev_pos];
                double xr = x_coords[x_start + x_prev_pos + 1];
                double yd = y_coords[y_start + y_prev_pos];
                double yu = y_coords[y_start + y_prev_pos + 1];

                double curr_data = bilinear_interpolation(xl, xr, 
		                                                  yd, yu, 
		                                                  x , y,
		                                                  data_ld, data_lu,
                                                          data_rd, data_ru);

				nodeit->Real(tag_to_interpolate) = 
                    (fabs(curr_data) > max_abs_value) ? invalid_value_fill : curr_data;
            }
        }
    }
    BARRIER;
};

void InterpolateScalarsCams(IceMesh& ice_mesh,
                            INMOST::Tag tag_to_interpolate,
                            double* x_coords,
                            double* y_coords,
                            double** data_local_double,
                            double no_extrapolation_fill,
                            double invalid_value_fill,
                            double max_abs_value,
                            size_t x_start,
                            size_t x_count,
                            size_t y_start,
                            size_t y_count,
                            double dx,
                            double dy)
{
    // get square variables     
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
        nodeit != ice_mesh.GetMesh()->EndNode();
        ++nodeit) 
    {
        double local_x_start = x_coords[x_start];
        double local_y_start = y_coords[y_start];
        double local_x_end = x_coords[x_start + x_count - 1];
        double local_y_end = y_coords[y_start + y_count - 1];

        if(nodeit->GetStatus() != Element::Ghost)
        {
            double x = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::geo])[0]*(180.0/M_PI);
            double y = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::geo])[1]*(180.0/M_PI);

            if ((x < local_x_start) or
                (x > local_x_end)   or
                (y < local_y_start) or
                (y > local_y_end))
            {
                nodeit->Real(tag_to_interpolate) = no_extrapolation_fill;
                continue;
            }
            else
            {
                size_t x_prev_pos = (size_t)((x - local_x_start)/dx);
                size_t y_prev_pos = (size_t)((y - local_y_start)/dy);
                
                double data_ld = data_local_double[y_prev_pos][x_prev_pos];
                double data_lu = data_local_double[y_prev_pos+1][x_prev_pos];
                double data_rd = data_local_double[y_prev_pos][x_prev_pos+1];
                double data_ru = data_local_double[y_prev_pos+1][x_prev_pos+1];
                
                double xl = x_coords[x_start + x_prev_pos];
                double xr = x_coords[x_start + x_prev_pos + 1];
                double yd = y_coords[y_start + y_prev_pos];
                double yu = y_coords[y_start + y_prev_pos + 1];

                double curr_data = bilinear_interpolation(xl, xr, 
		                                                  yd, yu, 
		                                                  x , y,
		                                                  data_ld, data_lu,
                                                          data_rd, data_ru);

				nodeit->Real(tag_to_interpolate) = 
                    (fabs(curr_data) > max_abs_value) ? invalid_value_fill : curr_data;
            }
        }
    }
    BARRIER;
}

void InterpolateVectorsTopaz(IceMesh& ice_mesh,
                             INMOST::Tag tag_to_interpolate,
                             double* x_coords,
                             double* y_coords,
                             vector<vector<double>> data_local_double1,
                             vector<vector<double>> data_local_double2,
                             double no_extrapolation_fill,
                             double invalid_value_fill,
                             double max_abs_value,
                             size_t x_start,
                             size_t x_count,
                             size_t y_start,
                             size_t y_count,
                             double dx,
                             double dy)
{
    // get square variables     
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
        nodeit != ice_mesh.GetMesh()->EndNode();
        ++nodeit) 
    {
        double local_x_start = x_coords[x_start];
        double local_y_start = y_coords[y_start];
        double local_x_end = x_coords[x_start + x_count - 1];
        double local_y_end = y_coords[y_start + y_count - 1];

        if(nodeit->GetStatus() != Element::Ghost)
        {
            double x = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::topaz_stereographic])[0];
            double y = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::topaz_stereographic])[1];

            if ((x < local_x_start) or
                (x > local_x_end)   or
                (y < local_y_start) or
                (y > local_y_end))
            {
                nodeit->RealArray(tag_to_interpolate)[0] = no_extrapolation_fill;
                nodeit->RealArray(tag_to_interpolate)[1] = no_extrapolation_fill;
                continue;
            }
            else
            {
                size_t x_prev_pos = (size_t)((x - local_x_start)/dx);
                size_t y_prev_pos = (size_t)((y - local_y_start)/dy);
                
                double data1_ld = data_local_double1[y_prev_pos][x_prev_pos];
                double data1_lu = data_local_double1[y_prev_pos+1][x_prev_pos];
                double data1_rd = data_local_double1[y_prev_pos][x_prev_pos+1];
                double data1_ru = data_local_double1[y_prev_pos+1][x_prev_pos+1];

                double data2_ld = data_local_double2[y_prev_pos][x_prev_pos];
                double data2_lu = data_local_double2[y_prev_pos+1][x_prev_pos];
                double data2_rd = data_local_double2[y_prev_pos][x_prev_pos+1];
                double data2_ru = data_local_double2[y_prev_pos+1][x_prev_pos+1];
                
                double xl = x_coords[x_start + x_prev_pos];
                double xr = x_coords[x_start + x_prev_pos + 1];
                double yd = y_coords[y_start + y_prev_pos];
                double yu = y_coords[y_start + y_prev_pos + 1];

                double curr_data1 = bilinear_interpolation(xl, xr, 
		                                                   yd, yu, 
		                                                   x , y,
		                                                   data1_ld, data1_lu,
                                                           data1_rd, data1_ru);
                
                double curr_data2 = bilinear_interpolation(xl, xr, 
		                                                   yd, yu, 
		                                                   x , y,
		                                                   data2_ld, data2_lu,
                                                           data2_rd, data2_ru);

				nodeit->RealArray(tag_to_interpolate)[0] = 
                    (fabs(curr_data1) > max_abs_value) ? invalid_value_fill : curr_data1;
                
                nodeit->RealArray(tag_to_interpolate)[1] = 
                    (fabs(curr_data2) > max_abs_value) ? invalid_value_fill : curr_data2;
            }
        }
    }
    BARRIER;
};

void InterpolateVectorsCams(IceMesh& ice_mesh,
                            INMOST::Tag tag_to_interpolate,
                            double* x_coords,
                            double* y_coords,
                            vector<vector<double>> data_local_double1,
                            vector<vector<double>> data_local_double2,
                            double no_extrapolation_fill,
                            double invalid_value_fill,
                            double max_abs_value,
                            size_t x_start,
                            size_t x_count,
                            size_t y_start,
                            size_t y_count,
                            double dx,
                            double dy)
{
    // get square variables     
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
        nodeit != ice_mesh.GetMesh()->EndNode();
        ++nodeit) 
    {
        double local_x_start = x_coords[x_start];
        double local_y_start = y_coords[y_start];
        double local_x_end = x_coords[x_start + x_count - 1];
        double local_y_end = y_coords[y_start + y_count - 1];

        if(nodeit->GetStatus() != Element::Ghost)
        {
            double x = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::geo])[0]*(180.0/M_PI);
            double y = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::geo])[1]*(180.0/M_PI);

            if ((x < local_x_start) or
                (x > local_x_end)   or
                (y < local_y_start) or
                (y > local_y_end))
            {
                nodeit->RealArray(tag_to_interpolate)[0] = no_extrapolation_fill;
                nodeit->RealArray(tag_to_interpolate)[1] = no_extrapolation_fill;
                continue;
            }
            else
            {
                size_t x_prev_pos = (size_t)((x - local_x_start)/dx);
                size_t y_prev_pos = (size_t)((y - local_y_start)/dy);
                
                double data1_ld = data_local_double1[y_prev_pos][x_prev_pos];
                double data1_lu = data_local_double1[y_prev_pos+1][x_prev_pos];
                double data1_rd = data_local_double1[y_prev_pos][x_prev_pos+1];
                double data1_ru = data_local_double1[y_prev_pos+1][x_prev_pos+1];

                double data2_ld = data_local_double2[y_prev_pos][x_prev_pos];
                double data2_lu = data_local_double2[y_prev_pos+1][x_prev_pos];
                double data2_rd = data_local_double2[y_prev_pos][x_prev_pos+1];
                double data2_ru = data_local_double2[y_prev_pos+1][x_prev_pos+1];
                
                double xl = x_coords[x_start + x_prev_pos];
                double xr = x_coords[x_start + x_prev_pos + 1];
                double yd = y_coords[y_start + y_prev_pos];
                double yu = y_coords[y_start + y_prev_pos + 1];

                double curr_data1 = bilinear_interpolation(xl, xr, 
		                                                   yd, yu, 
		                                                   x , y,
		                                                   data1_ld, data1_lu,
                                                           data1_rd, data1_ru);
                
                double curr_data2 = bilinear_interpolation(xl, xr, 
		                                                   yd, yu, 
		                                                   x , y,
		                                                   data2_ld, data2_lu,
                                                           data2_rd, data2_ru);

				nodeit->RealArray(tag_to_interpolate)[0] = 
                    (fabs(curr_data1) > max_abs_value) ? invalid_value_fill : curr_data1;
                
                nodeit->RealArray(tag_to_interpolate)[1] = 
                    (fabs(curr_data2) > max_abs_value) ? invalid_value_fill : curr_data2;
            }
        }
    }
    BARRIER;
};

void TopazScalarInterpolation(IceMesh& ice_mesh,
                              INMOST::Tag tag_to_interpolate,
                              const string& nc_variable_name,
                              const string& filename,
                              const string& xname,
                              const string& yname,
                              int serial_number,
                              const string& scale_factor_name,
                              const string& invalid_value_name,
                              const string& offset_name,
                              double invalid_value_fill,
                              double no_extrapolation_fill,
                              double max_abs_value,
                              bool is_depth)
{
    // open file
    int retval;
    int fileid;

    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &fileid)))
        INMOST_ICE_ERR(to_string(retval));

    //Initialize tag with no interpolation value
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                           nodeit != ice_mesh.GetMesh()->EndNode();
                           ++nodeit) 
	{
        nodeit->Real(tag_to_interpolate) = no_extrapolation_fill;
    }

    // find extremal coords
    vector<double> extremal_coords = FindExtremalCoords(ice_mesh, NodeCoordsNotation::topaz_stereographic);
    
    double min_topaz_x = extremal_coords[0];
    double min_topaz_y = extremal_coords[1];
    double max_topaz_x = extremal_coords[2];
    double max_topaz_y = extremal_coords[3];
    

    // get coords from netcdf
    auto x_coo = GetNetcdfCoords(fileid, xname);
    auto y_coo = GetNetcdfCoords(fileid, yname);

    size_t x_size = x_coo.first;
    size_t y_size = y_coo.first;

    double* x_coords = x_coo.second;
    double* y_coords = y_coo.second;
    
    double dx = x_coords[1] - x_coords[0];
    double dy = y_coords[1] - y_coords[0];

    // find x start & count
    vector<size_t> tmpx = FindXStartCount(x_coords, x_size, dx, min_topaz_x, max_topaz_x);
    size_t x_start = tmpx[0];
    size_t x_count = tmpx[1];

    // find y start & count
    vector<size_t> tmpy = FindYStartCount(y_coords, y_size, dy, min_topaz_y, max_topaz_y);
    size_t y_start = tmpy[0];
    size_t y_count = tmpy[1];


    // get data in parallel
    double** data_local_double = GetScalarData(fileid,
                                               serial_number,
                                               nc_variable_name,
                                               scale_factor_name,
                                               invalid_value_name,
                                               offset_name,
                                               is_depth,
                                               invalid_value_fill,
                                               x_start,
                                               x_count,
                                               y_start,
                                               y_count);
    

    // interpolate scalars
    InterpolateScalarsTopaz(ice_mesh,
                            tag_to_interpolate,
                            x_coords,
                            y_coords,
                            data_local_double,
                            no_extrapolation_fill,
                            invalid_value_fill,
                            max_abs_value,
                            x_start,
                            x_count,
                            y_start,
                            y_count,
                            dx,
                            dy);   

    // exchange data to ghost cells
    ice_mesh.GetMesh()->ExchangeData(tag_to_interpolate, NODE, 0);
    BARRIER

    // close file
    if ((retval = nc_close(fileid)))
        INMOST_ICE_ERR(to_string(retval));
    
    // clean data array
    for (size_t j = 0; j < y_count; ++j)
    {
        delete [] data_local_double[j]; 
    }
    BARRIER
};

void CamsScalarInterpolation(IceMesh& ice_mesh,
                             INMOST::Tag tag_to_interpolate,
                             const string& nc_variable_name,
                             const string& filename,
                             const string& xname,
                             const string& yname,
                             int serial_number,
                             const string& scale_factor_name,
                             const string& invalid_value_name,
                             const string& offset_name,
                             double invalid_value_fill,
                             double no_extrapolation_fill,
                             double max_abs_value,
                             bool is_depth)
{
    // open file
    int retval;
    int fileid;

    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &fileid)))
        INMOST_ICE_ERR(to_string(retval));

    //Initialize tag with no interpolation value
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                           nodeit != ice_mesh.GetMesh()->EndNode();
                           ++nodeit) 
	{
        nodeit->Real(tag_to_interpolate) = no_extrapolation_fill;
    }

    // find extremal coords
    vector<double> extremal_coords = FindExtremalCamsCoords(ice_mesh, NodeCoordsNotation::geo);
    
    double min_cams_x = extremal_coords[0];
    double min_cams_y = extremal_coords[1];
    double max_cams_x = extremal_coords[2];
    double max_cams_y = extremal_coords[3];

    //for (size_t i = 0; i < ice_mesh.GetMesh()->GetProcessorsNumber(); ++i)
    //{
    //    if (ice_mesh.GetMesh()->GetProcessorRank() == i)
    //    {
    //        cout << "pr num: " << ice_mesh.GetMesh()->GetProcessorRank() 
    //             << " min x: " << min_cams_x << " max x: " << max_cams_x << endl
    //             << " min y: " << min_cams_y << " max y: " << max_cams_y << endl; 
    //    }
    //    BARRIER
    //}
    

    // get coords from netcdf
    auto x_coo = GetNetcdfCoords(fileid, xname);
    auto y_coo = GetNetcdfCoords(fileid, yname);

    size_t x_size = x_coo.first;
    size_t y_size = y_coo.first;

    double* x_coords = x_coo.second;
    double* y_coords = y_coo.second;

    
    double dx = x_coords[1] - x_coords[0];
    double dy = y_coords[1] - y_coords[0];


    // find x start & count
    vector<size_t> tmpx = FindXStartCount(x_coords, x_size, dx, min_cams_x, max_cams_x);
    size_t x_start = tmpx[0];
    size_t x_count = tmpx[1];

    // find y start & count
    vector<size_t> tmpy = FindYStartCount(y_coords, y_size, dy, min_cams_y, max_cams_y);
    size_t y_start = tmpy[0];
    size_t y_count = tmpy[1];

    // get data in parallel
    double** data_local_double = GetScalarData(fileid,
                                               serial_number,
                                               nc_variable_name,
                                               scale_factor_name,
                                               invalid_value_name,
                                               offset_name,
                                               is_depth,
                                               invalid_value_fill,
                                               x_start,
                                               x_count,
                                               y_start,
                                               y_count);
    

    // interpolate scalars
    InterpolateScalarsCams(ice_mesh,
                           tag_to_interpolate,
                           x_coords,
                           y_coords,
                           data_local_double,
                           no_extrapolation_fill,
                           invalid_value_fill,
                           max_abs_value,
                           x_start,
                           x_count,
                           y_start,
                           y_count,
                           dx,
                           dy);   

    // exchange data to ghost cells
    ice_mesh.GetMesh()->ExchangeData(tag_to_interpolate, NODE, 0);
    BARRIER

    // close file
    if ((retval = nc_close(fileid)))
        INMOST_ICE_ERR(to_string(retval));
    
    // clean data array
    for (size_t j = 0; j < y_count; ++j)
    {
        delete [] data_local_double[j]; 
    }
    BARRIER
};

void TopazVectorInterpolation(IceMesh& ice_mesh,
                              INMOST::Tag tag_to_interpolate,
                              const string& nc_vec_variable_name1,
                              const string& nc_vec_variable_name2,
                              const string& filename,
                              const string& xname,
                              const string& yname,
                              const string& lonname,
                              const string& latname,
                              int serial_number,
                              const string& scale_factor_name,
                              const string& invalid_value_name,
                              const string& offset_name,
                              double invalid_value_fill,
                              double no_extrapolation_fill,
                              double max_abs_value,
                              bool is_depth
                              )
{
    // open file
    int retval;
    int fileid;

    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &fileid)))
        INMOST_ICE_ERR(to_string(retval));

    //Initialize tag with no interpolation value
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                           nodeit != ice_mesh.GetMesh()->EndNode();
                           ++nodeit) 
	{
        nodeit->RealArray(tag_to_interpolate)[0] = no_extrapolation_fill;
        nodeit->RealArray(tag_to_interpolate)[1] = no_extrapolation_fill;
        nodeit->RealArray(tag_to_interpolate)[2] = 0.0;
    }

    // find extremal coords
    vector<double> extremal_coords = FindExtremalCoords(ice_mesh, NodeCoordsNotation::topaz_stereographic);
    
    double min_topaz_x = extremal_coords[0];
    double min_topaz_y = extremal_coords[1];
    double max_topaz_x = extremal_coords[2];
    double max_topaz_y = extremal_coords[3];
    

    // get coords from netcdf
    auto x_coo = GetNetcdfCoords(fileid, xname);
    auto y_coo = GetNetcdfCoords(fileid, yname);

    size_t x_size = x_coo.first;
    size_t y_size = y_coo.first;

    double* x_coords = x_coo.second;
    double* y_coords = y_coo.second;
    
    double dx = x_coords[1] - x_coords[0];
    double dy = y_coords[1] - y_coords[0];

    // find x start & count
    vector<size_t> tmpx = FindXStartCount(x_coords, x_size, dx, min_topaz_x, max_topaz_x);
    size_t x_start = tmpx[0];
    size_t x_count = tmpx[1];

    // find y start & count
    vector<size_t> tmpy = FindYStartCount(y_coords, y_size, dy, min_topaz_y, max_topaz_y);
    size_t y_start = tmpy[0];
    size_t y_count = tmpy[1];

    // get data in parallel
    double** vector_data_local_double1 = GetScalarData(fileid,
                                                       serial_number,
                                                       nc_vec_variable_name1,
                                                       scale_factor_name,
                                                       invalid_value_name,
                                                       offset_name,
                                                       is_depth,
                                                       invalid_value_fill,
                                                       x_start,
                                                       x_count,
                                                       y_start,
                                                       y_count);
    
    // get data in parallel
    double** vector_data_local_double2 = GetScalarData(fileid,
                                                       serial_number,
                                                       nc_vec_variable_name2,
                                                       scale_factor_name,
                                                       invalid_value_name,
                                                       offset_name,
                                                       is_depth,
                                                       invalid_value_fill,
                                                       x_start,
                                                       x_count,
                                                       y_start,
                                                       y_count);

    // get lat lon
    int lon_id, lat_id;

    if ((retval = nc_inq_varid(fileid, lonname.c_str(), &lon_id)))
    {
        INMOST_ICE_ERR(to_string(retval));
    }

    if ((retval = nc_inq_varid(fileid, latname.c_str(), &lat_id)))
    {
        INMOST_ICE_ERR(to_string(retval));
    }

    double lon_data_local_double[y_count][x_count];
    double lat_data_local_double[y_count][x_count];

    size_t start_coords[] = {y_start, x_start};
    size_t count_coords[] = {y_count, x_count};

    if ((retval = nc_get_vara_double(fileid, lon_id, start_coords, count_coords, &lon_data_local_double[0][0])))
        INMOST_ICE_ERR(retval)

    if ((retval = nc_get_vara_double(fileid, lat_id, start_coords, count_coords, &lat_data_local_double[0][0])))
        INMOST_ICE_ERR(retval)

    // rotate vector to geo directions using lat, lon info
    vector<vector<double>> rotated_vector_x(y_count, vector<double>(x_count));
    vector<vector<double>> rotated_vector_y(y_count, vector<double>(x_count));

    double radian = M_PI/180.0;
    double dlon_ip1[y_count][x_count];
    double dlat_ip1[y_count][x_count];
    double theta_ip1;
    double theta_jp1;
    
    for (size_t j = 0; j < y_count; ++j)
    {
        for (size_t i = 0; i < x_count - 1; ++i)
        {
            dlon_ip1[j][i] = (lon_data_local_double[j][i+1] - lon_data_local_double[j][i])*
            cos(radian* 0.5*(lat_data_local_double[j][i+1] + lat_data_local_double[j][i]));

            dlat_ip1[j][i] = (lat_data_local_double[j][i+1] - lat_data_local_double[j][i]);
        }
        dlon_ip1[j][x_count - 1] = dlon_ip1[j][x_count - 2];
        dlat_ip1[j][x_count - 1] = dlat_ip1[j][x_count - 2];
    }

    for (size_t j = 0; j < y_count; ++j)
    {
        for (size_t i = 0; i < x_count; ++i)
        {
            theta_ip1 = atan2(dlat_ip1[j][i], dlon_ip1[j][i]);
            theta_jp1 = theta_ip1 + radian*90.0;

            rotated_vector_x[j][i] = vector_data_local_double1[j][i]*cos(theta_ip1) + 
                                     vector_data_local_double2[j][i]*cos(theta_jp1);
            rotated_vector_y[j][i] = vector_data_local_double1[j][i]*sin(theta_ip1) + 
                                     vector_data_local_double2[j][i]*sin(theta_jp1);
        }
    } 

    // rotate vector to model directions
    for (size_t j = 0; j < y_count; ++j)
    {
        for (size_t i = 0; i < x_count; ++i)
        {
            double vec_x = rotated_vector_x[j][i];
            double vec_y = rotated_vector_y[j][i];
            std::vector<double> model_rot_vec = from_geo_2_model_vec<double>(vec_x, vec_y,
                                                                             lon_data_local_double[j][i],
                                                                             lat_data_local_double[j][i]);
            rotated_vector_x[j][i] = model_rot_vec[0];
            rotated_vector_y[j][i] = model_rot_vec[1];
       }
    }



    // interpolate vectors
    InterpolateVectorsTopaz(ice_mesh,
                            tag_to_interpolate,
                            x_coords,
                            y_coords,
                            rotated_vector_x,
                            rotated_vector_y,
                            no_extrapolation_fill,
                            invalid_value_fill,
                            max_abs_value,
                            x_start,
                            x_count,
                            y_start,
                            y_count,
                            dx,
                            dy); 

    // exchange data to ghost cells
    ice_mesh.GetMesh()->ExchangeData(tag_to_interpolate, NODE, 0);

    // close file
    if ((retval = nc_close(fileid)))
        INMOST_ICE_ERR(retval);

    // clean data arrays
    for (size_t j = 0; j < y_count; ++j)
    {
        delete [] vector_data_local_double1[j];
        delete [] vector_data_local_double2[j]; 
    }
    BARRIER
};

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
                             )
{
    // open file
    int retval;
    int fileid;

    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &fileid)))
        INMOST_ICE_ERR(to_string(retval));

    //Initialize tag with no extrapolation value
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
                           nodeit != ice_mesh.GetMesh()->EndNode();
                           ++nodeit) 
	{
        nodeit->RealArray(tag_to_interpolate)[0] = no_extrapolation_fill;
        nodeit->RealArray(tag_to_interpolate)[1] = no_extrapolation_fill;
        nodeit->RealArray(tag_to_interpolate)[2] = 0.0;
    }

    // find extremal coords
    vector<double> extremal_coords = FindExtremalCamsCoords(ice_mesh, NodeCoordsNotation::geo);
    
    double min_cams_x = extremal_coords[0];
    double min_cams_y = extremal_coords[1];
    double max_cams_x = extremal_coords[2];
    double max_cams_y = extremal_coords[3];
    

    // get coords from netcdf
    auto x_coo = GetNetcdfCoords(fileid, lonname);
    auto y_coo = GetNetcdfCoords(fileid, latname);

    size_t x_size = x_coo.first;
    size_t y_size = y_coo.first;

    double* x_coords = x_coo.second;
    double* y_coords = y_coo.second;
    
    double dx = x_coords[1] - x_coords[0];
    double dy = y_coords[1] - y_coords[0];

    // find x start & count
    vector<size_t> tmpx = FindXStartCount(x_coords, x_size, dx, min_cams_x, max_cams_x);
    size_t x_start = tmpx[0];
    size_t x_count = tmpx[1];

    // find y start & count
    vector<size_t> tmpy = FindYStartCount(y_coords, y_size, dy, min_cams_y, max_cams_y);
    size_t y_start = tmpy[0];
    size_t y_count = tmpy[1];

    // get data in parallel
    double** vector_data_local_double1 = GetScalarData(fileid,
                                                       serial_number,
                                                       nc_vec_variable_name1,
                                                       scale_factor_name,
                                                       invalid_value_name,
                                                       offset_name,
                                                       is_depth,
                                                       invalid_value_fill,
                                                       x_start,
                                                       x_count,
                                                       y_start,
                                                       y_count);
    
    // get data in parallel
    double** vector_data_local_double2 = GetScalarData(fileid,
                                                       serial_number,
                                                       nc_vec_variable_name2,
                                                       scale_factor_name,
                                                       invalid_value_name,
                                                       offset_name,
                                                       is_depth,
                                                       invalid_value_fill,
                                                       x_start,
                                                       x_count,
                                                       y_start,
                                                       y_count);


    // rotate vector to model directions
    vector<vector<double>> rotated_vector_x(y_count, vector<double>(x_count));
    vector<vector<double>> rotated_vector_y(y_count, vector<double>(x_count));

    for (size_t j = 0; j < y_count; ++j)
    {
        for (size_t i = 0; i < x_count; ++i)
        {
            double vec_x = vector_data_local_double1[j][i];
            double vec_y = vector_data_local_double2[j][i];
            std::vector<double> model_rot_vec = from_geo_2_model_vec<double>(vec_x, vec_y,
                                                                             x_coords[x_start + i],
                                                                             y_coords[y_start + j]);
            rotated_vector_x[j][i] = model_rot_vec[0];
            rotated_vector_y[j][i] = model_rot_vec[1];
       }
    }

    // interpolate vectors
    InterpolateVectorsCams(ice_mesh,
                           tag_to_interpolate,
                           x_coords,
                           y_coords,
                           rotated_vector_x,
                           rotated_vector_y,
                           no_extrapolation_fill,
                           invalid_value_fill,
                           max_abs_value,
                           x_start,
                           x_count,
                           y_start,
                           y_count,
                           dx,
                           dy); 

    // exchange data to ghost cells
    ice_mesh.GetMesh()->ExchangeData(tag_to_interpolate, NODE, 0);

    // close file
    if ((retval = nc_close(fileid)))
        INMOST_ICE_ERR(retval);

    // clean data arrays
    for (size_t j = 0; j < y_count; ++j)
    {
        delete [] vector_data_local_double1[j];
        delete [] vector_data_local_double2[j]; 
    }
    BARRIER
};
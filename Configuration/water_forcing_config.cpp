#include "config.h"

using namespace std;

WaterForcingParams::WaterForcingParams()
{};

WaterForcingParams::WaterForcingParams(const std::string& json_path)
{
    string no_spaces_json = json_path;
    rtrim(no_spaces_json);
    if(no_spaces_json.substr(no_spaces_json.size()-5, 5) != ".json")
    {
        INMOST_ICE_ERR("input file shoud be ended by .json!");
    }

    ifstream ifs(no_spaces_json);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    BARRIER

    // Parse .nc file path
    if (!j_input["file path"].empty())
    {
        file_path = j_input["file path"];

        // check file extension
        if(file_path.substr(file_path.size()-3, 3) != ".nc")
        {
            INMOST_ICE_ERR("input file for water forcing should be ended by .nc!");
        }
    }
    else
    {
        INMOST_ICE_ERR(".nc file for water forcing should be given!");
    }
    BARRIER
    
    // Parse data time gap
    if (!j_input["time difference (h)"].empty())
    {
        time_diff = j_input["time difference (h)"];
    }
    else
    {
        INMOST_ICE_ERR("time difference for water forcing should be given!");
    }
    BARRIER

    // Parse netcdf coords
    if (!j_input["coordinates"].empty())
    {
        NetcdfCoordNames c;
        c.x_name = j_input["coordinates"]["x"];
        c.y_name = j_input["coordinates"]["y"];
        c.lon_name = j_input["coordinates"]["lon"];
        c.lat_name = j_input["coordinates"]["lat"];
        coords = c;
    }
    else
    {
        INMOST_ICE_ERR("netcdf coords for water forcing should be given!");
    }
    BARRIER

    // Parse model variables input list
    if (!j_input["variables"].empty())
    {
        if (!j_input["variables"]["water level"].empty())
        {
            NetcdfVarAtts v;
            v.nc_name.push_back(j_input["variables"]["water level"]["nc name"]);
            v.scale_factor_name = j_input["variables"]["water level"]["scale factor name"];
            v.offset_name = j_input["variables"]["water level"]["offset name"];
            v.invalid_value_name = j_input["variables"]["water level"]["invalid value name"];

            model_var_to_file_var[ModelVariableNameToNotation["water level"]] = v;
        }

        if (!j_input["variables"]["water velocity"].empty())
        {
            NetcdfVarAtts v;
            v.nc_name.push_back(j_input["variables"]["water velocity"]["nc name"][0]);
            v.nc_name.push_back(j_input["variables"]["water velocity"]["nc name"][1]);
            v.scale_factor_name = j_input["variables"]["water velocity"]["scale factor name"];
            v.offset_name = j_input["variables"]["water velocity"]["offset name"];
            v.invalid_value_name = j_input["variables"]["water velocity"]["invalid value name"];

            model_var_to_file_var[ModelVariableNameToNotation["water velocity"]] = v;
        }
    }
    else
    {
        INMOST_ICE_ERR("at least one water forcing variable should be given!")
    }
    BARRIER
};

double WaterForcingParams::GetTimeDiff() const
{
    return time_diff;
};

std::string WaterForcingParams::GetFilePath() const
{
    return file_path;
};

std::map<ModelVariableNotation, NetcdfVarAtts> WaterForcingParams::GetModelVarToFileVar() const
{
    return model_var_to_file_var;
};

NetcdfCoordNames WaterForcingParams::GetNetcdfCoordNames() const
{
    return coords;
}

void WaterForcingParams::Log() const
{
    cout << endl;
    cout << "#### Water Forcing config params ####" << endl;
    cout << "Water forcing file: " << file_path << endl;
    cout << "Time gap in data (h): " << time_diff << endl;
    cout << "data: ";
    for (const auto &[key, value]: model_var_to_file_var)
    {
        cout << ModelVariableNotationToName[key] << ", ";
    }
    cout << endl;
    cout << "###################################" << endl;
}
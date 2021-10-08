#include "config.h"

using namespace std;

AirForcingParams::AirForcingParams()
{};

AirForcingParams::AirForcingParams(const std::string& json_path)
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
            INMOST_ICE_ERR("input file for air forcing should be ended by .nc!");
        }
    }
    else
    {
        INMOST_ICE_ERR(".nc file for air forcing should be given!");
    }
    BARRIER
    
    // Parse data time gap
    if (!j_input["time difference (h)"].empty())
    {
        time_diff = j_input["time difference (h)"];
    }
    else
    {
        INMOST_ICE_ERR("time difference for air forcing should be given!");
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
        INMOST_ICE_ERR("netcdf coords for air forcing should be given!");
    }
    BARRIER

    // Parse model variables input list
    if (!j_input["variables"].empty())
    {
        if (!j_input["variables"]["air velocity"].empty())
        {
            NetcdfVarAtts v;
            v.nc_name.push_back(j_input["variables"]["air velocity"]["nc name"][0]);
            v.nc_name.push_back(j_input["variables"]["air velocity"]["nc name"][1]);
            v.scale_factor_name = j_input["variables"]["air velocity"]["scale factor name"];
            v.offset_name = j_input["variables"]["air velocity"]["offset name"];
            v.invalid_value_name = j_input["variables"]["air velocity"]["invalid value name"];

            model_var_to_file_var[ModelVariableNameToNotation["air velocity"]] = v;
        }

        if (!j_input["variables"]["air temperature"].empty())
        {
            NetcdfVarAtts v;
            v.nc_name.push_back(j_input["variables"]["air temperature"]["nc name"]);
            v.scale_factor_name = j_input["variables"]["air temperature"]["scale factor name"];
            v.offset_name = j_input["variables"]["air temperature"]["offset name"];
            v.invalid_value_name = j_input["variables"]["air temperature"]["invalid value name"];

            model_var_to_file_var[ModelVariableNameToNotation["air temperature"]] = v;
        }
    }
    else
    {
        INMOST_ICE_ERR("at least one air forcing variable should be given!")
    }
    BARRIER
};

double AirForcingParams::GetTimeDiff() const
{
    return time_diff;
};

std::string AirForcingParams::GetFilePath() const
{
    return file_path;
};

std::map<ModelVariableNotation, NetcdfVarAtts> AirForcingParams::GetModelVarToFileVar() const
{
    return model_var_to_file_var;
};

NetcdfCoordNames AirForcingParams::GetNetcdfCoordNames() const
{
    return coords;
}

void AirForcingParams::Log() const
{
    cout << endl;
    cout << "#### Air Forcing config params ####" << endl;
    cout << "Air forcing file: " << file_path << endl;
    cout << "Time gap in data (h): " << time_diff << endl;
    cout << "data: ";
    for (const auto &[key, value]: model_var_to_file_var)
    {
        cout << ModelVariableNotationToName[key] << ", ";
    }
    cout << endl;
    cout << "###################################" << endl;
}
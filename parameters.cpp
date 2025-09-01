//
// Created by Peter Berg Ammundsen on 06/08/2025.
//

#include "parameters.h"

// Access to the application's global console (defined in main.cpp)
class QPlainTextEdit;
extern QPlainTextEdit* g_consoleOutput;

#include <nlohmann/json.hpp>
void parameters::load_from_json(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open JSON file: " + filename);
    }
    using json = nlohmann::json;
    json data = json::parse(file);

    if (data.contains("nx")) nx = data["nx"];
    if (data.contains("ny")) ny = data["ny"];
    if (data.contains("xmin")) xmin = data["xmin"];
    if (data.contains("xmax")) xmax = data["xmax"];
    if (data.contains("ymin")) ymin = data["ymin"];
    if (data.contains("ymax")) ymax = data["ymax"];
    if (data.contains("rho")) rho = data["rho"];
    if (data.contains("lid_velocity")) lid_velocity = data["lid_velocity"];
    if (data.contains("Reynolds_number")) Re = data["Reynolds_number"];
    if (data.contains("u_iter_max")) u_iter_max = data["u_iter_max"];
    if (data.contains("v_iter_max")) v_iter_max = data["v_iter_max"];
    if (data.contains("pressure_iter_max")) pressure_iter_max = data["pressure_iter_max"];
    if (data.contains("beta")) beta = data["beta"];
    if (data.contains("alpha_u")) alpha_u = data["alpha_u"];
    if (data.contains("alpha_v")) alpha_v = data["alpha_v"];
    if (data.contains("alpha_p")) alpha_p = data["alpha_p"];
    if (data.contains("global_tol")) global_tol = data["global_tol"];
    int peter = 0;
}
//
// Created by Peter Berg Ammundsen on 06/08/2025.
//

#ifndef PARAMETER_H
#define PARAMETER_H
#include <vector>
#include <fstream>

/// \brief class to handle all parameters
///
/// Both general and specify fliud paramer
class parameters {

public:
    /// \brief input parameter

    ///geometry parameters
    /// Max number of points in x-direction
    unsigned int nx = 0.0;
    /// Max number of points in y-direction
    unsigned int ny = 0.0;
    ///Domain X min value
    double xmin = 0.0;
    ///Domain X max value
    double xmax = 0.0;
    ///Domain X min value
    double ymin = 0.0;
    ///Domain X max value
    double ymax = 0.0;    
    ///  \brief Fluid parameters
    /// Density of fluid[kg/m³]
    double rho = 0.0;
    /// Top plate velocity
    double lid_velocity = 0.0;
    /// Reynolds number
    double Re = 0.0;

    /// Max iterations of u-momentum solver
    int u_iter_max = 0;
    /// Max iterations of v-momentum solver
    int v_iter_max = 0;
    /// Max iterations of pressure correction solver
    int pressure_iter_max = 0;
    /// Set convection scheme blending factor
    double beta = 0.0;

    ///U-momentum relaxation factor
    double alpha_u = 0.0;
    ///V-momentum relaxation factor
    double alpha_v = 0.0;
    ///Pressure correction relaxation factor
    double alpha_p = 0.0;
    /// Tolerance to converge SIMPLE algorithm
    double global_tol = 0.0;

    //Generate Parameter
    /// Grid spacing x-dirction
    double dx = 0;
    /// Grid spacing y-dirction
    double dy = 0;
    /// Dynamic viscosity mu [m²/s]
    double mu = 0.0;

    /// U velocities in mesh
    std::vector<std::vector<double>> u;
    /// V velocities in mesh
    std::vector<std::vector<double>> v;
    /// Pressure in Mesh
    std::vector<std::vector<double>> p;
    ///
    std::vector<std::vector<double>> pprime;
    std::vector<std::vector<double>> pprime_old;
    std::vector<std::vector<double>> u_intp;
    std::vector<std::vector<double>> v_intp;
    std::vector<std::vector<double>> p_intp;
    std::vector<std::vector<double>> ap_u;
    std::vector<std::vector<double>> ap_v;
    std::vector<std::vector<double>> ap_p;
    std::vector<std::vector<double>> ae;
    std::vector<std::vector<double>> aw;
    std::vector<std::vector<double>> as;
    std::vector<std::vector<double>> an;

    std::vector<std::vector<double>> upw;
    std::vector<std::vector<double>> cds;
    std::vector<std::vector<double>> v_old;
    std::vector<std::vector<double>> u_old;
    std::vector<std::vector<double>> b_mass_imb;

    // --- strings til pæn visning ---
    std::string u_s;
    std::string v_s;
    std::string p_s;
    std::string pprime_s, pprime_old_s;
    std::string u_intp_s, v_intp_s, p_intp_s;
    std::string ap_u_s, ap_v_s, ap_p_s;
    std::string ae_s, aw_s, as_s, an_s;
    std::string upw_s, cds_s;
    std::string v_old_s, u_old_s;
    std::string mass_imb_s;
    double u_residual = 0.0;
    double global_res = 0.0;

    //functions
    void load_from_json(const std::string& filename);
};



#endif //PARAMETER_H

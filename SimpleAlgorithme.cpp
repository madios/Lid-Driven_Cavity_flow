//
// Created by Peter Berg Ammundsen on 06/08/2025.
//

#include "SimpleAlgorithme.h"
#include <sstream>
#include <iomanip>
#include <algorithm> // for std::max
#include "include/matrices_print.h"
#include <omp.h>
#include <chrono>

int global_threads = 6;
SimpleAlgorithme::SimpleAlgorithme(parameters* params)
{
    _params = params;
}

void SimpleAlgorithme::initializeParam() const
{
    //omp_set_num_threads(2);
    auto& p = *_params;

    // calculate grid spacing
    p.dx = (p.xmax - p.xmin) / p.nx;
    p.dy = (p.ymax - p.ymin) / p.ny;

    /// set governing parameters
    /// calculating dynamic viscosity mu = rho * Ulid * L / Re
    p.mu           = (p.rho * p.lid_velocity * (p.xmax - p.xmin)) / p.Re;


}

void iniMatrix(std::vector<std::vector<double>>& matrix,const size_t x,const size_t y,const double value)
{
    /// initialize matrix to value.
    /// Size is x and y
       matrix.assign(x,std::vector<double>(y,value));
}

void set_row_to_a_value_in_matrix(std::vector<std::vector<double>>& mat,
                                  size_t row,
                                  double value) {
    if (mat.empty()) return;

    if (row >= mat.size()) {
        throw std::out_of_range("Row index out of bounds");
    }

    size_t cols = mat[row].size();
    for (size_t j = 0; j < cols; ++j) {
        mat[row][j] = value;
    }
}
void set_column_to_a_value_in_matrix(std::vector<std::vector<double>>& mat,
                                     size_t col,
                                     double value) {
    if (mat.empty()) return;

    size_t cols = mat[0].size();
    if (col >= cols) {
        throw std::out_of_range("Column index out of bounds");
    }

    for (size_t i = 0; i < mat.size(); ++i) {
        mat[i][col] = value;
    }
}

void SimpleAlgorithme::Initialize_FVMSolver() const
{
    auto& p = *_params;

    // Initialisering
    iniMatrix(p.u,        p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.v,        p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.p,        p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.pprime,   p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.u_intp,   p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.v_intp,   p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.p_intp,   p.ny + 2, p.nx + 2, 0.0);

    iniMatrix(p.ap_u,     p.ny + 2, p.nx + 2, 1.0);
    iniMatrix(p.ap_v,     p.ny + 2, p.nx + 2, 1.0);
    iniMatrix(p.ap_p,     p.ny + 2, p.nx + 2, 1.0);

    iniMatrix(p.ae,       p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.aw,       p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.an,       p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.as,       p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.pprime_old, p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.upw,      p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.cds,      p.ny + 2, p.nx + 2, 0.0);
    iniMatrix(p.b_mass_imb, p.ny + 2, p.nx + 2, 0.0);

    // copy u and v to u_old and v_old
    p.u_old = p.u;
    p.v_old = p.v;

    // set bc - north boundary
    set_row_to_a_value_in_matrix(p.u,    p.ny + 1, p.lid_velocity);
    p.u[p.ny+1] [p.nx+1] = 0.0;
    // set bc - west boundary
    set_column_to_a_value_in_matrix(p.u,1,0 );

    //set bc - south boundary

    //u(0,:)    = double(0.0);  % MatLab style +1
    set_column_to_a_value_in_matrix(p.v,2,0 );
    //v( 1,:)    = double(0.0);  % MatLab style +1

    // set bc - east boundary
    //%u(:,nx+1) = double(0.0);
    set_column_to_a_value_in_matrix(p.u,p.nx+1,0 );
    //%v(:,nx+1) = double(0.0);

    // set bc - west boundary
    //u(:,2)    = double(0.0);
    //%v(:,0)    = double(0.0);

    p.u_old = p.u;
    p.v_old = p.v;

}


void SimpleAlgorithme::U_MOMENTUM_SOLVER() const
{
    omp_set_num_threads(global_threads);
    using std::max;
    auto& p = *_params;
    p.u_residual = 0.0;
    p.global_res = 0.0;


    // set coefficients for u-momentum discrete equation
    const double base_de = p.mu * p.dy / p.dx;
    const double base_dn = p.mu * p.dx / p.dy;

#pragma omp parallel for collapse(2) schedule(static)
    for (int y = 1; y <= p.ny; ++y) {
        for (int x = 2; x <= p.nx; ++x)
        {
            double de = base_de, dw = base_de;
            double dn = base_dn, ds = base_dn;

            //masseflow (convective mass flux. avg u velocity)
            const double fe = p.rho * p.dy * 0.50 * (p.u_old[y][x + 1] + p.u_old[y][x]);
            const double fw = p.rho * p.dy * 0.50 * (p.u_old[y][x - 1] + p.u_old[y][x]);
            const double fn = p.rho * p.dx * 0.50 * (p.v_old[y + 1][x] + p.v_old[y + 1][x - 1]);
            const double fs = p.rho * p.dx * 0.50 * (p.v_old[y][x]     + p.v_old[y][x - 1]);

            // Boundary conditions
            if (y == p.ny) dn *= 2.0;
            if (y == 1)    ds *= 2.0;

            //CDS  box under (i) page 145
            p.aw[y][x] = dw + fw/2;
            p.ae[y][x] = de - fe/2;
            p.as[y][x] = ds + fs/2;
            p.an[y][x] = dn - fn/2;

            p.ap_u[y][x] = p.aw[y][x] + p.ae[y][x]  + p.as[y][x] + p.an[y][x]
                         + fe - fw + fn - fs;
        }
    }

    //ap_u = ap_u / alpha_u;
    auto inv_alpha = 1.0/(p.alpha_u);
#pragma omp parallel for schedule(static)
    for (size_t y = 0; y < p.ap_u.size(); ++y) {
        for (size_t x = 0; x < p.ap_u[y].size(); ++x) {
            p.ap_u[y][x] *= inv_alpha;
        }
    }

    // iterating u-momentum equation
    // Step 1 in Book
    for (int k = 1; k <= p.u_iter_max; ++k)
        {
#pragma omp parallel for schedule(static)
        for (int y = 1; y <= p.ny; ++y)
        {
            for (int x = 2; x <= p.nx; ++x   )
            {
                // eq 6.36 page 189
                 auto sumAnbUnb =
                          p.ae[y][x] * p.u[y][x+1] +
                          p.aw[y][x] * p.u[y][x-1] +
                          p.an[y][x] * p.u[y+1][x] +
                          p.as[y][x] * p.u[y-1][x];
                auto A = p.dy ;
                p.u[y][x] = (1.0 / p.ap_u[y][x]) *
                            (sumAnbUnb +
                             (p.p[y][x-1] - p.p[y][x])*A
                             +(1.0 - p.alpha_u) * p.ap_u[y][x] * p.u_old[y][x]) ;
            }
        }
    }

    // compute u-momentum residual
    p.u_residual = 0.0;
#pragma omp parallel for schedule(static)
    for (int y = 1; y < p.ny; ++y   ) {
        for (int x = 2; x < p.nx; ++x   )
        {
            double diff = p.u[y][x] - p.u_old[y][x];
            p.u_residual += diff * diff;
        }
    }

    p.u_residual = sqrt(p.u_residual);
    p.global_res = p.u_residual;

}

void SimpleAlgorithme::V_MOMENTUM_SOLVER() const
{
    omp_set_num_threads(global_threads);
    using std::max;
    auto& p = *_params;
    double v_residual = 0.0;
    double global_res = 0.0;

    const double base_de = p.mu * p.dy / p.dx;
    const double base_dn = p.mu * p.dx / p.dy;

    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 2; j <= p.ny; ++j) {
        for (int i = 1; i <= p.nx; ++i) {

            double de = base_de, dw = base_de;
            double dn = base_dn, ds = base_dn;

            //masseflow (convective mass flux. avg v velocity)
            const double fe = p.rho * p.dy * 0.5 * (p.u_old[j][i + 1] + p.u_old[j - 1][i + 1]);
            const double fw = p.rho * p.dy * 0.5 * (p.u_old[j - 1][i] + p.u_old[j][i]);
            const double fn = p.rho * p.dx * 0.5 * (p.v_old[j + 1][i] + p.v_old[j][i]);
            const double fs = p.rho * p.dx * 0.5 * (p.v_old[j - 1][i] + p.v_old[j][i]);

            // Boundary conditions
            if (i == p.ny) de *= 2.0;
            if (i == 1)    dw *= 2.0;

            //CDS  box under (i) page 145
            p.ae[j][i] = de - fe/2;
            p.aw[j][i] = dw + fw/2;
            p.an[j][i] = dn - fn/2;
            p.as[j][i] = ds + fs/2;

            p.ap_v[j][i] = p.ae[j][i] + p.aw[j][i] + p.an[j][i] + p.as[j][i]
                          + fe - fw + fn - fs;

        }
    }

    // ap_v = ap_v / alpha_v
    auto inv_alpha = 1.0/(p.alpha_v);
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < p.ap_v.size(); ++j)
    {
        for (int i = 0; i < p.ap_v[j].size(); ++i)
        {
            p.ap_v[j][i] /= p.alpha_v;
        }
    }

// --- Iteration over v-momentum ligningen ---
// Step 1 in Book
#pragma omp parallel for schedule(static)
for (int k = 1; k < p.v_iter_max; ++k) {
    for (int j = 2; j <= p.ny; ++j) {
        for (int i = 1; i <= p.nx; ++i) {
            // eq 6.37 page 189
            auto sumAnbVnb =p.ae[j][i] * p.v[j][i + 1] +
                          p.aw[j][i] * p.v[j][i - 1] +
                          p.an[j][i] * p.v[j + 1][i] +
                          p.as[j][i] * p.v[j - 1][i] ;
            auto A = p.dy;
            p.v[j][i] = (1.0 / p.ap_v[j][i]) *
                        (sumAnbVnb +
                        (p.p[j - 1][i] - p.p[j][i]) * A)
                      + (1.0 - p.alpha_v) * p.v_old[j][i];
        }
    }
}

// --- Udregn residual ---
v_residual = 0.0;
#pragma omp parallel for schedule(static)
for (int j = 2; j < p.ny; ++j) {
    for (int i = 1; i < p.nx; ++i) {
        double diff = p.v[i][j] - p.v_old[i][j];
        v_residual += diff * diff;
    }
}
v_residual = std::sqrt(v_residual);

global_res += v_residual;


}

double SimpleAlgorithme::PRESSURE_CORRECTION_SOLVER() const
{
    omp_set_num_threads(global_threads);
    using std::max;
    auto& p = *_params;
    double de = 0.0, dw = 0.0, dn = 0.0, ds = 0.0;
    double me = 0.0, mw = 0.0, mn = 0.0, ms = 0.0;
    double v_residual = 0.0;

    // --- set coefficients for pressure correction discrete equation ---
#pragma omp parallel for schedule(static)
    for (int j = 1; j <= p.ny; ++j)
    {
        for (int i = 1; i <= p.nx; ++i)
        {
            // under eq 6.32 page 188
            auto Ay     =  p.dy;
            auto Ax     =  p.dx;
            auto duiJW   = Ay / p.ap_u[j][i];
            auto duip1JE = Ay / p.ap_u[j][i + 1];
            auto dviJp1N = Ax / p.ap_v[j + 1][i];
            auto dviJS   = Ax / p.ap_v[j][i];
            p.ae[j][i] = p.rho * Ay * duip1JE;
            p.aw[j][i] = p.rho * Ay * duiJW;
            p.an[j][i] = p.rho * Ax * dviJp1N;
            p.as[j][i] = p.rho * Ax * dviJS;
        }
    }


    // --- update coefficients along boundaries ---
#pragma omp parallel for schedule(static)
    for (int j = 0; j < p.ny + 2; ++j) p.ae[j][p.nx] = 0.0;
#pragma omp parallel for schedule(static)
    for (int j = 0; j < p.ny + 2; ++j) p.aw[j][1]    = 0.0;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < p.nx + 2; ++i) p.an[p.ny][i] = 0.0;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < p.nx + 2; ++i) p.as[1][i]    = 0.0;
    // ap_p = ae + aw + an + as
#pragma omp parallel for schedule(static)
    for (int j = 0; j < p.ny + 2; ++j)
    {
        for (int i = 0; i < p.nx + 2; ++i)
        {
            // eq 6.32 page 188
            p.ap_p[j][i] = p.ae[j][i] + p.aw[j][i] + p.an[j][i] + p.as[j][i];
        }
    }

    // Reference cell
    p.ap_p[1][1] = 1.0;
    // pprime = 0, pprime_old = pprime
#pragma omp parallel for schedule(static)
    for (int j = 0; j < p.ny + 2; ++j)
    {
        for (int i = 0; i < p.nx + 2; ++i)
        {
            p.pprime[j][i] = 0.0;
        }
    }
    p.pprime_old = p.pprime;

    // --- mass imbalance Init ---
#pragma omp parallel for schedule(static)
    for (int j = 0; j < p.ny; ++j)
    {
        for (int i = 0; i < p.nx; ++i)
        {
            p.b_mass_imb[j][i] = 0.0;
        }
    }
#pragma omp parallel for schedule(static)
    for (int j = 0; j <= p.ny; ++j)
    {
        for (int i = 0; i <= p.nx; ++i)
        {
            // under eq 6.32 page 188
            p.b_mass_imb[j][i] =
                p.rho * p.dy * (p.u[j][i] - p.u[j][i+1]) +
                p.rho * p.dx * (p.v[j][i] - p.v[j+1][i]);
        }
    }

    double mass_error = 0.0;
#pragma omp parallel for schedule(static)
    for (int j = 1; j < p.ny; ++j)
    {
        for (int i = 1; i < p.nx; ++i)
        {
            double t = p.b_mass_imb[j][i];
            mass_error += t * t;
        }
    }
    mass_error = std::sqrt(mass_error);

    // --- solve pressure correction (Jacobi, parallel over j,i) ---
    // Step 2 in Book
    for (int k = 1; k <= p.pressure_iter_max; ++k)
    {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j <= p.ny; ++j)
        {
            for (int i = 1; i <= p.nx; ++i)
            {
                //Step 2 page 190
                const double inv_ap = 1.0 / p.ap_p[j][i];
                p.pprime[j][i] = inv_ap * (
                      p.ae[j][i] * p.pprime_old[j][i + 1]
                    + p.aw[j][i] * p.pprime_old[j][i - 1]
                    + p.an[j][i] * p.pprime_old[j + 1][i]
                    + p.as[j][i] * p.pprime_old[j - 1][i]
                    + p.b_mass_imb[j][i]
                );
            }
        }
        // swap buffers to avoid deep copy; next iter reads from p.pprime_old
        p.pprime.swap(p.pprime_old);
    }
    // ensure p.pprime holds latest after an odd number of swaps
    if (p.pressure_iter_max % 2 != 0) {
        p.pprime.swap(p.pprime_old);
    }

    // apply pressure corrections to pressure (under-relaxed)
    // Step 3
#pragma omp parallel for schedule(static)
    for (int j = 1; j <= p.ny; ++j)
    {
        for (int i = 1; i <= p.nx; ++i)
        {
            p.p[j][i] += p.alpha_p *p.pprime[j][i]; // pnew = p* + alpha_p * p' (eq 6.33)
        }
    }

    // apply velocity corrections to u-velocity
#pragma omp parallel for schedule(static)
    for (int j = 1; j <= p.ny; ++j)
    {
        for (int i = 2; i <= p.nx; ++i)
        {
            // u(iJ) = u*(i,J) + d(i,J) *(p'(I-1,J) - p'(I,J))   Step 3 page 190
            auto d = p.dy / p.ap_u[j][i];
            p.u[j][i] = p.u[j][i] + d * (p.pprime[j][i - 1] - p.pprime[j][i]) ;
        }
    }
    // apply velocity corrections to v-velocity
#pragma omp parallel for schedule(static)
    for (int j = 2; j <= p.ny; ++j)
    {
        for (int i = 0; i <= p.nx; ++i)
        {
            // v(iJ) = u*(i,J) + d(i,J) *(p'(I-1,J) - p'(I,J))   Step 3 page 190
            p.v[j][i] += (p.dx / p.ap_v[j][i]) * (p.pprime[j - 1][i] - p.pprime[j][i]);
        }
    }
    // computing mass imbalance after pressure correction
#pragma omp parallel for schedule(static)
    for (int j = 0; j < p.ny + 2; ++j)
    {
        for (int i = 0; i < p.nx + 2; ++i)
        {
            p.b_mass_imb[j][i] = 0.0;
        }
    }
#pragma omp parallel for schedule(static)
    for (int j = 1; j <= p.ny; ++j)
    {
        for (int i = 1; i <= p.nx; ++i)
        {      
            p.b_mass_imb[j][i] =
                p.rho * p.dy * (p.u[j][i + 1] - p.u[j][i]) +
                p.rho * p.dx * (p.v[j + 1][i] - p.v[j][i]);
        }
    }
    mass_error = 0.0;
    #pragma omp parallel for collapse(2) reduction(+:mass_error) schedule(static)
    for (int j = 0; j < p.ny + 2; ++j)
    {
        for (int i = 0; i < p.nx + 2; ++i)
        {
            const double t = p.b_mass_imb[j][i];
            mass_error += t * t;
        }
    }
    mass_error = std::sqrt(mass_error);

    // global_res opdateres
    p.global_res += mass_error;

    // update velocity field
    p.v_old = p.v;
    p.u_old = p.u;

    // No Step 4, np temperatur

    return mass_error;
}

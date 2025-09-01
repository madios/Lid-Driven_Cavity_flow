//
// Created by Peter Berg Ammundsen on 08/08/2025.
//

#ifndef MATRICES_PRINT_H
#define MATRICES_PRINT_H


// matrices_print.hpp
#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include "../parameters.h"
#include <QPlainTextEdit>

inline bool g_print = false;

// Global pointer til console
extern QPlainTextEdit* g_consoleOutput;

namespace matprint {

    using Matrix = std::vector<std::vector<double>>;
    inline void printToGlobalConsole(const QString& text) {
        if (g_consoleOutput) {
            g_consoleOutput->appendPlainText(text);
            g_consoleOutput->repaint();              // tvinger immediate repaint af widget
            QCoreApplication::processEvents(QEventLoop::AllEvents);
        }
    }
    struct DiffPos {
        size_t row;
        size_t col;
        double cppValue;
        double matlabValue;
        double diff;
    };

    inline std::string matrixToString(
        const Matrix& mat,
        const std::string& name = "Matrix",
        const std::vector<DiffPos>& diffPos = {},
        int precision = 6
        )
    {
        const size_t rows = mat.size();
        const size_t cols = rows ? mat[0].size() : 0;

        // Byg en maske med markører for hurtige O(1) opslag
        std::vector<std::vector<char>> mark(rows, std::vector<char>(cols, ' '));
        for (const auto& d : diffPos) {
            if (d.row < rows && d.col < cols) {
                mark[d.row][d.col] = '*';
            }
        }

        // bredder: 'width' til tallet, +1 til markør
        const int width = precision + 3;    // fx 3 decimaler -> ca. 6-7 tegn, + lidt margin
        const int colWidth = width + 1;     // ekstra plads til markør

        std::ostringstream oss;
        oss << name << " (" << rows << "x" << cols << "):\n\n";

        // Header med kolonneindeks
        oss << std::setw(8) << " ";
        for (size_t j = 0; j < cols; ++j)
            oss << std::setw(colWidth) << j;
        oss << "\n";

        // Rækker
        oss << std::fixed << std::setprecision(precision);
        for (size_t i = 0; i < rows; ++i) {
            oss << std::setw(6) << i << " |";
            for (size_t j = 0; j < cols; ++j) {
                oss << std::setw(width) << mat[i][j]
                    << mark[i][j]; // ' ' eller '*'
            }
            oss << "\n";
        }

        return oss.str();
    }

    inline std::string matricesToString(const std::vector<std::pair<std::string, Matrix>>& mats)
    {
        std::ostringstream oss;
        bool first = true;
        for (const auto& [name, m] : mats) {
            if (!first) oss << "\n";
            first = false;
            oss << matrixToString(m, name);
        }
        return oss.str();
    }



inline void updateAllMatrix(parameters& p, int width = 10, int precision = 3)
{
    using matprint::matrixToString;

    p.p_s           = matrixToString(p.p,           "p"       );
    p.u_s           = matrixToString(p.u,           "u"       );
    p.v_s           = matrixToString(p.v,           "v"       );
    p.pprime_s      = matrixToString(p.pprime,      "pprime"   );
    p.pprime_old_s  = matrixToString(p.pprime_old,  "pprime_ld");

    p.u_intp_s      = matrixToString(p.u_intp,      "u_intp"   );
    p.v_intp_s      = matrixToString(p.v_intp,      "v_intp"   );
    p.p_intp_s      = matrixToString(p.p_intp,      "p_intp"   );

    p.ap_u_s        = matrixToString(p.ap_u,        "ap_u"    );
    p.ap_v_s        = matrixToString(p.ap_v,        "ap_v"    );
    p.ap_p_s        = matrixToString(p.ap_p,        "ap_p"    );

    p.ae_s          = matrixToString(p.ae,          "ae"      );
    p.aw_s          = matrixToString(p.aw,          "aw"      );
    p.as_s          = matrixToString(p.as,          "as"      );
    p.an_s          = matrixToString(p.an,          "an"      );

    p.upw_s         = matrixToString(p.upw,         "upw"     );
    p.cds_s         = matrixToString(p.cds,         "cds"     );

    p.v_old_s       = matrixToString(p.v_old,       "v_old"   );
    p.u_old_s       = matrixToString(p.u_old,       "u_old"   );

    p.mass_imb_s    = matrixToString(p.b_mass_imb,    "mass_im" );
}

// Læser en matrix fra tekstfil
inline std::vector<std::vector<double>> readMatrixFromFile(const std::string& filepath) {
    std::ifstream f(filepath);
    if (!f) throw std::runtime_error("Kunne ikke åbne: " + filepath);
    std::vector<std::vector<double>> M;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream iss(line);
        std::vector<double> row; double v;
        while (iss >> v) row.push_back(v);
        if (!row.empty()) M.push_back(std::move(row));
    }
    return M;
}










    inline std::vector<DiffPos> compareMatricesDetailed(
        const std::vector<std::vector<double>>& a,
        const std::vector<std::vector<double>>& b,
        const std::string& name,
        double tol = 1e-8)
    {
        std::vector<DiffPos> differences;
        if (a.size() != b.size()) {
            std::cerr << "Fejl i " << name << ": rækker " << a.size() << " vs " << b.size() << "\n";
            return differences;
        }

        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i].size() != b[i].size()) {
                std::cerr << "Fejl i " << name << ": række " << i
                          << " kolonner " << a[i].size() << " vs " << b[i].size() << "\n";
                return differences;
            }
            for (size_t j = 0; j < a[i].size(); ++j) {
                double diff = std::fabs(a[i][j] - b[i][j]);
                if (diff > tol) {
                   //// std::cerr << "Fejl i " << name << " ved (" << i << "," << j
                   //           << "): C++=" << a[i][j]
                   //           << " MATLAB=" << b[i][j]
                   //           << " diff=" << diff << "\n";
                    differences.push_back({i, j, a[i][j], b[i][j], diff});
                }
            }
        }

        return differences;
    }

// Tjekker alle matricer i p mod filer i TestFiler-mappen
inline void checkAllMatrices(const parameters& p, std::string name = "", double tol = 1e-9)
    {
    const std::string folder =
        "/Users/ri03jm/Library/Mobile Documents/com~apple~CloudDocs/"
        "DropboxFolder/Peter/skole/PhD/MatlabScript/CFD/LidcavityFlow/TestFiler";

    if (!std::filesystem::exists(folder)) {
        std::cerr << "Mappe findes ikke: " << folder << "\n";
        return;
    }
    auto check = [&](const auto& mat, const std::string& file) {
        try {
            auto fileMat = readMatrixFromFile(folder + "/" + file);
            auto diffPosition = compareMatricesDetailed(mat, fileMat, file, tol);
            if (diffPosition.size() != 0)
            {
                 //matprint::printToGlobalConsole(matrixToString(mat) + " /n " + matrixToString(fileMat));
                 std::string mat_s = file + "MatA" + matrixToString(mat,"MatA",diffPosition) + " /n MatB" + matrixToString(fileMat);
                 std::cerr  << mat_s;
                 //QString mat_q = QString::fromStdString(mat_s);
                 //matprint::printToGlobalConsole(mat_q);


            }

        } catch (const std::exception& e) {
            std::cerr << "Fejl ved " << file << ": " << e.what() << "\n";
        }
    };
        std::cerr  << name ;
        if (!g_print)
        {
            std::cerr  << " - No print";
            std::cerr  << std::endl; 
            return;
        }
        std::cerr  << std::endl;
        check(p.u,          "u.txt");
        check(p.u_old,          "u_old.txt");
        check(p.v_old,          "v_old.txt");
        check(p.p,          "p.txt");
        check(p.pprime,     "pprime.txt");
        check(p.u_intp,     "u_intp.txt");
        check(p.v_intp,     "v_intp.txt");
        check(p.p_intp,     "p_intp.txt");
        check(p.ap_u,       "ap_u.txt");
        check(p.ap_v,       "ap_v.txt");
        check(p.ap_p,       "ap_p.txt");
        check(p.ae,         "ae.txt");
        check(p.aw,         "aw.txt");
        check(p.as,         "as.txt");
        check(p.an,         "an.txt");
        check(p.upw,        "upw.txt");
        check(p.cds,        "cds.txt");
        check(p.b_mass_imb,   "mass_imb.txt");
}







} // namespace matprint

#endif //MATRICES_PRINT_H

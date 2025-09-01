//
// Created by Peter Berg Ammundsen on 05/08/2025.
//

#include "LidDrivenCavityFlow.h"
#include <vector>
#include <cmath>
#include <QThread>
#include "SimpleAlgorithme.h"
#include "parameters.h"
#include "include/matrices_print.h"



std::vector<double> LidDrivenCavityFlow::simulate(QPlainTextEdit *consoleOutput, std::function<void(int, double,std::vector<double>)> updateGraph)
{
    parameters params;
    auto& p = params;
    params.load_from_json("input_parameters.json");

    consoleOutput->appendPlainText(("lid_velocity: " + std::to_string(params.lid_velocity) + " m/s").c_str());
    consoleOutput->appendPlainText(("Reynolds number: " + std::to_string(params.Re)).c_str());
    consoleOutput->appendPlainText(("number of cell in x directions: " + std::to_string(params.nx)).c_str());
    consoleOutput->appendPlainText(("number of cell in y directions: " + std::to_string(params.ny)).c_str());
    matprint::updateAllMatrix(params, /*width=*/8, /*precision=*/4);

    SimpleAlgorithme Simple(&params);
    Simple.initializeParam();
    Simple.Initialize_FVMSolver();



    std::vector<double> values;
    for (int i = 0; i < 3000; i++) {
        double value = 0;
        double value1 = i;
        for (int j = 0; j < 40; j++) {
            Simple.U_MOMENTUM_SOLVER();
            Simple.V_MOMENTUM_SOLVER();
            value = Simple.PRESSURE_CORRECTION_SOLVER();
        }
        values.push_back(value);
        if (consoleOutput) {
            consoleOutput->appendPlainText(QString("Value %1: %2").arg(i).arg(value));
            if (updateGraph)
            {
                //std::vector uLine = p.u[35];
                size_t rows   = p.u.size();
                size_t cols   = rows ? p.u[0].size() : 0;    // antag rektangulær matrix
                size_t midCol = cols / 2;                    // eller: size_t midCol = (p.nx + 2) / 2;

                std::vector<double> uCol;
                uCol.reserve(rows);
                for (size_t r = 0; r < rows; ++r) {
                    uCol.push_back(p.u[r][midCol]);
                }

                updateGraph(i, value,uCol);
                QCoreApplication::processEvents();
            }
            QCoreApplication::processEvents();
        }
        QThread::msleep(1);
    }
    return values;
}



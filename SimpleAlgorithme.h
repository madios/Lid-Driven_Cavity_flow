//
// Created by Peter Berg Ammundsen on 06/08/2025.
//

#ifndef SIMPLEALGORITHME_H
#define SIMPLEALGORITHME_H
#include "parameters.h"

class SimpleAlgorithme {
public:
    explicit SimpleAlgorithme(parameters* params);
    void initializeParam() const;
    void Initialize_FVMSolver() const;
    void U_MOMENTUM_SOLVER() const;
    void V_MOMENTUM_SOLVER() const;
    double PRESSURE_CORRECTION_SOLVER() const;
private:
    parameters* _params = nullptr;
};



#endif //SIMPLEALGORITHME_H

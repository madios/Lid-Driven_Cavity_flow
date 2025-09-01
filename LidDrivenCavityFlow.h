//
// Created by Peter Berg Ammundsen on 05/08/2025.
//

#include <QPlainTextEdit>

#ifndef LIDDRIVENCAVITYFLOW_H
#define LIDDRIVENCAVITYFLOW_H



class LidDrivenCavityFlow {
public:
    std::vector<double> simulate(QPlainTextEdit *consoleOutput, std::function<void(int, double,std::vector<double>)> updateGraph = nullptr);
};



#endif //LIDDRIVENCAVITYFLOW_H

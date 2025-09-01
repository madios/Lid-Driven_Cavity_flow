#include <QMainWindow>
#include <QPushButton>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPainter>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include "LidDrivenCavityFlow.h"
#include <QtCharts/QChartView>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QLogValueAxis>
#include <algorithm>
#include <vector>
#include <limits>
#include <QCoreApplication>
#include <QApplication>
#include "include/matrices_print.h"
#include "Benchmark.h"

// Global console pointer so other translation units can write to the app console
QPlainTextEdit* g_consoleOutput = nullptr;
//QPlainTextEdit *consoleOutput = new QPlainTextEdit(centralWidget);
//consoleOutput->setReadOnly(true);

// Sæt global pointer
//g_consoleOutput = consoleOutput;



void printToConsole(QPlainTextEdit *console, const QString &text) {
    console->appendPlainText(text);
}

std::vector<double> linspace(double start, double stop, std::size_t n)
{
    std::vector<double> v;
    if (n == 0) return v;
    if (n == 1) { v.push_back(stop); return v; } // matcher MATLAB: linspace(a,b,1) = [b]

    v.resize(n);
    const double step = (stop - start) / static_cast<double>(n - 1);
    for (std::size_t i = 0; i < n; ++i)
        v[i] = start + step * static_cast<double>(i);
    v.back() = stop; // sikrer præcis slutværdi
    return v;
}

int main(int argc, char *argv[]) {
    QApplication a(argc, argv);
    QMainWindow mainWindow;
    mainWindow.resize(2500, 768);
    mainWindow.setWindowTitle("Simple MainWindow");

    QWidget *centralWidget = new QWidget(&mainWindow);
    mainWindow.setCentralWidget(centralWidget);

    QPlainTextEdit *consoleOutput = new QPlainTextEdit(centralWidget);
    consoleOutput->setReadOnly(true);
    QFont monoFont("Courier New"); // eller "Monospace" / "Consolas"
    monoFont.setStyleHint(QFont::Monospace);
    consoleOutput->setFont(monoFont);
    // expose console globally
    g_consoleOutput = consoleOutput;
    //consoleOutput->setGeometry(50, 30, 300, 50);
    consoleOutput->appendPlainText("Hello World!");
    printToConsole(consoleOutput, "Velkommen til programmet!");

    LidDrivenCavityFlow simulation;

    QPushButton *button = new QPushButton("Click Me", centralWidget);

    // === Figure 1 (øverst) ===
    QLineSeries *series1 = new QLineSeries();
    QChart *chart1 = new QChart();
    chart1->addSeries(series1);

    QValueAxis *axisX1 = new QValueAxis();
    QValueAxis *axisY1 = new QValueAxis();
    axisX1->setTitleText("Step");
    axisY1->setTitleText("Value");
    chart1->addAxis(axisX1, Qt::AlignBottom);
    chart1->addAxis(axisY1, Qt::AlignLeft);
    series1->attachAxis(axisX1);
    series1->attachAxis(axisY1);
    chart1->setTitle("Figure 1");

    // Benchmark-curve på Figure 1 (Re1000)
    Benchmark bm;
    QLineSeries *series1_benchmark = new QLineSeries();
    auto rownumber = 3;
    series1_benchmark->setName(QString("Benchmark %1").arg(bm.headers[rownumber])); // 3 = Re1000
    chart1->addSeries(series1_benchmark);
    series1_benchmark->attachAxis(axisX1);
    series1_benchmark->attachAxis(axisY1);
    for (const auto &row : bm.data)
    {
        double y = row[0];
        double u = row[rownumber];
        series1_benchmark->append(u, y);
    }
    axisX1->setTickCount(8);
    axisY1->setTickCount(11);
    axisX1->setRange(-0.4, 1);
    axisY1->setRange(0, 1);
    series1_benchmark->setPointsVisible(true);

    QChartView *chartView1 = new QChartView(chart1, centralWidget);
    chartView1->setRenderHint(QPainter::Antialiasing);

    // === Figure 2 (nederst) ===
    QLineSeries *series2 = new QLineSeries();
    QChart *chart2 = new QChart();
    chart2->addSeries(series2);

    QValueAxis *axisX2 = new QValueAxis();
    QLogValueAxis *axisY2 = new QLogValueAxis();
    axisX2->setTitleText("Step");
    axisY2->setBase(10);
    axisY2->setLabelFormat("%g");
    axisY2->setTitleText("Value (log)");
    chart2->addAxis(axisX2, Qt::AlignBottom);
    chart2->addAxis(axisY2, Qt::AlignLeft);
    series2->attachAxis(axisX2);
    series2->attachAxis(axisY2);
    chart2->setTitle("Figure 2");

    QChartView *chartView2 = new QChartView(chart2, centralWidget);
    chartView2->setRenderHint(QPainter::Antialiasing);

    // Container to accumulate all values during simulation
    std::vector<double> allValues;
    std::vector<double> allValues1;

    QHBoxLayout *mainLayout = new QHBoxLayout(centralWidget);

    // Left layout (console + button)
    QVBoxLayout *leftLayout = new QVBoxLayout();
    leftLayout->addWidget(consoleOutput);
    leftLayout->addWidget(button);
    leftLayout->setStretch(0, 10);
    leftLayout->setStretch(1, 0);

    // Right layout (to grafer under hinanden)
    QVBoxLayout *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(chartView1);
    rightLayout->addWidget(chartView2);
    rightLayout->setStretch(0, 1);
    rightLayout->setStretch(1, 1);

    // Add to main layout
    mainLayout->addLayout(leftLayout, 1);
    mainLayout->addLayout(rightLayout, 1);

    QObject::connect(button, &QPushButton::clicked, [&]()
    {
        // Create two extra buttons on first click
        static QPushButton* pauseButton = nullptr;
        static QPushButton* stopButton  = nullptr;
        if (!pauseButton) {
            pauseButton = new QPushButton("Figure 1", centralWidget);
            leftLayout->addWidget(pauseButton);
            QObject::connect(pauseButton, &QPushButton::clicked, [&](){
                consoleOutput->appendPlainText("[UI] Pause clicked");
            });
        }
        if (!stopButton) {
            stopButton = new QPushButton("Figure 2", centralWidget);
            leftLayout->addWidget(stopButton);
            QObject::connect(stopButton, &QPushButton::clicked, [&](){
                consoleOutput->appendPlainText("[UI] Stop clicked");
            });
        }

        allValues.clear();
        allValues1.clear();
        const int WINDOW = 1000; // sidste 1000 punkter i Figure 1

        auto values = simulation.simulate(consoleOutput, [&](int i, double value,std::vector<double> value1)
        {
            // Saml alle værdier
            allValues.push_back(value);
            //allValues1.push_back(value1);
            allValues1.clear();
            allValues1 = value1;
            auto y_vec = linspace(0.0, 1.0, allValues1.size());
            // === Opdatér Figure 1 (sidste 200) ===
            {
                series1->clear();
                for (int i = 0; i < static_cast<int>(allValues1.size()); ++i) {
                    double x = allValues1[i];
                    double y = y_vec[i];
                    series1->append(x, y);
                }
                axisX1->setRange(-0.4, 1);
                axisY1->setRange(0, 1);
                axisX1->setTickCount(8);
                axisY1->setTickCount(11);
                series1->setName(QString("model %1").arg(bm.headers[3])); // 3 = Re1000
            }

            // === Opdatér Figure 2 (hele historikken) ===
            {
                series2->clear();
                constexpr double EPS = 1e-6; // log-akse kræver > 0
                double yMin2 = std::numeric_limits<double>::max();
                double yMax2 = EPS;
                for (int x = 0; x < static_cast<int>(allValues.size()); ++x) {
                    double y = allValues[x];
                    double safeY = (y <= 0.0) ? EPS : y; // forskyd ikke-positive værdier
                    series2->append(x, safeY);
                    if (safeY < yMin2) yMin2 = safeY;
                    if (safeY > yMax2) yMax2 = safeY;
                }
                axisX2->setRange(0, std::max(0, i));
                if (!(yMax2 > yMin2)) { yMin2 = EPS; yMax2 = 10 * EPS; }
                axisY2->setRange(yMin2, yMax2);
            }

            QCoreApplication::processEvents();
        });
        // After simulate returns, leave the last window on screen (nothing else to do)
    });

    mainWindow.show();
    return QApplication::exec();
}
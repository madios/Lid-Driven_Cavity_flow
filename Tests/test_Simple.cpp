//
// Created by Peter Berg Ammundsen on 14/08/2025.
//


#include <gtest/gtest.h>
#include "../SimpleAlgorithme.h"

class SimpleTest : public ::testing::Test {
protected:
    parameters p;
    SimpleAlgorithme* solver;

    void SetUp() override {
        p.xmin = 0;
        p.xmax = 1;
        p.nx = 10;
        solver = new SimpleAlgorithme(&p);
        solver->initializeParam();
    }

    void TearDown() override {
        delete solver;
    }
};

TEST_F(SimpleTest, DxCorrect) {
    EXPECT_DOUBLE_EQ(p.dx, 0.1);
}


TEST(SimpleAlgorithmeTest, GridSpacingCalculation) {
    parameters p;
    p.xmin = 0;
    p.xmax = 1;
    p.nx = 10;

    SimpleAlgorithme simple(&p);
    simple.initializeParam();

    EXPECT_DOUBLE_EQ(p.dx, 0.1);
}
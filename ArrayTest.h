#ifndef ARRAYTEST_H
#define ARRAYTEST_H

#include "gtest/gtest.h"
#include "Array.h"

class ArrayTest : public ::testing::Test {
protected:
    bool contains(const Array &matrix, const double value, const double epsilon = 1e-10);
    void expectEqual(const Array& expected, const Array& actual, const std::string& message = "");
};
#endif
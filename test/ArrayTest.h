#ifndef ARRAYTEST_H
#define ARRAYTEST_H

#include "gtest/gtest.h"
#include "../src/Array.h"

class ArrayTest : public ::testing::Test {
protected:
    bool contains(const Array &matrix, const double value, const double epsilon = 1e-10);
    bool isEqual(const Array &expected, const Array &actual);
    void expectEqual(const Array &expected, const Array &actual, const std::string &message = "");
    static constexpr double EPSILON = 1e-14;
};
#endif
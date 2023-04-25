#ifndef ARRAYTEST_H
#define ARRAYTEST_H

#include "gtest/gtest.h"
#include "../src/Array.h"

class ArrayTest : public ::testing::Test {
protected:
    static constexpr double EPSILON = 1e-14;
    bool contains(const Array &matrix, const double value, const double epsilon = EPSILON);
    bool isEqual(const Array &expected, const Array &actual, const double epsilon = EPSILON);
    void expectEqual(const Array &expected, const Array &actual, const std::string &message = "", const double epsilon = EPSILON);
};
#endif
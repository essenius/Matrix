#ifndef ARRAYTEST_H
#define ARRAYTEST_H

#include "gtest/gtest.h"
#include "../src/Array.h"

class ArrayTest : public ::testing::Test {
protected:
    bool contains(const Array &matrix, const double value, const double epsilon = Array::EPSILON);
    bool isEqual(const Array &expected, const Array &actual, const double epsilon = Array::EPSILON);
    void expectEqual(const Array &expected, const Array &actual, const std::string &message = "", const double epsilon = Array::EPSILON);
};
#endif
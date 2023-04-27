#include <gtest/gtest.h>
#include "../src/Matrix.h"
#include "ArrayTest.h"

class MatrixTest : public ArrayTest {
protected:
    void expectNormalizedEqual(const Matrix& expected, const Matrix& actual, const std::string& message = "", const double epsilon = Array::EPSILON);
};
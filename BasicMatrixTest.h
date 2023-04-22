#include <gtest/gtest.h>
#include "BasicMatrix.h"
#include "ArrayTest.h"

class BasicMatrixTest : public ArrayTest {
protected:
    void expectNormalizedEqual(const BasicMatrix& expected, const BasicMatrix& actual, const std::string& message = "");
};
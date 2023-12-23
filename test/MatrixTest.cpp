// Copyright 2023 Rik Essenius
// 
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file
// except in compliance with the License. You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software distributed under the License
// is distributed on an "AS IS" BASIS WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and limitations under the License.

#include <gtest/gtest.h>
#include "MatrixTest.h"

namespace RixMatrixTest {
    using RixMatrix::Matrix;

    void MatrixTest::expectNormalizedEqual(const Matrix& expected, const Matrix& actual, const std::string& message, const double epsilon) {
        expectEqual(expected.normalized(), actual.normalized(), message, epsilon);
    }

    TEST_F(MatrixTest, add) {
        const Matrix m({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
        const Matrix n({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
        const Matrix expected({ {2, 4, 6}, {8, 10, 12}, {14, 16, 18} });
        const Matrix actual = m + n;
        expectEqual(expected, actual, "add");
    }

    TEST_F(MatrixTest, multiplyScalar) {
        Matrix m({ {1, 2}, {3, 4} });
        m *= 2;
        expectEqual(Matrix({ {2, 4}, {6, 8} }), m);
        const Matrix n = 2 * m;
        expectEqual(Matrix({ {4, 8}, {12, 16} }), n);
    }

    TEST_F(MatrixTest, multiply2dMatrix) {
        const Matrix m({ {1, 2}, {3, 4} });
        const Matrix n({ {1, 2}, {3, 4} });
        const Matrix o = m * n;
        expectEqual(Matrix({ {7, 10}, {15, 22} }), o);
    }

    TEST_F(MatrixTest, multiplyDifferentSizeMatrix) {
        const Matrix m({ {1, 2, 3}, {3, -1, 0} });
        const Matrix n({ {1, 2}, {3, 4}, {5, 6} });
        const Matrix o = m * n;
        expectEqual(Matrix({ {22, 28}, {0, 2} }), o, "m * n");
        const Matrix p = n * m;
        expectEqual(Matrix({ {7, 0, 3}, {15, 2, 9}, {23, 4, 15} }), p, "n * m");
    }

    TEST_F(MatrixTest, determinant) {
        const Matrix m({ {1, 2}, {3, 4} });
        EXPECT_EQ(-2, m.getDeterminant());
        const Matrix n({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
        EXPECT_EQ(0, n.getDeterminant());
        const Matrix o({ {1, 3, 5, 9}, {1, 3, 1, 7}, {4, 3, 9, 7}, {5, 2, 0, 9} });
        EXPECT_EQ(-376, o.getDeterminant());
    }

    TEST_F(MatrixTest, cofactor2d) {
        const Matrix m({ {1, 2}, {3, 4} });
        EXPECT_EQ(4, m.getCofactor(0, 0));
        EXPECT_EQ(-3, m.getCofactor(0, 1));
        EXPECT_EQ(-2, m.getCofactor(1, 0));
        EXPECT_EQ(1, m.getCofactor(1, 1));
    }

    TEST_F(MatrixTest, cofactor3d) {
        const Matrix n({ {1, 2, 1}, {6, -1, 0}, {-1, -2, -1} });
        EXPECT_EQ(1, n.getCofactor(0, 0));
        EXPECT_EQ(6, n.getCofactor(0, 1));
        EXPECT_EQ(-13, n.getCofactor(0, 2));
        EXPECT_EQ(0, n.getCofactor(1, 0));
        EXPECT_EQ(0, n.getCofactor(1, 1));
        EXPECT_EQ(0, n.getCofactor(1, 2));
        EXPECT_EQ(1, n.getCofactor(2, 0));
        EXPECT_EQ(6, n.getCofactor(2, 1));
        EXPECT_EQ(-13, n.getCofactor(2, 2));
    }

    TEST_F(MatrixTest, adjugate2d) {
        const Matrix m({ {1, 2}, {3, 4} });
        const Matrix expected({ {4, -2}, {-3, 1} });
        const Matrix actual = m.getAdjugate();
        expectEqual(expected, actual);
    }

    TEST_F(MatrixTest, adjugate3d) {
        const Matrix m({ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} });
        const Matrix expected({ {-24, 18, 5}, {20, -15, -4}, {-5, 4, 1} });
        const Matrix actual = m.getAdjugate();
        expectEqual(expected, actual);
    }

    TEST_F(MatrixTest, inverse2d) {
        const Matrix m({ {1, 2}, {3, 4} });
        const Matrix actual = m.inverted();
        expectEqual(Matrix({ {-2, 1}, {1.5, -0.5} }), actual, "inverse");
        expectEqual(m, actual.inverted(), "inverse of inverse");
    }

    TEST_F(MatrixTest, inverse3d) {
        const Matrix m({ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} });
        const Matrix expected({ {-24, 18, 5}, {20, -15, -4}, {-5, 4, 1} });
        const Matrix actual = m.inverted();
        expectEqual(expected, actual, "inverse");
        expectEqual(m, actual.inverted(), "inverse of inverse");
    }

    TEST_F(MatrixTest, normalize) {
        const Matrix m({ {3, 4} });
        const Matrix expected({ {0.6, 0.8} });
        const Matrix actual = m.normalized();
        expectEqual(expected, actual, "normalize all positive");
        expectEqual(Matrix({ {0} }), Matrix({ {0} }).normalized(), "normalize zero");
        const Matrix m2({ {-1, 4, -2, 2} });
        const Matrix expected2({ {0.2, -0.8, 0.4, -0.4} });
        expectEqual(expected2, m2.normalized(), "normalize first negative");
    }

    TEST_F(MatrixTest, toArray) {
	    constexpr std::initializer_list<std::initializer_list<double>> Input = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
        const Matrix m(Input);
        const Array expected(Input);
        const auto actual = m.toArray();
        expectEqual(expected, actual, "toArray");
    }

    TEST_F(MatrixTest, minor3d) {
        const Matrix m({ {3, 5, 0}, {2, -1, -7}, {6, -1, 5} });
        expectEqual(Matrix({ {-1, -7}, {-1, 5} }), m.getMinor(0, 0));
        expectEqual(Matrix({ {2, -7}, {6, 5} }), m.getMinor(0, 1));
        expectEqual(Matrix({ {2, -1}, {6, -1} }), m.getMinor(0, 2));
        expectEqual(Matrix({ {5, 0}, {-1, 5} }), m.getMinor(1, 0));
        expectEqual(Matrix({ {3, 0}, {6, 5} }), m.getMinor(1, 1));
        expectEqual(Matrix({ {3, 5}, {6, -1} }), m.getMinor(1, 2));
        expectEqual(Matrix({ {5, 0}, {-1, -7} }), m.getMinor(2, 0));
        expectEqual(Matrix({ {3, 0}, {2, -7} }), m.getMinor(2, 1));
        expectEqual(Matrix({ {3, 5}, {2, -1} }), m.getMinor(2, 2));
    }

    TEST_F(MatrixTest, cofactor3dbis) {
        const Matrix m({ {3, 5, 0}, {2, -1, -7}, {6, -1, 5} });
        EXPECT_EQ(-12, m.getCofactor(0, 0));
        EXPECT_EQ(-52, m.getCofactor(0, 1));
        EXPECT_EQ(4, m.getCofactor(0, 2));
        EXPECT_EQ(-25, m.getCofactor(1, 0));
        EXPECT_EQ(15, m.getCofactor(1, 1));
        EXPECT_EQ(33, m.getCofactor(1, 2));
        EXPECT_EQ(-35, m.getCofactor(2, 0));
        EXPECT_EQ(21, m.getCofactor(2, 1));
        EXPECT_EQ(-13, m.getCofactor(2, 2));
    }

    TEST_F(MatrixTest, assignArray) {
        const Array a({ {1, 2}, {3, 4} });
        const auto m = Matrix(a);
        expectEqual(a, m);
    }

    TEST_F(MatrixTest, assertTest) {
        const Matrix m({ {1, 2} });
        EXPECT_DEATH(m.getCofactor(2, 2), "Assertion failed");
        EXPECT_DEATH(m.getDeterminant(), "Assertion failed");
        EXPECT_DEATH(m.getMinor(2, 2), "Assertion failed");
        EXPECT_DEATH(m.getMinor(0, 0), "Assertion failed");
        EXPECT_DEATH(m.getTrace(), "Assertion failed");
        EXPECT_DEATH(m.inverted(), "Assertion failed");
        EXPECT_FALSE(m.isInvertible());
    }
}
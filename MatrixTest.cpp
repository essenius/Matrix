#include <gtest/gtest.h>
#include "Matrix.h"

class MatrixTest : public ::testing::Test {
protected:
    void expectEqual(const Matrix& expected, const Matrix& actual) {

        EXPECT_EQ(expected.rows(), actual.rows());
        EXPECT_EQ(expected.columns(), actual.columns());
        for (Dimension i = 0; i < expected.rows(); ++i) {
            for (Dimension j = 0; j < expected.columns(); ++j) {
                constexpr double EPSILON = 0.000000000000001;
                EXPECT_NEAR(expected(i, j), actual(i, j), EPSILON);
            }
        }
    }

    void expectNormalizedEqual(const Matrix& expected, const Matrix& actual)  {
        expectEqual(expected.normalize(), actual.normalize());
    }
};

TEST_F(MatrixTest, initMatrix) {
    Matrix m(2, 2);
    m(0, 0) = 1;
    EXPECT_EQ(1, m(0, 0));
    m(0, 1) = 2;
    EXPECT_EQ(2, m(0, 1));
    m(1, 0) = 3;
}

TEST_F(MatrixTest, copyAddMatrix) {
    Matrix m({ {1, 2}, {3,4} });
    EXPECT_TRUE(m.isSquare());
    Matrix n = m;
    EXPECT_EQ(1, n(0, 0));
    EXPECT_EQ(2, n(0, 1));
    EXPECT_EQ(3, n(1, 0));
    EXPECT_EQ(4, n(1, 1));

    n += m;

    EXPECT_TRUE(n == Matrix({ {2, 4}, {6, 8} }));

    n -= m;

    EXPECT_TRUE(n == m);
}

TEST_F(MatrixTest, addMatrix) {
    Matrix m({ {1, 2}, {3,4} });
    Matrix n({ {1, 2}, {3,4} });
    Matrix o = m + n;
    EXPECT_TRUE(o == Matrix({ {2, 4}, {6, 8} }));
}

TEST_F(MatrixTest, subMatrix) {
    Matrix m({ {1, 2}, {3,4} });
    Matrix n({ {1, 2}, {3,4} });
    Matrix o = m - n;
    EXPECT_TRUE(o == Matrix({ {0, 0}, {0, 0} }));
}

TEST_F(MatrixTest, addScalar) {
    Matrix m({ {1, 2}, {3,4} });
    m += 1;
    EXPECT_TRUE(m == Matrix({ {2, 3}, {4, 5} }));
}

TEST_F(MatrixTest, subScalar) {
    Matrix m({ {1, 2}, {3,4} });
    m -= 1;
    EXPECT_TRUE(m == Matrix({ {0, 1}, {2, 3} }));
}

TEST_F(MatrixTest, multiply2dMatrix) {
    Matrix m({ {1, 2}, {3, 4} });
    Matrix n({ {1, 2}, {3, 4} });
    Matrix o = m * n;
    EXPECT_TRUE(o == Matrix({ {7, 10}, {15, 22} }));
}

TEST_F(MatrixTest, multiplyDifferentSizeMatrix) {
    Matrix m({ {1, 2, 3}, {3, -1, 0} });
    Matrix n({ {1, 2}, {3, 4}, {5, 6} });
    Matrix o = m * n;
    EXPECT_TRUE(o == Matrix({ {22, 28}, {0, 2} }));
    Matrix p = n * m;
    EXPECT_TRUE(p == Matrix({ {7, 0, 3}, {15, 2, 9}, {23, 4, 15} }));
}

TEST_F(MatrixTest, divideScalar) {
    Matrix m({ {1, 2}, {3,4} });
    m /= 2;
    EXPECT_TRUE(m == Matrix({ {0.5, 1}, {1.5, 2} }));
}

TEST_F(MatrixTest, determinantTest) {
    Matrix m({ {1, 2}, {3, 4} });
    EXPECT_EQ(-2, m.getDeterminant());
    Matrix n({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
    EXPECT_EQ(0, n.getDeterminant());
    Matrix o({ {1, 3, 5, 9}, {1, 3, 1, 7}, {4, 3, 9, 7}, {5, 2, 0, 9} });
    EXPECT_EQ(-376, o.getDeterminant());
}

TEST_F(MatrixTest, cofactor2dTest) {
    Matrix m({ {1, 2}, {3, 4} });
    EXPECT_EQ(4, m.cofactor(0, 0));
    EXPECT_EQ(-3, m.cofactor(0, 1));
    EXPECT_EQ(-2, m.cofactor(1, 0));
    EXPECT_EQ(1, m.cofactor(1, 1));
}

TEST_F(MatrixTest, cofactor3dTest) {
    Matrix n({ {1, 2, 1}, {6, -1, 0}, {-1, -2, -1} });
    EXPECT_EQ(1, n.cofactor(0, 0));
    EXPECT_EQ(6, n.cofactor(0, 1));
    EXPECT_EQ(-13, n.cofactor(0, 2));
    EXPECT_EQ(0, n.cofactor(1, 0));
    EXPECT_EQ(0, n.cofactor(1, 1));
    EXPECT_EQ(0, n.cofactor(1, 2));
    EXPECT_EQ(1, n.cofactor(2, 0));
    EXPECT_EQ(6, n.cofactor(2, 1));
    EXPECT_EQ(-13, n.cofactor(2, 2));
}

TEST_F(MatrixTest, adjugate2dTest) {
    Matrix m({ {1, 2}, {3, 4} });
    Matrix expected({ {4, -2}, {-3, 1} });
    Matrix actual = m.adjugate();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, adjugate3dTest) {
    Matrix m({ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} });
    Matrix expected({ {-24, 18, 5}, {20, -15, -4}, {-5, 4, 1} });
    Matrix actual = m.adjugate();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, inverse2dTest) {
    Matrix m({ {1, 2}, {3, 4} });
    Matrix actual = m.inverse();

    EXPECT_EQ(-2, actual(0, 0));
    EXPECT_EQ(1, actual(0, 1));
    EXPECT_EQ(1.5, actual(1, 0));
    EXPECT_EQ(-0.5, actual(1, 1));
}

TEST_F(MatrixTest, inverse3dTest) {
    Matrix m({ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} });
    const Matrix expected({ {-24, 18, 5}, {20, -15, -4}, {-5, 4, 1} });
    const Matrix actual = m.inverse();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, eigenvalues2dTest) {
    Matrix m({ {6, -1}, {2, 3} });
    const Matrix expected({ {5}, {4} });
    const auto actual = m.eigenvalues();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, eigenvalues3dDupTest) {
    Matrix m({ {2, 0, 0}, {1, 2, 1}, {-1, 0, 1} });
    const Matrix expected({ {1}, {2}, {2} });
    const auto actual = m.eigenvalues();
    EXPECT_EQ(2, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(actual.contains(2));
    EXPECT_TRUE(actual.contains(1));
}

TEST_F(MatrixTest, eigenvalues3dComplexTest) {
    Matrix m({ {33, -23, 9}, {22, 33, -23}, {19, 14, 50} });
    m /= 29;
    const Matrix expected({ {2} });
    const auto actual = m.eigenvalues();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, eigenvalues3dRealTest) {
    Matrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
    const auto actual = m.eigenvalues();
    EXPECT_EQ(3, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(actual.contains(3.0));
    EXPECT_TRUE(actual.contains(-5.0));
    EXPECT_TRUE(actual.contains(6.0));

    const auto actualVector1 = m.eigenvectorFor(3.0);
    expectNormalizedEqual(Matrix({ { -2 }, { 3 }, { 1 } }), actualVector1);
    const auto actualVector2 = m.eigenvectorFor(-5.0);
    expectNormalizedEqual(Matrix({ { -2 }, { -1 }, { 1 } }), actualVector2);
    const auto actualVector3 = m.eigenvectorFor(6.0);
    expectNormalizedEqual(Matrix({ { 1 }, { 6 }, { 16 } }), actualVector3);
}

TEST_F(MatrixTest, setRowTest) {
    Matrix m({ {1, 2}, {3, 4} });

    m.setRow(1, Matrix({{ 5, 6 }}));
    const Matrix expected = { {1, 2}, {5, 6} };
    expectEqual(expected, m);
}

TEST_F(MatrixTest, getRowTest) {
    Matrix m({ {1,2},{3,4} });
    const auto row = m.getRow(0);
    expectEqual(Matrix({ {1, 2} }), row);
}


TEST_F(MatrixTest, setColumnTest) {
    Matrix m({ {1, 2}, {3, 4} });

    m.setColumn(1, Matrix({ { 5 }, { 7 }}));
    const Matrix expected = { {1, 5}, {3, 7} };
    expectEqual(expected, m);
}

TEST_F(MatrixTest, getColumnTest) {
    Matrix m({ {1,2},{3,4} });
    const auto column = m.getColumn(1);
    expectEqual(Matrix({{2}, {4}}), column);
}

TEST_F(MatrixTest, toRowEchelonFormTest1) {
    Matrix m({ {1, 2, 3, 0}, {3, 4, 5, 0}, {4, 5, 6, 0} });
    m.toRowEchelonForm();
    const Matrix expected({ {1, 0, -1, 0}, {0, 1, 2, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
}

TEST_F(MatrixTest, toRowEchelonFormTest2) {
    Matrix m({ {-5, -4, 2, 0}, {-2, -2, 2, 0}, {4, 2, 2, 0} });
    m.toRowEchelonForm();
    const Matrix expected({ {1, 0, 2, 0}, {0, 1, -3, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
    m.toRowEchelonForm();
    expectEqual(expected, m);
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

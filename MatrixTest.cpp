#include <gtest/gtest.h>
#include "Matrix.h"

class MatrixTest : public ::testing::Test {
protected:
    void expectEqual(const Matrix& m, const Matrix& n) {
        EXPECT_EQ(m.rows(), n.rows());
        EXPECT_EQ(m.columns(), n.columns());
        for (Dimension i = 0; i < m.rows(); ++i) {
            for (Dimension j = 0; j < m.columns(); ++j) {
                EXPECT_EQ(m(i, j), n(i, j));
            }
        }
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
    Matrix m({{1, 2}, {3,4}});
    EXPECT_TRUE(m.isSquare());
    Matrix n = m;
    EXPECT_EQ(1, n(0, 0));
    EXPECT_EQ(2, n(0, 1));
    EXPECT_EQ(3, n(1, 0));
    EXPECT_EQ(4, n(1, 1));

    n += m;

    EXPECT_TRUE(n == Matrix({{2, 4}, {6, 8}}));

    n -= m;

    EXPECT_TRUE(n == m);
}

TEST_F(MatrixTest, addMatrix) {
    Matrix m({{1, 2}, {3,4}});
    Matrix n({{1, 2}, {3,4}});
    Matrix o = m + n;
    EXPECT_TRUE(o == Matrix({{2, 4}, {6, 8}}));
}

TEST_F(MatrixTest, subMatrix) {
    Matrix m({{1, 2}, {3,4}});
    Matrix n({{1, 2}, {3,4}});
    Matrix o = m - n;
    EXPECT_TRUE(o == Matrix({{0, 0}, {0, 0}}));
}

TEST_F(MatrixTest, addScalar) {
    Matrix m({{1, 2}, {3,4}});
    m += 1;
    EXPECT_TRUE(m == Matrix({{2, 3}, {4, 5}}));
}

TEST_F(MatrixTest, subScalar) {
    Matrix m({{1, 2}, {3,4}});
    m -= 1;
    EXPECT_TRUE(m == Matrix({{0, 1}, {2, 3}}));
}

TEST_F(MatrixTest, multiply2dMatrix) {
    Matrix m({{1, 2}, {3, 4}});
    Matrix n({{1, 2}, {3, 4}});
    Matrix o = m * n;
    EXPECT_TRUE(o == Matrix({{7, 10}, {15, 22}}));
}

TEST_F(MatrixTest, multiplyDifferentSizeMatrix) {
    Matrix m({{1, 2, 3}, {3, -1, 0}});
    Matrix n({{1, 2}, {3, 4}, {5, 6}});
    Matrix o = m * n;
    EXPECT_TRUE(o == Matrix({{22, 28}, {0, 2}}));
    Matrix p = n * m;
    EXPECT_TRUE(p == Matrix({{7, 0, 3}, {15, 2, 9}, {23, 4, 15}}));
}

TEST_F(MatrixTest, divideScalar) {
    Matrix m({{1, 2}, {3,4}});
    m /= 2;
    EXPECT_TRUE(m == Matrix({{0.5, 1}, {1.5, 2}}));
}

TEST_F(MatrixTest, determinantTest) {
    Matrix m({{1, 2}, {3, 4}});
    EXPECT_EQ(-2, m.determinant());
    Matrix n({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    EXPECT_EQ(0, n.determinant());
    Matrix o({{1, 3, 5, 9}, {1, 3, 1, 7}, {4, 3, 9, 7}, {5, 2, 0, 9}});
    EXPECT_EQ(-376, o.determinant());
}

TEST_F(MatrixTest, cofactor2dTest) {
    Matrix m({{1, 2}, {3, 4}});
    EXPECT_EQ(4, m.cofactor(0, 0));
    EXPECT_EQ(-3, m.cofactor(0, 1));
    EXPECT_EQ(-2, m.cofactor(1, 0));
    EXPECT_EQ(1, m.cofactor(1, 1));
}

TEST_F(MatrixTest, cofactor3dTest) {
    Matrix n({{1, 2, 1}, {6, -1, 0}, {-1, -2, -1}});
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
    Matrix m({{1, 2}, {3, 4}});
    Matrix expected({{4, -2}, {-3, 1}});
    Matrix actual = m.adjugate();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, adjugate3dTest) {
    Matrix m({{1, 2, 3}, {0, 1, 4}, {5, 6, 0}});
    Matrix expected({{-24, 18, 5}, {20, -15, -4}, {-5, 4, 1}});
    Matrix actual = m.adjugate();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, inverse2dTest) {
    Matrix m({{1, 2}, {3, 4}});
    Matrix actual = m.inverse();

    EXPECT_EQ(-2, actual(0, 0));
    EXPECT_EQ(1, actual(0, 1));
    EXPECT_EQ(1.5, actual(1, 0));
    EXPECT_EQ(-0.5, actual(1, 1));
}

TEST_F(MatrixTest, inverse3dTest) {
    Matrix m({{1, 2, 3}, {0, 1, 4}, {5, 6, 0}});
    Matrix expected({{-24, 18, 5}, {20, -15, -4}, {-5, 4, 1}});
    Matrix actual = m.inverse();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, eigenvalues2dTest) {
    const double EPSILON = 0.000000000000001;
    Matrix m({{6, -1}, {2, 3}});
    Matrix expected({{5}, {4}});
    auto actual = m.eigenvalues();
    EXPECT_EQ(expected.rows(), actual.rows());
    for (int rows = 0;rows < expected.rows(); rows++) {
        EXPECT_NEAR(expected(rows,0), actual(rows,0), EPSILON);
    }
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

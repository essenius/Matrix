#include <gtest/gtest.h>
#include "MatrixTest.h"

void MatrixTest::expectNormalizedEqual(const Matrix& expected, const Matrix& actual, const std::string& message, const double epsilon) const {
        expectEqual(expected.normalized(), actual.normalized(), message, epsilon);
};

TEST_F(MatrixTest, add) {
    Matrix m({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
    Matrix n({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
    Matrix expected({ {2, 4, 6}, {8, 10, 12}, {14, 16, 18} });
    Matrix actual = m + n;
    expectEqual(expected, actual, "add");
}

TEST_F(MatrixTest, multiplyScalar) {
    Matrix m({ {1, 2}, {3, 4} });
    m *= 2;
    expectEqual(Matrix({ {2, 4}, {6, 8} }), m);
    Matrix n = 2 * m;
    expectEqual(Matrix({ {4, 8}, {12, 16} }), n);
}

TEST_F(MatrixTest, multiply2dMatrix) {
    Matrix m({ {1, 2}, {3, 4} });
    Matrix n({ {1, 2}, {3, 4} });
    Matrix o = m * n;
    expectEqual(Matrix({ {7, 10}, {15, 22} }), o);
}

TEST_F(MatrixTest, multiplyDifferentSizeMatrix) {
    Matrix m({ {1, 2, 3}, {3, -1, 0} });
    Matrix n({ {1, 2}, {3, 4}, {5, 6} });
    Matrix o = m * n;
    expectEqual(Matrix({ {22, 28}, {0, 2} }), o, "m * n");
    Matrix p = n * m;
    expectEqual(Matrix({ {7, 0, 3}, {15, 2, 9}, {23, 4, 15} }), p, "n * m");
}

TEST_F(MatrixTest, determinant) {
    Matrix m({ {1, 2}, {3, 4} });
    EXPECT_EQ(-2, m.getDeterminant());
    Matrix n({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
    EXPECT_EQ(0, n.getDeterminant());
    Matrix o({ {1, 3, 5, 9}, {1, 3, 1, 7}, {4, 3, 9, 7}, {5, 2, 0, 9} });
    EXPECT_EQ(-376, o.getDeterminant());
}

TEST_F(MatrixTest, cofactor2d) {
    Matrix m({ {1, 2}, {3, 4} });
    EXPECT_EQ(4, m.getCofactor(0, 0));
    EXPECT_EQ(-3, m.getCofactor(0, 1));
    EXPECT_EQ(-2, m.getCofactor(1, 0));
    EXPECT_EQ(1, m.getCofactor(1, 1));
}

TEST_F(MatrixTest, cofactor3d) {
    Matrix n({ {1, 2, 1}, {6, -1, 0}, {-1, -2, -1} });
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
    Matrix m({ {1, 2}, {3, 4} });
    Matrix expected({ {4, -2}, {-3, 1} });
    Matrix actual = m.getAdjugate();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, adjugate3d) {
    Matrix m({ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} });
    Matrix expected({ {-24, 18, 5}, {20, -15, -4}, {-5, 4, 1} });
    Matrix actual = m.getAdjugate();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, inverse2d) {
    Matrix m({ {1, 2}, {3, 4} });
    Matrix actual = m.inverted();
    expectEqual(Matrix({ {-2, 1}, {1.5, -0.5} }), actual, "inverse");
    expectEqual(m, actual.inverted(), "inverse of inverse");
}

TEST_F(MatrixTest, inverse3d) {
    Matrix m({ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} });
    const Matrix expected({ {-24, 18, 5}, {20, -15, -4}, {-5, 4, 1} });
    const Matrix actual = m.inverted();
    expectEqual(expected, actual, "inverse");
    expectEqual(m, actual.inverted(), "inverse of inverse");
}

TEST_F(MatrixTest, normalize) {
    Matrix m({ {3, 4} });
    Matrix expected({ {0.6, 0.8} });
    Matrix actual = m.normalized();
    expectEqual(expected, actual, "normalize all positive");
    expectEqual(Matrix({{0}}), Matrix({{0}}).normalized(), "normalize zero");
    const Matrix m2({{-1, 4, -2, 2}});
    const Matrix expected2({{0.2, -0.8, 0.4, -0.4}});
    expectEqual(expected2, m2.normalized(), "normalize first negative");
}   

TEST_F(MatrixTest, transpose) {
    Matrix m({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
    const Matrix expected({ {1, 4, 7}, {2, 5, 8}, {3, 6, 9} });
    const Matrix actual = m.transposed();
    expectEqual(expected, actual, "Transpose");
    expectEqual(m, actual.transposed(), "Transposed transpose");
}

TEST_F(MatrixTest, toArray) {
    std::initializer_list<std::initializer_list<double>> input = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
    Matrix m(input);
    Array expected(input);
    auto actual = m.toArray();
    expectEqual(expected, actual, "toArray");
}

TEST_F(MatrixTest, minor3d) {
    Matrix m({ {3, 5, 0}, {2, -1, -7}, {6, -1, 5} });
    expectEqual(Matrix({{-1, -7}, {-1, 5}}), m.getMinor(0, 0));
    expectEqual(Matrix({{2, -7}, {6, 5}}), m.getMinor(0, 1));
    expectEqual(Matrix({{2, -1}, {6, -1}}), m.getMinor(0, 2));
    expectEqual(Matrix({{5, 0}, {-1, 5}}), m.getMinor(1, 0));
    expectEqual(Matrix({{3, 0}, {6, 5}}), m.getMinor(1, 1));
    expectEqual(Matrix({{3, 5}, {6, -1}}), m.getMinor(1, 2));
    expectEqual(Matrix({{5, 0}, {-1, -7}}), m.getMinor(2, 0));
    expectEqual(Matrix({{3, 0}, {2, -7}}), m.getMinor(2, 1));
    expectEqual(Matrix({{3, 5}, {2, -1}}), m.getMinor(2, 2));
}

TEST_F(MatrixTest, cofactor3dbis) {
    Matrix m({ {3, 5, 0}, {2, -1, -7}, {6, -1, 5} });
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
    auto m = Matrix(a);
    expectEqual(a, m);
}
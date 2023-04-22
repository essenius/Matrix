#include <gtest/gtest.h>
#include "BasicMatrixTest.h"

void BasicMatrixTest::expectNormalizedEqual(const BasicMatrix& expected, const BasicMatrix& actual, const std::string& message) {
        expectEqual(expected.normalize(), actual.normalize(), message);
};

TEST_F(BasicMatrixTest, multiply2dMatrix) {
    BasicMatrix m({ {1, 2}, {3, 4} });
    BasicMatrix n({ {1, 2}, {3, 4} });
    BasicMatrix o = m * n;
    expectEqual(BasicMatrix({ {7, 10}, {15, 22} }), o);
}

TEST_F(BasicMatrixTest, multiplyDifferentSizeMatrix) {
    BasicMatrix m({ {1, 2, 3}, {3, -1, 0} });
    BasicMatrix n({ {1, 2}, {3, 4}, {5, 6} });
    BasicMatrix o = m * n;
    expectEqual(BasicMatrix({ {22, 28}, {0, 2} }), o, "m * n");
    BasicMatrix p = n * m;
    expectEqual(BasicMatrix({ {7, 0, 3}, {15, 2, 9}, {23, 4, 15} }), p, "n * m");
}

TEST_F(BasicMatrixTest, determinant) {
    BasicMatrix m({ {1, 2}, {3, 4} });
    EXPECT_EQ(-2, m.getDeterminant());
    BasicMatrix n({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
    EXPECT_EQ(0, n.getDeterminant());
    BasicMatrix o({ {1, 3, 5, 9}, {1, 3, 1, 7}, {4, 3, 9, 7}, {5, 2, 0, 9} });
    EXPECT_EQ(-376, o.getDeterminant());
}

TEST_F(BasicMatrixTest, cofactor2d) {
    BasicMatrix m({ {1, 2}, {3, 4} });
    EXPECT_EQ(4, m.cofactor(0, 0));
    EXPECT_EQ(-3, m.cofactor(0, 1));
    EXPECT_EQ(-2, m.cofactor(1, 0));
    EXPECT_EQ(1, m.cofactor(1, 1));
}

TEST_F(BasicMatrixTest, cofactor3d) {
    BasicMatrix n({ {1, 2, 1}, {6, -1, 0}, {-1, -2, -1} });
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

TEST_F(BasicMatrixTest, adjugate2d) {
    BasicMatrix m({ {1, 2}, {3, 4} });
    BasicMatrix expected({ {4, -2}, {-3, 1} });
    BasicMatrix actual = m.adjugate();
    expectEqual(expected, actual);
}

TEST_F(BasicMatrixTest, adjugate3d) {
    BasicMatrix m({ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} });
    BasicMatrix expected({ {-24, 18, 5}, {20, -15, -4}, {-5, 4, 1} });
    BasicMatrix actual = m.adjugate();
    expectEqual(expected, actual);
}

TEST_F(BasicMatrixTest, inverse2d) {
    BasicMatrix m({ {1, 2}, {3, 4} });
    BasicMatrix actual = m.inverse();
    expectEqual(BasicMatrix({ {-2, 1}, {1.5, -0.5} }), actual, "inverse");
    expectEqual(m, actual.inverse(), "inverse of inverse");
}

TEST_F(BasicMatrixTest, inverse3d) {
    BasicMatrix m({ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} });
    const BasicMatrix expected({ {-24, 18, 5}, {20, -15, -4}, {-5, 4, 1} });
    const BasicMatrix actual = m.inverse();
    expectEqual(expected, actual, "inverse");
    expectEqual(m, actual.inverse(), "inverse of inverse");
}

TEST_F(BasicMatrixTest, transpose) {
    BasicMatrix m({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
    const BasicMatrix expected({ {1, 4, 7}, {2, 5, 8}, {3, 6, 9} });
    const BasicMatrix actual = m.transpose();
    expectEqual(expected, actual, "Transpose");
    expectEqual(m, actual.transpose(), "Transposed transpose");
}
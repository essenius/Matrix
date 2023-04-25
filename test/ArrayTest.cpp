#include <gtest/gtest.h>
#include "ArrayTest.h"

bool ArrayTest::isEqual(const Array& expected, const Array& actual, const double epsilon) {
    if (expected.rows() != actual.rows()) return false;
    if (expected.columns() != actual.columns()) return false;

    for (Dimension row = 0; row < expected.rows(); row++) {
        for (Dimension column = 0; column < expected.columns(); column++) {
            auto difference = expected(row, column) - actual(row, column);
            if (abs(difference) > epsilon) return false;
        }
    }
    return true;
}

// repeat of above, but now asserting equality and showing details
void ArrayTest::expectEqual(const Array& expected, const Array& actual, const std::string& message, const double epsilon) {
    EXPECT_EQ(expected.rows(), actual.rows()) << message << " rows";
    EXPECT_EQ(expected.columns(), actual.columns()) << message << " columns";
    for (Dimension row = 0; row < expected.rows(); row++) {
        for (Dimension column = 0; column < expected.columns(); column++) {
            EXPECT_NEAR(expected(row, column), actual(row, column), epsilon) << message <<  "(" << row << ", " << column << ")";
        }
    }
}

bool ArrayTest::contains(const Array& matrix, const double value, const double epsilon) {
    for (Dimension row = 0; row < matrix.rows(); row++) {
        for (Dimension column = 0; column < matrix.columns(); column++) {
            if (fabs(matrix(row, column) - value) <= epsilon) return true;
        }
    }
    return false;
}

TEST_F(ArrayTest, initArray) {
    Array m(2, 2);
    m(0, 0) = 1;
    EXPECT_EQ(1, m(0, 0));
    m(0, 1) = 2;
    EXPECT_EQ(2, m(0, 1));
    m(1, 0) = 3;
}

TEST_F(ArrayTest, copyAddArray) {
    Array m({ {1, 2}, {3,4} });
    EXPECT_TRUE(m.isSquare());
    Array n = m;
    EXPECT_EQ(1, n(0, 0));
    EXPECT_EQ(2, n(0, 1));
    EXPECT_EQ(3, n(1, 0));
    EXPECT_EQ(4, n(1, 1));

    n += m;
    EXPECT_TRUE(n == Array({ {2, 4}, {6, 8} }));
    EXPECT_FALSE(n == m);

    n -= m;
    EXPECT_TRUE(n == m);
}

TEST_F(ArrayTest, addArray) {
    Array m({ {1, 2}, {3,4} });
    Array n({ {1, 2}, {3,4} });
    Array o = m + n;
    EXPECT_TRUE(o == Array({ {2, 4}, {6, 8} }));
}

TEST_F(ArrayTest, subtractArray) {
    Array m({ {1, 2}, {3,4} });
    Array n({ {1, 2}, {3,4} });
    Array o = m - n;
    EXPECT_TRUE(o == Array({ {0, 0}, {0, 0} }));
}

TEST_F(ArrayTest, multiplyArray) {
    Array m({ {1, 2}, {3,4} });
    Array n({ {1, 2}, {3,4} });
    Array o = m * n;
    EXPECT_TRUE(o == Array({ {1, 4}, {9, 16} }));
}

TEST_F(ArrayTest, addScalar) {
    Array m({ {1, 2}, {3,4} });
    m += 1;
    EXPECT_TRUE(m == Array({ {2, 3}, {4, 5} }));
}

TEST_F(ArrayTest, subtractScalar) {
    Array m({ {1, 2}, {3,4} });
    m -= 1;
    EXPECT_TRUE(m == Array({ {0, 1}, {2, 3} }));
}

TEST_F(ArrayTest, multiplyScalar) {
    Array m({ {1, 2}, {3,4} });
    m *= 2;
    EXPECT_TRUE(m == Array({ {2, 4}, {6, 8} }));
}

TEST_F(ArrayTest, divideScalar) {
    Array m({ {1, 2}, {3,4} });
    m /= 2;
    expectEqual(Array({ {0.5, 1}, {1.5, 2} }), m);
}

TEST_F(ArrayTest, getColumn) {
    Array m({ {1,2},{3,4} });
    const auto column = m.getColumn(1);
    expectEqual(Array({{2}, {4}}), column);
}

TEST_F(ArrayTest, getRow) {
    Array m({ {1,2},{3,4} });
    const auto row = m.getRow(0);
    expectEqual(Array({ {1, 2} }), row);
}

TEST_F(ArrayTest, setColumn) {
    Array m({ {1, 2}, {3, 4} });

    m.setColumn(1, Array({ { 5 }, { 7 }}));
    const Array expected = { {1, 5}, {3, 7} };
    expectEqual(expected, m);
}

TEST_F(ArrayTest, setColumnCount) {
    Array m({ {1, 2}, {3, 4} });
    auto n = m;
    n.setColumnCount(4);
    const Array expected = { {1, 2, 0, 0}, {3, 4, 0, 0} };
    expectEqual(expected, n, "setColumnCount larger");
    n.setColumnCount(2);
    expectEqual(m, n, "setColumnCount smaller");
}


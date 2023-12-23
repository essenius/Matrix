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
#include "ArrayTest.h"

namespace RixMatrixTest {
    using RixMatrix::Array;
    using RixMatrix::Dimension;

	bool ArrayTest::isEqual(const Array& expected, const Array& actual, const double epsilon) {
        if (expected.rowCount() != actual.rowCount()) return false;
        if (expected.columnCount() != actual.columnCount()) return false;

        for (Dimension row = 0; row < expected.rowCount(); row++) {
            for (Dimension column = 0; column < expected.columnCount(); column++) {
                const auto difference = expected(row, column) - actual(row, column);
                if (abs(difference) > epsilon) return false;
            }
        }
        return true;
    }

    // repeat of above, but now asserting equality and showing details
    void ArrayTest::expectEqual(const Array& expected, const Array& actual, const std::string& message, const double epsilon) {
        EXPECT_EQ(expected.rowCount(), actual.rowCount()) << message << " rows";
        EXPECT_EQ(expected.columnCount(), actual.columnCount()) << message << " columns";
        for (Dimension row = 0; row < expected.rowCount(); row++) {
            for (Dimension column = 0; column < expected.columnCount(); column++) {
                EXPECT_NEAR(expected(row, column), actual(row, column), epsilon) << message << "(" << row << ", " << column << ")";
            }
        }
    }

    bool ArrayTest::contains(const Array& matrix, const double value, const double epsilon) {
        for (Dimension row = 0; row < matrix.rowCount(); row++) {
            for (Dimension column = 0; column < matrix.columnCount(); column++) {
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
        EXPECT_EQ(3, m(1, 0));

    }

    TEST_F(ArrayTest, copyAddArray) {
        const Array m({ {1, 2}, {3,4} });
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
        const Array m({ {1, 2}, {3,4} });
        const Array n({ {1, 2}, {3,4} });
        const Array o = m + n;
        EXPECT_TRUE(o == Array({ {2, 4}, {6, 8} }));
    }

    TEST_F(ArrayTest, subtractArray) {
        const Array m({ {1, 2}, {3,4} });
        const Array n({ {1, 2}, {3,4} });
        const Array o = m - n;
        EXPECT_TRUE(o == Array({ {0, 0}, {0, 0} }));
    }

    TEST_F(ArrayTest, multiplyArray) {
        const Array m({ {1, 2}, {3,4} });
        const Array n({ {1, 2}, {3,4} });
        const Array o = m * n;
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
        const Array n = 2 * m;
        EXPECT_TRUE(n == Array({ {4, 8}, {12, 16} }));
    }

    TEST_F(ArrayTest, divideScalar) {
        Array m({ {1, 2}, {3,4} });
        m /= 2;
        expectEqual(Array({ {0.5, 1}, {1.5, 2} }), m);
    }

    TEST_F(ArrayTest, getColumn) {
        const Array m({ {1,2},{3,4} });
        const auto column = m.getColumn(1);
        expectEqual(Array({ {2}, {4} }), column);
    }

    TEST_F(ArrayTest, getRow) {
        const Array m({ {1,2},{3,4} });
        const auto row = m.getRow(0);
        expectEqual(Array({ {1, 2} }), row);
    }

    TEST_F(ArrayTest, pow2) {
        const Array m({ {1,2},{3,4} });
        const auto actual = m.pow2();
        expectEqual(Array({ {1, 4}, {9, 16} }), actual);
    }
    TEST_F(ArrayTest, setColumn) {
        Array m({ {1, 2}, {3, 4} });

        m.setColumn(1, Array({ { 5 }, { 7 } }));
        const Array expected({ {1, 5}, {3, 7} });
        expectEqual(expected, m);
    }

    TEST_F(ArrayTest, setColumnCount) {
        const Array m({ {1, 2}, {3, 4} });
        auto n = m;
        n.setColumnCount(4);
        const Array expected({ {1, 2, 0, 0}, {3, 4, 0, 0} });
        expectEqual(expected, n, "setColumnCount larger");
        n.setColumnCount(2);
        expectEqual(m, n, "setColumnCount smaller");
    }

    TEST_F(ArrayTest, setRow) {
        Array m({ {1, 2}, {3, 4} });

        m.setRow(1, Array({ { 5, 7 } }));
        const Array expected({ {1, 2}, {5, 7} });
        expectEqual(expected, m);
    }

    TEST_F(ArrayTest, setRowCount) {
        const Array m({ {1, 2}, {3, 4} });
        auto n = m;
        n.setRowCount(4);
        const Array expected({ {1, 2}, {3, 4}, {0, 0}, {0, 0} });
        expectEqual(expected, n, "setRowCount larger");
        n.setRowCount(2);
        expectEqual(m, n, "setRowCount smaller");
    }
    TEST_F(ArrayTest, DifferentSizesNotEqual) {
        const Array m({ {1, 2}, {3, 4} });
        const Array n({ {1} });
        EXPECT_FALSE(m == n);
    }

    TEST_F(ArrayTest, SetColumnScalar) {
        Array m({ {1, 2}, {3, 4} });
        m.setColumn(1, 5);
        const Array expected({ {1, 5}, {3, 5} });
        expectEqual(expected, m);
    }

    TEST_F(ArrayTest, SetRowScalar) {
        Array m({ {1, 2}, {3, 4} });
        m.setRow(1, 5);
        const Array expected({ {1, 2}, {5, 5} });
        expectEqual(expected, m);
    }

    TEST_F(ArrayTest, GetEpsilon) {
        const Array m({ {Array::Epsilon / 2} });
        const Array n({ {Array::Epsilon / 3} });
        expectEqual(m, n);
    }

    TEST_F(ArrayTest, SwapRows) {
        Array m({ {1, 2}, {3, 4}, {5, 6} });
        EXPECT_FALSE(m.isSquare()) << "Not square";
        m.swapRows(0, 2);
        const Array expected({ {5,6}, {3, 4}, {1, 2} });
        expectEqual(expected, m);
        m.swapRows(1, 1);
        expectEqual(expected, m);
    }

    TEST_F(ArrayTest, SwapColumns) {
        Array m({ {1, 2}, {3, 4}, {5, 6} });
        m.swapColumns(0, 1);
        const Array expected({ {2, 1}, {4, 3}, {6, 5} });
        expectEqual(expected, m);
        m.swapColumns(1, 1);
        expectEqual(expected, m);
    }

    TEST_F(ArrayTest, transpose) {
        const Array m({ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} });
        const Array expected({ {1, 4, 7}, {2, 5, 8}, {3, 6, 9} });
        const Array actual = m.transposed();
        expectEqual(expected, actual, "Transpose");
        expectEqual(m, actual.transposed(), "Transposed transpose");
    }

    TEST_F(ArrayTest, AssertColumns) {
        Array m({ {1, 2} });
        EXPECT_DEATH(m.swapColumns(0, 2), "Assertion failed: .*column1 < _columns && column2 < _columns");
        EXPECT_DEATH(m.swapColumns(2, 0), "Assertion failed: .*column1 < _columns && column2 < _columns");
        EXPECT_DEATH(m.getColumn(2), "Assertion failed: .*column < _columns");

        // testing the non-const operator()
        EXPECT_DEATH(m(0, 2), "Assertion failed: .*row < _rows && column < _columns");

        EXPECT_DEATH(m.me(0, 2), "Assertion failed: .*row < _rows && column < _columns");

        const Array n({ {1, 2}, {3, 4} });

        // testing the const operator() and operator[]
        EXPECT_DEATH(n(0, 3), "Assertion failed: .*row < _rows && column < _columns");

        EXPECT_DEATH(m.setColumn(2, n), "Assertion failed: .*column < _columns");
        EXPECT_DEATH(m.setColumn(1, n), "Assertion failed: .*input.rowCount\\(\\) == _rows");
        EXPECT_DEATH(m.setColumn(2, 1), "Assertion failed: .*column < _columns");
    }

    TEST_F(ArrayTest, AssertRows) {
        Array m({ {1, 2} });
        EXPECT_DEATH(m.swapRows(0, 2), "Assertion failed: .*row1 < _rows && row2 < _rows");
        EXPECT_DEATH(m.swapRows(2, 0), "Assertion failed: .*row1 < _rows && row2 < _rows");
        EXPECT_DEATH(m.getRow(2), "Assertion failed: .*row < _rows");
        // testing the non-const operator()
        EXPECT_DEATH(m(2, 0), "Assertion failed: .*row < _rows && column < _columns");

        EXPECT_DEATH(m.me(2, 0), "Assertion failed: .*row < _rows && column < _columns");
        const Array n({ {1, 2}, {3, 4} });

        // testing the const operator() and operator[]
        EXPECT_DEATH(n(3, 0), "Assertion failed: .*row < _rows && column < _columns");

        const Array o({ {1} });
        EXPECT_DEATH(m.setRow(2, n), "Assertion failed: .*row < _rows");
        EXPECT_DEATH(m.setRow(1, o), "Assertion failed: .*input.columnCount\\(\\) == _columns");
        EXPECT_DEATH(m.setRow(2, 1), "Assertion failed: .*row < _rows");
    }

    TEST_F(ArrayTest, AssertOther) {
        Array m({ {1, 2} });
        EXPECT_DEATH(m[2], "Assertion failed: .*cell < _arraySize");

        const Array n({ {1, 2}, {3, 4} });

        EXPECT_DEATH(m += n, "Assertion failed: .*other.sizeIsEqual\\(\\*this\\)");
        EXPECT_DEATH(m -= n, "Assertion failed: .*other.sizeIsEqual\\(\\*this\\)");
        EXPECT_DEATH(m *= n, "Assertion failed: .*other.sizeIsEqual\\(\\*this\\)");

        EXPECT_DEATH(n[4], "Assertion failed: .*cell < _arraySize");
    }
}
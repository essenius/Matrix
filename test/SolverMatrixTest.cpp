// Copyright 2023-2024 Rik Essenius
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
#include <SolverMatrix.h>

namespace RixMatrixTest {
    using RixMatrix::SolverMatrix;
    using RixMatrix::Dimension;

    class SolverMatrixTest : public MatrixTest {
    protected:
	    static bool normalizedContains(const Matrix& expected, const Matrix& actual) {
            const auto expectedNormalized = expected.normalized();
            for (Dimension column = 0; column < actual.columnCount(); column++) {
                // actual should already be normalized
                Array vector = actual.getColumn(column);
                if (isEqual(expectedNormalized, vector, SolverMatrix::EigenEpsilon)) return true;
            }
            return false;
        }
    };

    // *** Row echelon format ***

    TEST_F(SolverMatrixTest, toRowEchelonForm10x) {
        SolverMatrix m({ {10, 20, 30, 0}, {30, 40, 50, 0}, {40, 50, 60, 0} });
        const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
        const Array expected({ {1, 0, 0.5, 0}, {0, 1, 0.5, 0}, {0, 0, 0, 0} });
        expectEqual(expected, m);
        expectEqual(Array({ {0,1,0,0}, {0,0,1,0}, {1,0,0,0}, {0,0,0,1} }), permutation);
    }


    TEST_F(SolverMatrixTest, toRowEchelonForm1) {
        SolverMatrix m({ {1, 2, 3, 0}, {3, 4, 5, 0}, {4, 5, 6, 0} });
        const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
        const Array expected({ {1, 0, 0.5, 0}, {0, 1, 0.5, 0}, {0, 0, 0, 0} });
        expectEqual(expected, m);
        expectEqual(Array({ {0,1,0,0}, {0,0,1,0}, {1,0,0,0}, {0,0,0,1} }), permutation);
        m.toReducedRowEchelonFormWithPivot();
        expectEqual(expected, m);
    }

    TEST_F(SolverMatrixTest, toRowEchelonForm2) {
        SolverMatrix m({ {-5, -4, 2, 0}, {-2, -2, 2, 0}, {4, 2, 2, 0} });
        const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
        const Array expected({ {1, 0, 2.0 / 3.0, 0}, {0, 1, -1.0 / 3.0, 0}, {0, 0, 0, 0} });
        expectEqual(expected, m);
        expectEqual(Array({ {1,0,0,0}, {0,0,1,0}, {0,1,0,0}, {0,0,0,1} }), permutation);
    }

    TEST_F(SolverMatrixTest, toRowEchelonForm3) {
        SolverMatrix m({ {0, 0, 0}, {1, 0, 1}, {0, 1, -1} });
        const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
        const Array expected({ {1, 0, 1}, {0, 1, -1}, {0, 0, 0} });
        expectEqual(expected, m);
        expectEqual(Array({ {1, 0, 0}, {0, 1, 0}, {0, 0, 1} }), permutation);
    }

    TEST_F(SolverMatrixTest, toRowEchelonForm4) {
        SolverMatrix m({ {1, 0, 0}, {0, 2, 1}, {0, 0, 2} });
        const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
        const Array expected({ {1, 0, 0}, {0, 1, 0}, {0, 0, 1} });
        expectEqual(expected, m);
        expectEqual(Array({ {0, 0, 1}, {1, 0, 0}, {0, 1, 0} }), permutation);
    }

    TEST_F(SolverMatrixTest, toRowEchelonFormSpecial) {
        SolverMatrix m({ {8192.10277, -0.0341611, 8192.03422}, {0.0686703, -6.3222615e-7, 0.0686707}, {8192.03422, -0.03416085, 8192.10277} });
        const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
        const Array expected({ {1, 0, -3.9e-6}, {0, 1, -2.62e-7}, {0, 9.74e-7, -3.46e-7} });
        expectEqual(expected, m, "toRowEchelonFormSpecial", 1e-7);
        expectEqual(Array({ {1, 0, 0}, {0, 0, 1}, {0, 1, 0} }), permutation);
    }

    // *** Null spaces and free variables ***

    TEST_F(SolverMatrixTest, oneFreeVariableColumn) {
        const SolverMatrix m({ {0,1,1}, {0,0,1}, {0,0,0} });
        const auto actual = m.getFreeVariables();
        EXPECT_EQ(1, actual.size());
        EXPECT_EQ(0, actual[0]);
    }

    TEST_F(SolverMatrixTest, twoFreeVariableColumns) {
        const SolverMatrix m({ {0,1,0}, {0,0,0}, {0,0,0} });
        const auto actual = m.getFreeVariables();
        EXPECT_EQ(2, actual.size());
        EXPECT_EQ(0, actual[0]);
        EXPECT_EQ(2, actual[1]);
    }

    TEST_F(SolverMatrixTest, noFreeVariableColumns) {
        const SolverMatrix m({ {1,0,0}, {0,1,0}, {0,0,1} });
        const auto actual = m.getFreeVariables();
        EXPECT_EQ(0, actual.size());
    }

    TEST_F(SolverMatrixTest, nullSpaceNoFreeVariables) {
        const SolverMatrix m({ {1, 0, 0}, {0, 1, 0}, {0, 0, 1} });
        const auto actual = m.getNullSpace();
        EXPECT_EQ(0, actual.columnCount());
        EXPECT_EQ(3, actual.rowCount());
    }

    TEST_F(SolverMatrixTest, nullSpaceForOneFreeVariable) {
        const SolverMatrix m({ {1, 0, 0}, {0, 1, 0}, {0, 0, 0} });
        const auto actual = m.getNullSpace();
        EXPECT_EQ(1, actual.columnCount());
        expectNormalizedEqual(Matrix({ { 0 }, { 0 }, { 1 } }), actual);
    }

    TEST_F(SolverMatrixTest, nullSpaceForTwoFreeVariable) {
        const SolverMatrix m({ { 0, 1, 0}, {0, 0, 0}, {0, 0, 0} });
        const auto actual = m.getNullSpace();

        EXPECT_EQ(2, actual.columnCount());
        const auto expected1 = Matrix({ { 1 }, { 0 }, { 0 } });
        const auto expected2 = Matrix({ { 0 }, { 0 }, { 1 } });
        expectNormalizedEqual(expected1, Matrix(actual.getColumn(0)), "expected1");
        expectNormalizedEqual(expected2, Matrix(actual.getColumn(1)), "expected2");
    }

    TEST_F(SolverMatrixTest, nullSpaceSpecial) {
        const SolverMatrix m({ {1, -3.9e-6, 0}, {0, -2.62e-7, 1}, {0, -3.46e-7, 9.74e-7} });
        const auto actual = m.getNullSpace();

        EXPECT_EQ(1, actual.columnCount());
        const auto outcome = m * actual;
        const auto expected = Matrix({ { 3.9e-6 }, { 1 }, { 3.46e-7 } });
        expectNormalizedEqual(expected, actual, "null space", SolverMatrix::EigenEpsilon);
        expectEqual(Array({ {0},{0},{0} }), outcome, "matrix * null space", SolverMatrix::EigenEpsilon);
    }

    // *** Eigenvalues, eigen vectors ***

    TEST_F(SolverMatrixTest, getEigenvalues1d) {
        const SolverMatrix m({ {6} });
        const auto actual = m.getEigenvalues();
        expectEqual(m, actual);
    }

    TEST_F(SolverMatrixTest, getEigenvalues2d) {
        const SolverMatrix m({ {6, -1}, {2, 3} });
        const Matrix expected({ {5}, {4} });
        const auto actual = m.getEigenvalues();
        expectEqual(expected, actual);
    }

    TEST_F(SolverMatrixTest, getEigenvalues2dNoSolution) {
        const SolverMatrix m({ {0, 1}, {-1, 0} });
        const Matrix expected(0, 0);
        const auto actual = m.getEigenvalues();
        expectEqual(expected, actual);
    }

    TEST_F(SolverMatrixTest, getEigenvalues3dDup) {
        const SolverMatrix m({ {2, 0, 0}, {1, 2, 1}, {-1, 0, 1} });
        const auto actual = m.getEigenvalues();
        EXPECT_EQ(2, actual.rowCount());
        EXPECT_EQ(1, actual.columnCount());
        EXPECT_TRUE(contains(actual, 2));
        EXPECT_TRUE(contains(actual, 1));
    }

    TEST_F(SolverMatrixTest, getEigenvalues3dOne) {
        // trace(M) = 0, trace(M^2) = 0, determinant = 0
        const SolverMatrix m({ {0, 0, 1}, {0, 0, -1}, {1, 1, 0} });
        const auto actual = m.getEigenvalues();
        EXPECT_EQ(1, actual.rowCount());
        EXPECT_EQ(1, actual.columnCount());
        EXPECT_TRUE(contains(actual, 0));
    }

    TEST_F(SolverMatrixTest, getEigenvalues3dComplex) {
        SolverMatrix m({ {33, -23, 9}, {22, 33, -23}, {19, 14, 50} });
        m /= 29;
        const Array expected({ {2} });
        const auto actual = m.getEigenvalues();
        expectEqual(expected, actual);
    }

    TEST_F(SolverMatrixTest, getEigenvectorsFor3dReal) {
        const SolverMatrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
        const auto actual = m.getEigenvalues();
        EXPECT_EQ(3, actual.rowCount());
        EXPECT_EQ(1, actual.columnCount());
        EXPECT_TRUE(contains(actual, 3.0));
        EXPECT_TRUE(contains(actual, -5.0));
        EXPECT_TRUE(contains(actual, 6.0));

        const auto actualVector1 = m.getEigenvectorFor(3.0);
        expectNormalizedEqual(Matrix({ { -2 }, { 3 }, { 1 } }), actualVector1, "eigenvector 1");
        const auto actualVector2 = m.getEigenvectorFor(-5.0);
        expectNormalizedEqual(Matrix({ { -2 }, { -1 }, { 1 } }), actualVector2, "eigenvector 2");
        const auto actualVector3 = m.getEigenvectorFor(6.0);
        expectNormalizedEqual(Matrix({ { 1 }, { 6 }, { 16 } }), actualVector3, "eigenvector 3");
    }

    TEST_F(SolverMatrixTest, eigenvectorsForTwoFreeVariables) {
        const SolverMatrix m({ {1, 0, 0}, {0, 0, 0}, {0, 0, 1} });
        const auto actual = m.getEigenvalues();
        EXPECT_EQ(2, actual.rowCount());
        EXPECT_EQ(1, actual.columnCount());
        EXPECT_TRUE(contains(actual, 1.0));
        EXPECT_TRUE(contains(actual, 0.0));

        const auto actualVectors1 = m.getEigenvectorFor(1.0);
        EXPECT_EQ(2, actualVectors1.columnCount());
        const auto expected1 = Matrix({ { 1 }, { 0 }, { 0 } });
        const auto expected2 = Matrix({ { 0 }, { 0 }, { 1 } });
        expectNormalizedEqual(expected1, Matrix(actualVectors1.getColumn(0)), "actualVector1-1");
        expectNormalizedEqual(expected2, Matrix(actualVectors1.getColumn(1)), "actualVector1-2");

        const auto actualVector2 = m.getEigenvectorFor(0.0);
        expectNormalizedEqual(Matrix({ { 0 }, { 1 }, { 0 } }), actualVector2, "actualVector2");
    }

    TEST_F(SolverMatrixTest, getEigenVectorsTwoFreeVariables) {
        const SolverMatrix m({ {1, 0, 0}, {0, 0, 0}, {0, 0, 1} });
        const auto actual = m.getEigenvectors();
        EXPECT_EQ(3, actual.rowCount());
        EXPECT_EQ(3, actual.columnCount());

        const auto expected = Matrix({ {0, 1, 0}, {1, 0, 0}, {0, 0, 1} }).transposed();
        for (int i = 0; i < 3; ++i) {
            expectNormalizedEqual(Matrix(expected.getColumn(i)), Matrix(actual.getColumn(i)), "column " + std::to_string(i));
        }
    }

    TEST_F(SolverMatrixTest, getEigenvectors3dReal) {
        const SolverMatrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
        const auto actual = m.getEigenvectors();
        EXPECT_EQ(3, actual.rowCount());
        EXPECT_EQ(3, actual.columnCount());

        EXPECT_TRUE(normalizedContains(Matrix({ { -2 }, { 3 }, { 1 } }), actual)) << "eigenvector 1";
        EXPECT_TRUE(normalizedContains(Matrix({ { -2 }, { -1 }, { 1 } }), actual)) << "eigenvector 2";
        EXPECT_TRUE(normalizedContains(Matrix({ { 1 }, { 6 }, { 16 } }), actual)) << "eigenvector 3";
    }

    TEST_F(SolverMatrixTest, assignMatrix) {
        const Matrix m({ {1, 2}, {3, 4} });
        const auto sm = SolverMatrix(m);
        expectEqual(m, sm);
    }

#ifdef DEBUG
    TEST_F(SolverMatrixTest, assertTest) {
        SolverMatrix m({ {1, 2} });
        ASSERT_DEATH(m *= SolverMatrix({ {1} }), "Assertion.*failed");
        ASSERT_DEATH(m.getEigenvalues(), "Assertion.*failed");
        ASSERT_DEATH(m.getEigenvectorFor(1), "Assertion.*failed");
    }
#endif
}
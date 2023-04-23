#include <gtest/gtest.h>
#include "MatrixTest.h"
#include "SolverMatrix.h"

class SolverMatrixTest : public MatrixTest {
protected:
    bool normalizedContains(const Matrix& expected, const Matrix& actual) {
        auto expectedNormalized = expected.normalize();
        for (Dimension column = 0; column < actual.columns(); column++) {
            // actual should already be normalized
            Array vector = actual.getColumn(column);
            if (isEqual(expectedNormalized, vector)) return true;
        }
        return false;
    }
};

// *** Row echelon format ***

TEST_F(SolverMatrixTest, toRowEchelonForm1) {
    SolverMatrix m({ {1, 2, 3, 0}, {3, 4, 5, 0}, {4, 5, 6, 0} });
    m.toRowEchelonForm();
    const Array expected({ {1, 0, -1, 0}, {0, 1, 2, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
}

TEST_F(SolverMatrixTest, toRowEchelonForm2) {
    SolverMatrix m({ {-5, -4, 2, 0}, {-2, -2, 2, 0}, {4, 2, 2, 0} });
    m.toRowEchelonForm();
    const Array expected({ {1, 0, 2, 0}, {0, 1, -3, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
    m.toRowEchelonForm();
    expectEqual(expected, m);
}

TEST_F(SolverMatrixTest, toRowEchelonForm3) {
    SolverMatrix m({ {0, 0, 0}, {1, 0, 1}, {0, 1, -1} });
    m.toRowEchelonForm();
    const Array expected({ {1, 0, 1}, {0, 1, -1}, {0, 0, 0} });
    expectEqual(expected, m);
}

TEST_F(SolverMatrixTest, toRowEchelonForm4) {
    SolverMatrix m({ {1, 0, 0}, {0, 2, 1}, {0, 0, 2} });
    m.toRowEchelonForm();
    const Array expected({ {1, 0, 0}, {0, 1, 0}, {0, 0, 1} });
    expectEqual(expected, m);
}

// *** Null spaces and free variables ***

TEST_F(SolverMatrixTest, oneFreeVariableColumn) {
    SolverMatrix m({{0,1,1}, {0,0,1}, {0,0,0}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(1, actual.size());
    EXPECT_EQ(0, actual[0]);
}

TEST_F(SolverMatrixTest, twoFreeVariableColumns) {
    SolverMatrix m({{0,1,0}, {0,0,0}, {0,0,0}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(2, actual.size());
    EXPECT_EQ(0, actual[0]);
    EXPECT_EQ(2, actual[1]);
}

TEST_F(SolverMatrixTest, noFreeVariableColumns) {
    SolverMatrix m({{1,0,0}, {0,1,0}, {0,0,1}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(0, actual.size());
}

TEST_F(SolverMatrixTest, nullSpaceForOneFreeVariable) {
    SolverMatrix m({ {1, 0, 0}, {0, 1, 0}, {0, 0, 0} });
    const auto actual = m.getNullSpace();
    EXPECT_EQ(1, actual.columns());
    expectNormalizedEqual(Matrix({ { 0 }, { 0 }, { 1 } }), actual);
}

TEST_F(SolverMatrixTest, nullSpaceForTwoFreeVariable) {
    SolverMatrix m({ { 0, 1, 0}, {0, 0, 0}, {0, 0, 0} });
    const auto actual = m.getNullSpace();
    
    EXPECT_EQ(2, actual.columns());
    const auto expected1 = Matrix({ { 1 }, { 0 }, { 0 } });
    const auto expected2 = Matrix({ { 0 }, { 0 }, { 1 } });
    expectNormalizedEqual(expected1, actual.getColumn(0), "expected1");
    expectNormalizedEqual(expected2, actual.getColumn(1), "expected2"); 
}

// *** Eigenvalues, eigen vectors ***

TEST_F(SolverMatrixTest, getEigenvalues2d) {
    SolverMatrix m({ {6, -1}, {2, 3} });
    const Matrix expected({ {5}, {4} });
    const auto actual = m.getEigenvalues();
    expectEqual(expected, actual);
}

TEST_F(SolverMatrixTest, getEigenvalues3dDup) {
    SolverMatrix m({ {2, 0, 0}, {1, 2, 1}, {-1, 0, 1} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(2, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(contains(actual, 2));
    EXPECT_TRUE(contains(actual, 1));
}

TEST_F(SolverMatrixTest, getEigenvalues3dComplex) {
    SolverMatrix m({ {33, -23, 9}, {22, 33, -23}, {19, 14, 50} });
    m /= 29;
    const Array expected({ {2} });
    const auto actual = m.getEigenvalues();
    expectEqual(expected, actual);
}

TEST_F(SolverMatrixTest, getEigenvectorsFor3dReal) {
    SolverMatrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(3, actual.rows());
    EXPECT_EQ(1, actual.columns());
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
    SolverMatrix m({ {1, 0, 0}, {0, 0, 0}, {0, 0, 1} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(2, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(contains(actual, 1.0));
    EXPECT_TRUE(contains(actual, 0.0));

    const auto actualVectors1 = m.getEigenvectorFor(1.0);
    EXPECT_EQ(2, actualVectors1.columns());
    const auto expected1 = Matrix({ { 1 }, { 0 }, { 0 } });
    const auto expected2 = Matrix({ { 0 }, { 0 }, { 1 } });
    expectNormalizedEqual(expected1, actualVectors1.getColumn(0), "actualVector1-1");
    expectNormalizedEqual(expected2, actualVectors1.getColumn(1), "actualVector1-2"); 

    const auto actualVector2 = m.getEigenvectorFor(0.0);
    expectNormalizedEqual(Matrix({ { 0 }, { 1 }, { 0 } }), actualVector2, "actualVector2");
}

TEST_F(SolverMatrixTest, getEigenVectorsTwoFreeVariables) {
    SolverMatrix m({ {1, 0, 0}, {0, 0, 0}, {0, 0, 1} });
    const auto actual = m.getEigenvectors();
    EXPECT_EQ(3, actual.rows());
    EXPECT_EQ(3, actual.columns());

    const auto expected = Matrix({ {0, 1, 0}, {1, 0, 0}, {0, 0, 1} }).transpose();
    for (int i = 0; i < 3; ++i) {
        expectNormalizedEqual(expected.getColumn(i), actual.getColumn(i), "column " + std::to_string(i));
    }
}

TEST_F(SolverMatrixTest, getEigenvectors3dReal) {
    SolverMatrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
    const auto actual = m.getEigenvectors();
    EXPECT_EQ(3, actual.rows());
    EXPECT_EQ(3, actual.columns());

    EXPECT_TRUE(normalizedContains(Matrix({ { -2 }, { 3 }, { 1 } }), actual)) << "eigenvector 1";
    EXPECT_TRUE(normalizedContains(Matrix({ { -2 }, { -1 }, { 1 } }), actual)) << "eigenvector 2";
    EXPECT_TRUE(normalizedContains(Matrix({ { 1 }, { 6 }, { 16 } }), actual)) << "eigenvector 3";
}
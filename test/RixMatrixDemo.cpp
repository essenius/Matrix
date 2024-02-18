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

#include <SolverMatrix.h>

using namespace RixMatrix;

void printArray(const char* caption, const Array& array) {
    printf("%s:\n", caption);
    for (Dimension row = 0; row < array.rowCount(); row++) {
        printf("  ");
        for (Dimension column = 0; column < array.columnCount(); column++) {
            printf("%.3f ", array(row, column));
        }
        printf("\n");
    }
}

void setup() {
#ifdef ESP32
    Serial.begin(115200);
    delay(500);
#endif
    printf("** Array operations **\n\n");

    const Array a({ {1, 2}, {3, 4} });
    printArray("Array a", a);
    // 1.000 2.000 
    // 3.000 4.000 
    printArray("a squared", a.pow2());
    // 1.000 4.000
    // 9.000 16.000
    printArray("a * 3", a * 3);
    // 3.000 6.000
    // 9.000 12.000
    printArray("a / 2", a / 2);
    // 0.500 1.000
    // 1.500 2.000

    Array b{ {0, 2}, {3, 0} };
    printArray("Array b", b);
    // 0.000 2.000 
    // 3.000 0.000 
    printArray("a + b", a + b);
    // 1.000 4.000
    // 6.000 4.000
    printArray("a - b", a - b);
    // 1.000 0.000
    // 0.000 4.000  
    printArray("a * b", a * b);
    // 0.000 4.000
    // 9.000 0.000 

    const Array c(b.getColumn(1));
    printArray("Array c (= second column of b)", c);
    // 2.000
    // 0.000 
    b.setColumn(1, Array{ {7}, {11} });
    printArray("b with new second column (7, 11)", b);
    // 0.000 7.000
    // 3.000 11.000  
    b.setRow(1, c.transposed());
    printArray("b with new second row c. transposed()", b);
    // 0.000 7.000
    // 2.000 0.000

    printf("\n** Matrix operations **\n\n");

    const Matrix m(a);
    printArray("Matrix m (= a)", m);
    // 1.000 2.000
    // 3.000 4.000
    printf("Determinant of m:\n  %.3f\n", m.getDeterminant());
    // -2.000
    printArray("m squared", m.squared());
    // 7.000 10.000 
    // 15.000 22.000 
    printArray("inverse of m", m.inverted());
    // -2.000 1.000
    // 1.500 -0.500
    printArray("Normalized m", m.normalized());
    // 0.183 0.365
    // 0.548 0.730
    printArray("Adjugate of m", m.getAdjugate());
    // 4.000 -2.000
    // -3.000 1.000

    const Matrix n({ {1, 2, 3}, {3, -1, 0} });
    printArray("Matrix n", n);  // 1.000 2.000 
    // 1.000 2.000 3.000 
    // 3.000 -1.000 0.000 

    const Matrix o({ {1, 2}, {3, 4}, {5, 6} });
    printArray("Matrix o", o);
    // 3.000 4.000 
    // 5.000 6.000 

    printArray("n * o", n * o);
    // 22.000 28.000
    // 0.000 2.000
    printArray("o * n", o * n);
    // 7.000 0.000 3.000
    // 15.000 2.000 9.000
    // 23.000 4.000 15.000

    const Matrix mtm = m.transposed<Matrix>() * m;
    printArray("mT * m", mtm);
    // 10.000 14.000 
    // 14.000 20.000   

    printf("\n** SolverMatrix operations **\n\n");

    const SolverMatrix solver({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
    printArray("Matrix to be solved", solver);
    // -2.000 -4.000 2.000 
    // -2.000 1.000 2.000 
    // 4.000 2.000 5.000
    printArray("Eigenvalues", solver.getEigenvalues());
    // 6.000
    // -5.000
    // 3.000

    printArray("Eigenvectors", solver.getEigenvectors());
    // 0.058 0.816 0.535 
    // 0.351 0.408 -0.802 
    // 0.935 -0.408 -0.267 
}

void loop() {
    // wait forever
}
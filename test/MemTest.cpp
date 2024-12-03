// Copyright 2024 Rik Essenius
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
#include "Matrix.h"
#include <Windows.h>
#include <psapi.h>

#include "SolverMatrix.h"

namespace RixMatrixTest {
    using RixMatrix::Matrix;

    long long getMemoryUsage() {
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
            return pmc.WorkingSetSize;
        }
        return 0;
    }

    TEST(MemTest, test1) {
        constexpr int tests = 200000;
        long long startMem = getMemoryUsage();
        for (int i=0; i<tests;i++)
        {
            const RixMatrix::SolverMatrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
            const auto actual = m.getEigenvectors();
        }
        long long endMem = getMemoryUsage();

        std::cout << "Start: " << startMem << " End: " << endMem << std::endl;
        long long difference = std::llabs(startMem - endMem);
        std::cout << "Difference: " << difference << std::endl;
        EXPECT_TRUE(false) << "End";

    }
}

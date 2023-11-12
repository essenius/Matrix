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

#ifndef ARRAYTEST_H
#define ARRAYTEST_H

#include "gtest/gtest.h"
#include "../src/Array.h"

class ArrayTest : public testing::Test {
protected:
    bool contains(const Array &matrix, double value, double epsilon = Array::Epsilon) const;
    bool isEqual(const Array &expected, const Array &actual, double epsilon = Array::Epsilon) const;
    void expectEqual(const Array &expected, const Array &actual, const std::string &message = "", double epsilon = Array::Epsilon) const;
};
#endif
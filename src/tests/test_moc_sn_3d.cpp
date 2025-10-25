/*
   Copyright 2024

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "UnitTest++/UnitTest++.h"

#include "driver.hpp"

namespace {
const std::vector<std::string> base_args = {"integration_test"};

void run_case(const std::string &input)
{
    auto args = base_args;
    args.push_back(input);
    run(args);
}
}

TEST(MoC3DCase)
{
    run_case("3d_moc.xml");
}

TEST(Sn3DCase)
{
    run_case("3d_sn.xml");
}

int main()
{
    return UnitTest::RunAllTests();
}

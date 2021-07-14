#define CATCH_CONFIG_MAIN CommonToolsTest
#include "catch.hpp"

#include "common_tools.hpp"

TEST_CASE("bin_vec_test")
{
    SECTION("constructor")
    {
        bin_vec vec(12);
        REQUIRE(vec.size() == 12);
    }

    SECTION("insert() and at()")
    {
        bin_vec vec(3);

        REQUIRE_THROWS_AS(vec.insert(3, 1), std::out_of_range);
        vec.insert(0, 0);
        vec.insert(2, 0);

        REQUIRE_THROWS_AS(vec.at(3), std::out_of_range);
        REQUIRE(vec.at(0) == 0);
        REQUIRE(vec.at(1) == 1);
        REQUIRE(vec.at(2) == 0);

        vec.insert(2, 1);
        REQUIRE(vec.at(0) == 0);
        REQUIRE(vec.at(1) == 1);
        REQUIRE(vec.at(2) == 1);
    }

    SECTION("and_()")
    {
        bin_vec vec1(3);

        vec1.insert(0, 0);
        vec1.insert(1, 1);
        vec1.insert(2, 1);

        bin_vec vec2(3);
        vec2.insert(0, 1);
        vec2.insert(1, 1);
        vec2.insert(2, 0);

        auto result = bin_vec::and_(vec1, vec2);
        REQUIRE(result == std::vector<bool>{false, true, false});
    }
}

TEST_CASE("mean_var_cov_test")
{
    SECTION("double mean(std::vector<value> X)")
    {
        std::vector<int> vec_i = {479, 720, 63, 435, 234, 967, 608, 1, 270, 991,
                                  220, 279, 836, 117, 498, 300, 101, 269, 983, 575,
                                  923, 401, 37, 197, 398, 355, 228, 865, 307, 897};
        auto res = mean(vec_i);
        REQUIRE(res == 451.8);

        std::vector<double> vec_d = {0.971, 0.564, 0.329, 0.962, 0.803, 0.894, 0.178, 0.709, 0.633, 0.607,
                                     0.917, 0.810, 0.410, 0.952, 0.234, 0.708, 0.168, 0.648, 0.066, 0.902,
                                     0.616, 0.426, 0.839, 0.397, 0.682, 0.398, 0.940, 0.190, 0.852, 0.822};
        res = mean(vec_d);
        REQUIRE(res == Approx(0.6209).margin(0.0001));
    }

    SECTION("double var(std::vector<value> X, double meanX)")
    {
        std::vector<int> vec_i = {479, 720, 63, 435, 234, 967, 608, 1, 270, 991,
                                  220, 279, 836, 117, 498, 300, 101, 269, 983, 575,
                                  923, 401, 37, 197, 398, 355, 228, 865, 307, 897};
        auto meanX = mean(vec_i);
        auto res = var(vec_i, meanX);
        REQUIRE(res == Approx(94151.293333).margin(0.000001));

        std::vector<double> vec_d = {0.971, 0.564, 0.329, 0.962, 0.803, 0.894, 0.178, 0.709, 0.633, 0.607,
                                     0.917, 0.810, 0.410, 0.952, 0.234, 0.708, 0.168, 0.648, 0.066, 0.902,
                                     0.616, 0.426, 0.839, 0.397, 0.682, 0.398, 0.940, 0.190, 0.852, 0.822};
        meanX = mean(vec_d);
        res = var(vec_d, meanX);
        REQUIRE(res == 0.07458149);
    }

    SECTION("double cov_X_Y(std::vector<value1> X_vec, double meanX, std::vector<value2> Y_vec, double meanY)")
    {
        std::vector<int> vec_i = {479, 720, 63, 435, 234, 967, 608, 1, 270, 991,
                                  220, 279, 836, 117, 498, 300, 101, 269, 983, 575,
                                  923, 401, 37, 197, 398, 355, 228, 865, 307, 897};

        std::vector<double> vec_d = {0.971, 0.564, 0.329, 0.962, 0.803, 0.894, 0.178, 0.709, 0.633, 0.607,
                                     0.917, 0.810, 0.410, 0.952, 0.234, 0.708, 0.168, 0.648, 0.066, 0.902,
                                     0.616, 0.426, 0.839, 0.397, 0.682, 0.398, 0.940, 0.190, 0.852, 0.822};
        auto meanX = mean(vec_i);
        auto meanY = mean(vec_d);
        auto res = cov_X_Y(vec_i, meanX, vec_d, meanY);
        REQUIRE(res == Approx(-18.231953).margin(0.000001));
    }
}

TEST_CASE("regr_r2_test")
{
    SECTION("std::array<double, 3> linear_regres_X_Y(std::vector<std::array<value, 2>> X_Y_vec)")
    {
        std::vector<double> vec_i = {479, 720, 63, 435, 234, 967, 608, 1, 270, 991,
                                    220, 279, 836, 117, 498, 300, 101, 269, 983, 575,
                                    923, 401, 37, 197, 398, 355, 228, 865, 307, 897};

        std::vector<double> vec_d = {0.971, 0.564, 0.329, 0.962, 0.803, 0.894, 0.178, 0.709, 0.633, 0.607,
                                    0.917, 0.810, 0.410, 0.952, 0.234, 0.708, 0.168, 0.648, 0.066, 0.902,
                                    0.616, 0.426, 0.839, 0.397, 0.682, 0.398, 0.940, 0.190, 0.852, 0.822};

        std::vector<std::array<double, 2>> vec_i_d(vec_i.size());
        auto vec_i_itr = vec_i.begin();
        auto vec_d_itr = vec_d.begin();

        for (auto &val : vec_i_d)
        {
            val = {*vec_i_itr++, *vec_d_itr++};
        }

        auto res = linear_regres_X_Y(vec_i_d, false);
        REQUIRE(res.at(0) == Approx(-0.0001936).margin(0.0000001));
        REQUIRE(res.at(1) == Approx(0.7083889).margin(0.0000001));
        REQUIRE(res.at(2) == Approx(0.0473379).margin(0.0000001));
    }
}

TEST_CASE("combination_test")
{
    REQUIRE(combination(2, 6) == 15);
}
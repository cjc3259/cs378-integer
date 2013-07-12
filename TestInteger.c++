// --------------------------------
// projects/integer/TestInteger.c++
// Copyright (C) 2013
// Glenn P. Downing
// --------------------------------

/*
To test the program:
    % ls /usr/include/cppunit/
    ...
    TestFixture.h
    ...
    % locate libcppunit.a
    /usr/lib/libcppunit.a
    % g++ -pedantic -std=c++0x -Wall Integer.h TestInteger.c++ -o TestInteger -lcppunit -ldl
    % valgrind TestInteger > TestInteger.out
*/

// --------
// includes
// --------

#include <algorithm> // equal
#include <cstring>   // strcmp
#include <sstream>   // ostringstream
#include <stdexcept> // invalid_argument
#include <string>    // ==
#include <iostream>

#include "cppunit/extensions/HelperMacros.h" // CPPUNIT_TEST, CPPUNIT_TEST_SUITE, CPPUNIT_TEST_SUITE_END
#include "cppunit/TestFixture.h"             // TestFixture
#include "cppunit/TextTestRunner.h"          // TextTestRunner

#include "Integer.h"

using namespace std;
// -----------
// TestInteger
// -----------

struct TestInteger : CppUnit::TestFixture {
    // -----------------
    // shift_left_digits
    // -----------------

    void test_shift_left_digits () {
        const int a[] = {2, 3, 4};
        const int b[] = {2, 3, 4, 0, 0};
              int x[10];
        const int* p = shift_left_digits(a, a + 3, 2, x);
        CPPUNIT_ASSERT((p - x) == 5);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, b));}

    // ------------------
    // shift_right_digits
    // ------------------

    void test_shift_right_digits () {
        const int a[] = {2, 3, 4};
        const int b[] = {2};
              int x[10];
        const int* p = shift_right_digits(a, a + 3, 2, x);
        // std::cout << p - x;
        CPPUNIT_ASSERT((p - x) == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, b));}

    void test_shift_right_digits_2 () {
        const int a[] = {2, 3, 4};
        const int b[] = {0};
              int x[10];
        const int* p = shift_right_digits(a, a + 3, 3, x);
        // std::cout << p - x;
        CPPUNIT_ASSERT((p - x) == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, b));}

    void test_shift_right_digits_3 () {
        const int a[] = {2, 3, 4};
        const int b[] = {0};
              int x[10];
        const int* p = shift_right_digits(a, a + 3, 4, x);
        // std::cout << p - x;
        CPPUNIT_ASSERT((p - x) == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, b));}

    // -----------
    // plus_digits
    // -----------

    void test_plus_digits () {
        const int a[] = {2, 3, 4};
        const int b[] = {5, 6, 7};
        const int c[] = {8, 0, 1};
              int x[10];
        const int* p = plus_digits(a, a + 3, b, b + 3, x);
        CPPUNIT_ASSERT(p - x == 3);
        // cout << x[0] << x[1] << x[2] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_plus_digits_2 () {
        const int a[] = {2, 3, 4};
        const int b[] = {6, 7};
        const int c[] = {3, 0, 1};
              int x[10];
        const int* p = plus_digits(a, a + 3, b, b + 2, x);
        CPPUNIT_ASSERT(p - x == 3);
        // cout << x[0] << x[1] << x[2] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_plus_digits_3 () {
        const int a[] = {6, 7};
        const int b[] = {2, 3, 4};
        const int c[] = {3, 0, 1};
              int x[10];
        const int* p = plus_digits(a, a + 2, b, b + 3, x);
        // cout << p - x << endl;
        CPPUNIT_ASSERT(p - x == 3);
        // cout << x[0] << x[1] << x[2] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_plus_digits_4 () {
        const int a[] = {1};
        const int b[] = {9, 9, 9};
        const int c[] = {1, 0, 0, 0};
              int x[10];
        const int* p = plus_digits(a, a + 1, b, b + 3, x);
        // cout << p - x << endl;
        CPPUNIT_ASSERT(p - x == 4);
        // cout << x[0] << x[1] << x[2] << x[3] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_plus_digits_5 () {
        const int a[] = {9, 9, 9};
        const int b[] = {9, 9, 9};
        const int c[] = {1, 9, 9, 8};
              int x[10];
        const int* p = plus_digits(a, a + 3, b, b + 3, x);
        // cout << p - x << endl;
        CPPUNIT_ASSERT(p - x == 4);
        // cout << x[0] << x[1] << x[2] << x[3] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_plus_digits_6 () {
        const int a[] = {1, 6, 3, 8};
        const int b[] = {1, 4, 0, 4, 0};
        const int c[] = {1, 5, 6, 7, 8};
              int x[10];
        const int* p = plus_digits(a, a + 4, b, b + 5, x);
        // cout << p - x << endl;
        CPPUNIT_ASSERT(p - x == 5);
        // cout << x[0] << x[1] << x[2] << x[3] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_plus_digits_7 () {
        const int a[] = {9, 0, 0};
        const int b[] = {1};
        const int c[] = {9, 0, 1};
              int x[10];
        const int* p = plus_digits(a, a + 3, b, b + 1, x);
        // cout << p - x << endl;
        CPPUNIT_ASSERT(p - x == 3);
        // cout << x[0] << x[1] << x[2] << x[3] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    // ------------
    // minus_digits
    // ------------

    void test_minus_digits () {
        const int a[] = {8, 0, 1};
        const int b[] = {5, 6, 7};
        const int c[] = {2, 3, 4};
              int x[10];
        const int* p = minus_digits(a, a + 3, b, b + 3, x);
        CPPUNIT_ASSERT(p - x == 3);
        // cout << x[0] << x[1] << x[2] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_minus_digits_2 () {
        const int a[] = {1, 0, 1};
        const int b[] = {1};
        const int c[] = {1, 0, 0};
              int x[10];
        const int* p = minus_digits(a, a + 3, b, b + 1, x);
        // cout << p - x << endl;
        // cout << x[0] << x[1] << x[2] << x[3] << endl;
        CPPUNIT_ASSERT(p - x == 3);
        // cout << x[0] << x[1] << x[2] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_minus_digits_3 () {
        const int a[] = {1, 0, 0};
        const int b[] = {1};
        const int c[] = {9, 9};
              int x[10];
        const int* p = minus_digits(a, a + 3, b, b + 1, x);
        // cout << p - x << endl;
        // cout << x[0] << x[1] << x[2] << x[3] << endl;
        CPPUNIT_ASSERT(p - x == 2);
        // cout << x[0] << x[1] << x[2] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_minus_digits_4 () {
        const int a[] = {1, 0, 0};
        const int b[] = {9, 9};
        const int c[] = {1};
              int x[10];
        const int* p = minus_digits(a, a + 3, b, b + 2, x);
        // cout << p - x << endl;
        // cout << x[0] << x[1] << x[2] << x[3] << endl;
        CPPUNIT_ASSERT(p - x == 1);
        // cout << x[0] << x[1] << x[2] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_minus_digits_5 () {
        const int a[] = {1};
        const int b[] = {1};
        const int c[] = {0};
              int x[10];
        const int* p = minus_digits(a, a + 1, b, b + 1, x);
        // cout << p - x << endl;
        // cout << x[0] << endl;
        CPPUNIT_ASSERT(p - x == 1);
        // cout << x[0] << x[1] << x[2] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_minus_digits_6 () {
        const int a[] = {1, 3, 2, 1, 1, 1};
        const int b[] = {5, 6, 7};
        const int c[] = {1, 3, 1, 5, 4, 4};
              int x[10];
        const int* p = minus_digits(a, a + 6, b, b + 3, x);
        // cout << p - x << endl;
        // cout << x[0] << endl;
        CPPUNIT_ASSERT(p - x == 6);
        cout << x[0] << x[1] << x[2] << x[3] << x[4] << x[5] << endl;
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    // -----------------
    // multiplies_digits
    // -----------------

    void test_multiplies_digits () {
        const int a[] = {2, 3, 4};
        const int b[] = {5, 6, 7};
        const int c[] = {1, 3, 2, 6, 7, 8};
              int x[10];
        const int* p = multiplies_digits(a, a + 3, b, b + 3, x);
        CPPUNIT_ASSERT(p - x == 6);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_2 () {
        const int a[] = {2, 3, 4};
        const int b[] = {7};
        const int c[] = {1, 6, 3, 8};
              int x[10];
        const int* p = multiplies_digits(a, a + 3, b, b + 1, x);
        CPPUNIT_ASSERT(p - x == 4);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_3 () {
        const int a[] = {7};
        const int b[] = {2, 3, 4};
        const int c[] = {1, 6, 3, 8};
              int x[10];
        const int* p = multiplies_digits(a, a + 1, b, b + 3, x);
        CPPUNIT_ASSERT(p - x == 4);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_4 () {
        const int a[] = {2, 3, 4};
        const int b[] = {0};
        const int c[] = {0};
              int x[10];
        const int* p = multiplies_digits(a, a + 3, b, b + 1, x);
        CPPUNIT_ASSERT(p - x == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_5 () {
        const int a[] = {0};
        const int b[] = {2, 3, 4};
        const int c[] = {0};
              int x[10];
        const int* p = multiplies_digits(a, a + 1, b, b + 3, x);
        CPPUNIT_ASSERT(p - x == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_6 () {
        const int a[] = {9, 0};
        const int b[] = {1};
        const int c[] = {9, 0};
              int x[10];
        const int* p = multiplies_digits(a, a + 2, b, b + 1, x);
        CPPUNIT_ASSERT(p - x == 2);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_7 () {
        const int a[] = {8, 5};
        const int b[] = {2};
        const int c[] = {1, 7, 0};
              int x[10];
        const int* p = multiplies_digits(a, a + 2, b, b + 1, x);
        CPPUNIT_ASSERT(p - x == 3);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_8 () {
        const int a[] = {8, 0};
        const int b[] = {1};
        const int c[] = {8, 0};
              int x[10];
        const int* p = multiplies_digits(a, a + 2, b, b + 1, x);
        CPPUNIT_ASSERT(p - x == 2);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_9 () {
        const int a[] = {1, 5};
        const int b[] = {2};
        const int c[] = {3, 0};
              int x[10];
        const int* p = multiplies_digits(a, a + 2, b, b + 1, x);
        cout << "p - x" << p - x << endl;
        CPPUNIT_ASSERT(p - x == 2);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_multiplies_digits_10 () {
        const int a[] = {1, 3, 4, 2, 1, 7, 7, 2, 8};
        const int b[] = {2};
        const int c[] = {2, 6, 8, 4, 3, 5, 4, 5, 6};
              int x[20];
        const int* p = multiplies_digits(a, a + 9, b, b + 1, x);
        // cout << "p - x" << p - x << endl;
        CPPUNIT_ASSERT(p - x == 9);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    // --------------
    // divides_digits
    // --------------

    void test_divides_digits () {
        const int a[] = {1, 3, 2, 6, 7, 8};
        const int b[] = {5, 6, 7};
        const int c[] = {2, 3, 4};
              int x[10];
        const int* p = divides_digits(a, a + 6, b, b + 3, x);
        cout << "Test x size = " << (p - x) << endl;
        CPPUNIT_ASSERT(p - x == 3);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_divides_digits_2 () {
        const int a[] = {1, 3};
        const int b[] = {5, 6, 7};
        const int c[] = {0};
              int x[10];
        const int* p = divides_digits(a, a + 2, b, b + 3, x);
        CPPUNIT_ASSERT(p - x == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_divides_digits_3 () {
        const int a[] = {9};
        const int b[] = {1};
        const int c[] = {9};
              int x[10];
        const int* p = divides_digits(a, a + 1, b, b + 1, x);
        CPPUNIT_ASSERT(p - x == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_divides_digits_4 () {
        const int a[] = {9};
        const int b[] = {2};
        const int c[] = {4};
              int x[10];
        const int* p = divides_digits(a, a + 1, b, b + 1, x);
        CPPUNIT_ASSERT(p - x == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_divides_digits_5 () {
        const int a[] = {0};
        const int b[] = {2};
        const int c[] = {0};
              int x[10];
        const int* p = divides_digits(a, a + 1, b, b + 1, x);
        CPPUNIT_ASSERT(p - x == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}

    void test_divides_digits_6 () {
        const int a[] = {4, 5, 3, 4, 5};
        const int b[] = {5, 3, 4, 5};
        const int c[] = {8};
              int x[10];
        const int* p = divides_digits(a, a + 5, b, b + 4, x);
        CPPUNIT_ASSERT(p - x == 1);
        CPPUNIT_ASSERT(std::equal(const_cast<const int*>(x), p, c));}


    // -----------
    // constructor
    // -----------

    void test_constructor_1 () {
        try {
            const Integer<int> x("abc");
            CPPUNIT_ASSERT(false);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(strcmp(e.what(), "Integer()"));}}

    void test_constructor_2 () {
        try {
            const Integer<int> x("-abc");
            CPPUNIT_ASSERT(false);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(strcmp(e.what(), "Integer()"));}}

    void test_constructor_3 () {
        try {
            const Integer<int> x("2");}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    void test_constructor_4 () {
        try {
            const Integer<int> x("-2");}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    void test_constructor_5 () {
        try {
            const Integer<int> x(2);}
            // cout << x << endl;
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    void test_constructor_6 () {
        try {
            const Integer<int> x(-442);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    void test_constructor_7 () {
        try {
            const Integer<int> x("-1c");
            CPPUNIT_ASSERT(false);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(strcmp(e.what(), "Integer()"));}}

    void test_constructor_8 () {
        try {
            const Integer<int> x("");
            CPPUNIT_ASSERT(false);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(strcmp(e.what(), "Integer()"));}}

    void test_constructor_9 () {
        try {
            const Integer<int> x(0);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    // ---
    // abs
    // ---

    void test_abs_1 () {
        try {
            Integer<int>  x = -98765;
            Integer<int>& y = x.abs();
            CPPUNIT_ASSERT( x == 98765);
            CPPUNIT_ASSERT(&x == &y);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    void test_abs_2 () {
        try {
            const Integer<int> x = -98765;
            const Integer<int> y = abs(x);
            CPPUNIT_ASSERT(x == -98765);
            CPPUNIT_ASSERT(y ==  98765);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    // --------
    // equal_to
    // --------

    void test_equal_to () {
        try {
            const Integer<int> x = 98765;
            const Integer<int> y = 98765;
            CPPUNIT_ASSERT(      x == y);
            CPPUNIT_ASSERT(      x == 98765);
            CPPUNIT_ASSERT(  98765 == x);
            CPPUNIT_ASSERT(    !(x != y));
            CPPUNIT_ASSERT(    !(x != 98765));
            CPPUNIT_ASSERT(!(98765 != y));}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    // --------
    // negation
    // --------

    void test_negation () {
        try {
            const Integer<int> x = 98765;
            const Integer<int> y = -x;
            cout << y << endl;
            CPPUNIT_ASSERT(y == -98765);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    // ------
    // output
    // ------

    void test_output () {
        try {
            const Integer<int> x = 98765;
            std::ostringstream out;
            out << x;
            CPPUNIT_ASSERT(out.str() == "98765");}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    // ---
    // pow
    // ---

    void test_pow_1 () {
        try {
            cout << "test_pow" << endl;
            Integer<int>       x = 10;
            const int          e =  7;
            Integer<int>&      y = x.pow(e);
            CPPUNIT_ASSERT( e == 7);
            CPPUNIT_ASSERT( x ==    10000000);
            CPPUNIT_ASSERT(&x ==   &y);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    void test_pow_2 () {
        try {
            cout << "test_pow" << endl;
            Integer<int>       x = 2;
            const int          e =  19;
            Integer<int>&      y = x.pow(e);
            CPPUNIT_ASSERT( e == 19);
            CPPUNIT_ASSERT( x ==    524288);
            CPPUNIT_ASSERT(&x ==   &y);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    void test_pow_3 () {
        try {
            Integer<int>       x = 9;
            const int          e =  2;
            Integer<int>&      y = x.pow(e);
            CPPUNIT_ASSERT( e == 2);
            CPPUNIT_ASSERT( x ==    81);
            CPPUNIT_ASSERT(&x ==   &y);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}

    void test_pow_4 () {
        try {
            Integer<int>       x = 9;
            const int          e =  3;
            Integer<int>&      y = x.pow(e);
            CPPUNIT_ASSERT( e == 3);
            CPPUNIT_ASSERT( x ==    729);
            CPPUNIT_ASSERT(&x ==   &y);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}}



    void test_Integer_plus() {
        try {
            const Integer<int> x = 900;
            const Integer<int> y = 1;
            const Integer<int> z = x + y;
            // cout << z << endl;
            CPPUNIT_ASSERT(x == 900);
            CPPUNIT_ASSERT(y == 1);
            CPPUNIT_ASSERT(z == 901);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_plus_2() {
        try {
            const Integer<int> x = 85;
            const Integer<int> y = 900;
            const Integer<int> z = x + y;
            // cout << z << endl;
            CPPUNIT_ASSERT(x == 85);
            CPPUNIT_ASSERT(y == 900);
            CPPUNIT_ASSERT(z == 985);}
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_multiply() {
        try {
            const Integer<int> x = 85;
            const Integer<int> y = 2;
            const Integer<int> z = x * y;
            CPPUNIT_ASSERT(x == 85);
            CPPUNIT_ASSERT(y == 2);
            CPPUNIT_ASSERT(z == 170);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_multiply_2() {
        try {
            const Integer<int> x = 234;
            const Integer<int> y = 567;
            const Integer<int> z = x * y;
            CPPUNIT_ASSERT(x == 234);
            CPPUNIT_ASSERT(y == 567);
            CPPUNIT_ASSERT(z == 132678);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_multiply_3() {
        try {
            const Integer<int> x = 85;
            const Integer<int> y = 2;
            const Integer<int> z = x * y;
            CPPUNIT_ASSERT(x == 85);
            CPPUNIT_ASSERT(y == 2);
            CPPUNIT_ASSERT(z == 170);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_multiply_4() {
        try {
            const Integer<int> x = 16;
            const Integer<int> y = 2;
            const Integer<int> z = x * y;
            CPPUNIT_ASSERT(x == 16);
            CPPUNIT_ASSERT(y == 2);
            CPPUNIT_ASSERT(z == 32);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_minus() {
        try {
            const Integer<int> x = 16;
            const Integer<int> y = 2;
            const Integer<int> z = x - y;
            CPPUNIT_ASSERT(x == 16);
            CPPUNIT_ASSERT(y == 2);
            CPPUNIT_ASSERT(z == 14);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_minus_2() {
        try {
            const Integer<int> x = 16;
            const Integer<int> y = 1;
            const Integer<int> z = x - y;
            CPPUNIT_ASSERT(x == 16);
            CPPUNIT_ASSERT(y == 1);
            CPPUNIT_ASSERT(z == 15);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_minus_3() {
        try {
            const Integer<int> x = 16;
            const Integer<int> y = 0;
            const Integer<int> z = x - y;
            CPPUNIT_ASSERT(x == 16);
            CPPUNIT_ASSERT(y == 0);
            CPPUNIT_ASSERT(z == 16);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_divide() {
        try {
            const Integer<int> x = 16;
            const Integer<int> y = 2;
            const Integer<int> z = x/y;
            cout << "z " << z << endl;
            CPPUNIT_ASSERT(x == 16);
            CPPUNIT_ASSERT(y == 2);
            CPPUNIT_ASSERT(z == 8);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_divide_2() {
        try {
            const Integer<int> x = 16;
            const Integer<int> y = 5;
            const Integer<int> z = x/y;
            cout << "z " << z << endl;
            CPPUNIT_ASSERT(x == 16);
            CPPUNIT_ASSERT(y == 5);
            CPPUNIT_ASSERT(z == 3);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_divide_3() {
        try {
            const Integer<int> x = 44999;
            const Integer<int> y = 5000;
            const Integer<int> z = x/y;
            cout << "z " << z << endl;
            CPPUNIT_ASSERT(x == 44999);
            CPPUNIT_ASSERT(y == 5000);
            CPPUNIT_ASSERT(z == 8);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_divide_4() {
        try {
            const Integer<int> x = 45345;
            const Integer<int> y = 5345;
            const Integer<int> z = x/y;
            cout << "z " << z << endl;
            CPPUNIT_ASSERT(x == 45345);
            CPPUNIT_ASSERT(y == 5345);
            CPPUNIT_ASSERT(z == 8);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }


    void test_Integer_shift_left() {
        try {
            Integer<int> x = 45345;
            x <<= 2;
            cout << "x " << x << endl;
            CPPUNIT_ASSERT(x == 4534500);
            // CPPUNIT_ASSERT(y == 5345);
            // CPPUNIT_ASSERT(z == 8);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }

    void test_Integer_shift_right() {
        try {
            Integer<int> x = 45345;
            x >>= 2;
            cout << "x " << x << endl;
            CPPUNIT_ASSERT(x == 453);
            // CPPUNIT_ASSERT(y == 5345);
            // CPPUNIT_ASSERT(z == 8);
        }
        catch (std::invalid_argument& e) {
            CPPUNIT_ASSERT(false);}
    }
    // -----
    // suite
    // -----

    CPPUNIT_TEST_SUITE(TestInteger);
    CPPUNIT_TEST(test_shift_left_digits);
    CPPUNIT_TEST(test_shift_right_digits);
    CPPUNIT_TEST(test_shift_right_digits_2);
    CPPUNIT_TEST(test_shift_right_digits_3);
    CPPUNIT_TEST(test_plus_digits);
    CPPUNIT_TEST(test_plus_digits_2);
    CPPUNIT_TEST(test_plus_digits_3);
    CPPUNIT_TEST(test_plus_digits_4);
    CPPUNIT_TEST(test_plus_digits_5);
    CPPUNIT_TEST(test_plus_digits_6);
    CPPUNIT_TEST(test_plus_digits_7);
    CPPUNIT_TEST(test_minus_digits);
    CPPUNIT_TEST(test_minus_digits_2);
    CPPUNIT_TEST(test_minus_digits_3);
    CPPUNIT_TEST(test_minus_digits_4);
    CPPUNIT_TEST(test_minus_digits_5);
    CPPUNIT_TEST(test_minus_digits_6);

    CPPUNIT_TEST(test_multiplies_digits);
    CPPUNIT_TEST(test_multiplies_digits_2);
    CPPUNIT_TEST(test_multiplies_digits_3);
    CPPUNIT_TEST(test_multiplies_digits_4);
    CPPUNIT_TEST(test_multiplies_digits_5);
    CPPUNIT_TEST(test_multiplies_digits_6);
    CPPUNIT_TEST(test_multiplies_digits_7);
    CPPUNIT_TEST(test_multiplies_digits_8);
    CPPUNIT_TEST(test_multiplies_digits_9);
    CPPUNIT_TEST(test_multiplies_digits_10);


    CPPUNIT_TEST(test_divides_digits);
    CPPUNIT_TEST(test_divides_digits_2);
    CPPUNIT_TEST(test_divides_digits_3);
    CPPUNIT_TEST(test_divides_digits_4);
    CPPUNIT_TEST(test_divides_digits_5);
    CPPUNIT_TEST(test_divides_digits_6);



    CPPUNIT_TEST(test_constructor_1);
    CPPUNIT_TEST(test_constructor_2);
    CPPUNIT_TEST(test_constructor_3);
    CPPUNIT_TEST(test_constructor_4);
    CPPUNIT_TEST(test_constructor_5);
    CPPUNIT_TEST(test_constructor_6);
    CPPUNIT_TEST(test_constructor_7);
    CPPUNIT_TEST(test_constructor_8);
    CPPUNIT_TEST(test_constructor_9);

    // CPPUNIT_TEST(test_abs_1);
    // CPPUNIT_TEST(test_abs_2);
    CPPUNIT_TEST(test_equal_to);
    // CPPUNIT_TEST(test_negation);
    CPPUNIT_TEST(test_output);
    CPPUNIT_TEST(test_pow_1);
    CPPUNIT_TEST(test_pow_2);
    CPPUNIT_TEST(test_pow_3);
    CPPUNIT_TEST(test_pow_4);

    CPPUNIT_TEST(test_Integer_plus);
    CPPUNIT_TEST(test_Integer_plus_2);
    CPPUNIT_TEST(test_Integer_multiply);
    CPPUNIT_TEST(test_Integer_multiply_2);
    CPPUNIT_TEST(test_Integer_multiply_3);
    CPPUNIT_TEST(test_Integer_multiply_4);
    CPPUNIT_TEST(test_Integer_minus);
    CPPUNIT_TEST(test_Integer_minus_2);
    CPPUNIT_TEST(test_Integer_minus_3);
    CPPUNIT_TEST(test_Integer_divide);
    CPPUNIT_TEST(test_Integer_divide_2);
    CPPUNIT_TEST(test_Integer_divide_3);
    CPPUNIT_TEST(test_Integer_divide_4);
    CPPUNIT_TEST(test_Integer_shift_left);
    CPPUNIT_TEST(test_Integer_shift_right);









    CPPUNIT_TEST_SUITE_END();};

// ----
// main
// ----

int main () {
    using namespace std;
    cout << "TestInteger.c++" << endl;

    CppUnit::TextTestRunner tr;
    tr.addTest(TestInteger::suite());
    tr.run();

    cout << "Done." << endl;
    return 0;}

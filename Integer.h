// --------------------------
// projects/integer/Integer.h
// Copyright (C) 2013
// Glenn P. Downing
// --------------------------

#ifndef Integer_h
#define Integer_h

// --------
// includes
// --------

#include <algorithm> // reverse
#include <cassert>   // assert
#include <iostream>  // ostream
#include <iterator>  // not sure if we need this or not..
#include <stdexcept> // invalid_argument
#include <string>    // string
#include <vector>    // vector

using namespace std;

// -----------------
// shift_left_digits
// -----------------

/**
 * @param b an iterator to the beginning of an input  sequence (inclusive)
 * @param e an iterator to the end       of an input  sequence (exclusive)
 * @param x an iterator to the beginning of an output sequence (inclusive)
 * @return  an iterator to the end       of an output sequence (exclusive)
 * the sequences are of decimal digits
 * output the shift left of the input sequence into the output sequence
 * ([b, e) << n) => x
 */
template <typename II, typename OI>
OI shift_left_digits (II b, II e, int n, OI x) {
    // <your code>
    while(b != e) {
        *x = *b;
        x++;
        b++;
    }
    for(int i = 0; i < n; i++) {
        *x = 0;
        x++;
    }
    return x;}

// ------------------
// shift_right_digits
// ------------------

/**
 * @param b an iterator to the beginning of an input  sequence (inclusive)
 * @param e an iterator to the end       of an input  sequence (exclusive)
 * @param x an iterator to the beginning of an output sequence (inclusive)
 * @return  an iterator to the end       of an output sequence (exclusive)
 * the sequences are of decimal digits
 * output the shift right of the input sequence into the output sequence
 * ([b, e) >> n) => x
 */
template <typename II, typename OI>
OI shift_right_digits (II b, II e, int n, OI x) {
    // <your code>
    if (n >= (e - b)) {
        *x = 0;
        x++;
    }
    else {
        while (b != (e - n)) {
            *x = *b;
            x++;
            b++;
        }
    }
    return x;}

// -----------
// plus_digits
// -----------

/**
 * @param b  an iterator to the beginning of an input  sequence (inclusive)
 * @param e  an iterator to the end       of an input  sequence (exclusive)
 * @param b2 an iterator to the beginning of an input  sequence (inclusive)
 * @param e2 an iterator to the end       of an input  sequence (exclusive)
 * @param x  an iterator to the beginning of an output sequence (inclusive)
 * @return   an iterator to the end       of an output sequence (exclusive)
 * the sequences are of decimal digits
 * output the sum of the two input sequences into the output sequence
 * ([b1, e1) + [b2, e2)) => x
 */
template <typename II1, typename II2, typename OI>
OI plus_digits (II1 b1, II1 e1, II2 b2, II2 e2, OI x) {
    // <your code>
    int sum = 0;
    int carry = 0;
    int count = 0;
    int s1 = distance(b1, e1);
    int s2 = distance(b2, e2);

    // if first number is bigger or equal to second number
    if (s1 >= s2) {
        for (int i = 0; i < s2; i++) {
            e1--;
            e2--;
            sum = carry + *e1 + *e2;
            if (sum > 9) {            
                sum = sum % 10;
                carry = 1;
            } 
            else carry = 0;
            *x = sum;
            sum = 0;
            count++;
            x++;
            
        }
        while (b1 != e1) {
            sum = carry + *(e1 - 1);
            if (sum > 9) {
                sum = sum % 10;
                carry = 1;
            } 
            else carry = 0;
            *x = sum;
            sum = 0;
            e1--;
            count++;
            x++;
        }
    }

    // if second number is bigger than first number
    else {
        for (int i = 0; i < s1; i++) {
            e1--;
            e2--;
            sum = carry + *e1 + *e2;
            if (sum > 9) {
                sum = sum%10;
                carry = 1;
            } 
            else carry = 0;
            *x = sum;
            sum = 0;
            count++;
            x++;
            
        }
        while (b2 != e2) {
            sum = carry + *(e2 - 1);
            if (sum > 9) {
                sum = sum%10;
                carry = 1;
            } 
            else carry = 0;
            *x = sum;
            sum = 0;
            e2--;
            count++;
            x++;
        }
    }

    // if there is still a carry after adding all 'columns', add the 1
    if (carry == 1) {
        *x = 1;
        count++;
        x++;
    }

    // reverse the order of the numbers for output
    vector<int> y(x - count, x);
    reverse(y.begin(), y.end());
    x = x - count;
    for(int i = 0; i < y.size(); i++){
        *x = y[i];
        x++;
    }
    return x;}

// ------------
// minus_digits
// ------------

/**
 * @param b  an iterator to the beginning of an input  sequence (inclusive)
 * @param e  an iterator to the end       of an input  sequence (exclusive)
 * @param b2 an iterator to the beginning of an input  sequence (inclusive)
 * @param e2 an iterator to the end       of an input  sequence (exclusive)
 * @param x  an iterator to the beginning of an output sequence (inclusive)
 * @return   an iterator to the end       of an output sequence (exclusive)
 * the sequences are of decimal digits
 * output the difference of the two input sequences into the output sequence
 * ([b1, e1) - [b2, e2)) => x
 */
template <typename II1, typename II2, typename OI>
OI minus_digits (II1 b1, II1 e1, II2 b2, II2 e2, OI x) {
    // <your code>
    int s2 = distance(b2, e2);
    int count = 0;
    bool borrow = false;

    // perform minus for all columns in smallest number 
    for (int i = 0; i < s2; i++){
        e1--;
        e2--;
        int minuend = *e1;
        if (borrow) {
            if (*e1 == 0) {
                minuend = 9;
            }
            else {
                minuend = *e1 - 1;
                borrow = false;
            }
        }
        if (minuend < *e2) {
            minuend = *e1 + 10;
            borrow = true;
        }
        *x = minuend - *e2;
        x++;
        count++;
    }

    // continue subtracting for the rest of the minuend
    while (b1 != e1) {

        int minuend = *(e1 - 1);
        if (borrow) {
            if (minuend == 0) {
                minuend = 9;
            }
            else {
                minuend = minuend - 1;
                borrow = false;
            }
        }
        *x = minuend;
        e1--;
        x++;
        count++;
    }

    // clearing leading zeros (e.g.--0001 = 1).
    vector<int> y(x - count, x);
    for(int i = y.size() - 1; i > 0; --i) {
        if (y[i] == 0) {
            y.pop_back();
        }
        else break;
    }

    // reverse order and output
    reverse(y.begin(), y.end());
    x = x - count;
    for(int i = 0; i < y.size(); i++){
        *x = y[i];
        x++;
    }
    return x;}

// -----------------
// multiplies_digits
// -----------------

/**
 * @param b  an iterator to the beginning of an input  sequence (inclusive)
 * @param e  an iterator to the end       of an input  sequence (exclusive)
 * @param b2 an iterator to the beginning of an input  sequence (inclusive)
 * @param e2 an iterator to the end       of an input  sequence (exclusive)
 * @param x  an iterator to the beginning of an output sequence (inclusive)
 * @return   an iterator to the end       of an output sequence (exclusive)
 * the sequences are of decimal digits
 * output the product of the two input sequences into the output sequence
 * ([b1, e1) * [b2, e2)) => x
 */
template <typename II1, typename II2, typename OI>
OI multiplies_digits (II1 b1, II1 e1, II2 b2, II2 e2, OI x) {
    // <your code>
    int s1 = distance(b1, e1);
    int s2 = distance(b2, e2);
    vector<int> m1;
    vector<int> m2;
    vector<vector<int> > w;

    // if either number is 0, output 0
    if((*b1 == 0 && (e1 - b1) == 1) || (*b2 == 0 && (e2 - b2) == 1)) {
        *x = 0;
        return x + 1;
    }
    
    // order the two numbers, 'biggest' first
    if (s1 >= s2) {
        for(int i = 0; i < s1; ++i) {
            m1.push_back(*b1 + i);
        }
        for(int i = 0; i < s2; ++i) {
            m2.push_back(*b2 + i);
        }
    }
    else {
        for(int i = 0; i < s1; ++i) {
            m2.push_back(*b1 + i);
        }
        for(int i = 0; i < s2; ++i) {
            m1.push_back(*b2 + i);
        }
    }

    // perform multiplication on digits to get matrix of numbers to eventually add up for the product
    for(int i = 0; i < m2.size(); ++i) {
        int carry = 0;
        vector<int> d;
        w.push_back(d);
        for(int j = 0; j < i; ++j) {
            w[i].push_back(0);
        }
        for(int k = 0; k < m1.size(); ++k) {
            int p = m2[m2.size() - 1 - i]*m1[m1.size() - 1 - k] + carry;
            if(p > 9) {
                carry = p/10;
                p = p%10;
            }
            w[i].push_back(p);
        }
        if (carry != 0) {
            w[i] .push_back(carry);
        }
    }

    // add up multiple vectors of numbers to get the final product
    vector<int> product;
    product.push_back(0);
    *x = 0;
    int* product_end = x + 1;
    for (int i = 0; i < w.size(); ++i) {
        reverse(w[i].begin(), w[i].end());

        product_end = plus_digits(product.begin(), product.end(), w[i].begin(), w[i].end(), x);

        product.clear();
        for (int k = 0; k < (product_end - x); ++k) {
            product.push_back(*(x + k));
        }
    }
    return product_end;}

// -------------- 
// divides_digits
// --------------

/**
 * @param b  an iterator to the beginning of an input  sequence (inclusive)
 * @param e  an iterator to the end       of an input  sequence (exclusive)
 * @param b2 an iterator to the beginning of an input  sequence (inclusive)
 * @param e2 an iterator to the end       of an input  sequence (exclusive)
 * @param x  an iterator to the beginning of an output sequence (inclusive)
 * @return   an iterator to the end       of an output sequence (exclusive)
 * the sequences are of decimal digits
 * output the division of the two input sequences into the output sequence
 * ([b1, e1) / [b2, e2)) => x
 */
template <typename II1, typename II2, typename OI>
OI divides_digits (II1 b1, II1 e1, II2 b2, II2 e2, OI x) {
    // <your code>
    return x;}

// -------
// Integer
// -------

template < typename T, typename C = std::vector<T> >
class Integer {
    // -----------
    // operator ==
    // -----------

    /**
     * <your documentation>
     */
    friend bool operator == (const Integer& lhs, const Integer& rhs) {
        // <your code>
        return false;}

    // -----------
    // operator !=
    // -----------

    /**
     * <your documentation>
     */
    friend bool operator != (const Integer& lhs, const Integer& rhs) {
        return !(lhs == rhs);}

    // ----------
    // operator <
    // ----------

    /**
     * <your documentation>
     */
    friend bool operator < (const Integer& lhs, const Integer& rhs) {
        // <your code>
        return false;}

    // -----------
    // operator <=
    // -----------

    /**
     * <your documentation>
     */
    friend bool operator <= (const Integer& lhs, const Integer& rhs) {
        return !(rhs < lhs);}

    // ----------
    // operator >
    // ----------

    /**
     * <your documentation>
     */
    friend bool operator > (const Integer& lhs, const Integer& rhs) {
        return (rhs < lhs);}

    // -----------
    // operator >=
    // -----------

    /**
     * <your documentation>
     */
    friend bool operator >= (const Integer& lhs, const Integer& rhs) {
        return !(lhs < rhs);}

    // ----------
    // operator +
    // ----------

    /**
     * <your documentation>
     */
    friend Integer operator + (Integer lhs, const Integer& rhs) {
        return lhs += rhs;}

    // ----------
    // operator -
    // ----------

    /**
     * <your documentation>
     */
    friend Integer operator - (Integer lhs, const Integer& rhs) {
        return lhs -= rhs;}

    // ----------
    // operator *
    // ----------

    /**
     * <your documentation>
     */
    friend Integer operator * (Integer lhs, const Integer& rhs) {
        return lhs *= rhs;}

    // ----------
    // operator /
    // ----------

    /**
     * <your documentation>
     * @throws invalid_argument if (rhs == 0)
     */
    friend Integer operator / (Integer lhs, const Integer& rhs) {
        return lhs /= rhs;}

    // ----------
    // operator %
    // ----------

    /**
     * <your documentation>
     * @throws invalid_argument if (rhs <= 0)
     */
    friend Integer operator % (Integer lhs, const Integer& rhs) {
        return lhs %= rhs;}

    // -----------
    // operator <<
    // -----------

    /**
     * <your documentation>
     * @throws invalid_argument if (rhs < 0)
     */
    friend Integer operator << (Integer lhs, int rhs) {
        return lhs <<= rhs;}

    // -----------
    // operator >>
    // -----------

    /**
     * <your documentation>
     * @throws invalid_argument if (rhs < 0)
     */
    friend Integer operator >> (Integer lhs, int rhs) {
        return lhs >>= rhs;}

    // -----------
    // operator <<
    // -----------

    /**
     * <your documentation>
     */
    friend std::ostream& operator << (std::ostream& lhs, const Integer& rhs) {
        // <your code>
        return lhs << "0";}

    // ---
    // abs
    // ---

    /**
     * absolute value
     * does NOT modify the argument
     * <your documentation>
     */
    friend Integer abs (Integer x) {
        return x.abs();}

    // ---
    // pow
    // ---

    /**
     * power
     * does NOT modify the argument
     * <your documentation>
     * @throws invalid_argument if (x == 0) && (e == 0)
     * @throws invalid_argument if (e < 0)
     */
    friend Integer pow (Integer x, int e) {
        return x.pow(e);}

    private:
        // ----
        // data
        // ----

        // <your data>

    private:
        // -----
        // valid
        // -----

        bool valid () const {
            // <your code>
            return true;}

    public:
        // ------------
        // constructors
        // ------------

        /**
         * <your documentation>
         */
        Integer (int value) {
            // <your code>
            assert(valid());}

        /**
         * <your documentation>
         * @throws invalid_argument if value is not a valid representation of an Integer
         */
        explicit Integer (const std::string& value) {
            // <your code>
            if (!valid())
                throw std::invalid_argument("Integer::Integer()");}

        // Default copy, destructor, and copy assignment.
        // Integer (const Integer&);
        // ~Integer ();
        // Integer& operator = (const Integer&);

        // ----------
        // operator -
        // ----------

        /**
         * <your documentation>
         */
        Integer operator - () const {
            // <your code>
            return Integer(0);}

        // -----------
        // operator ++
        // -----------

        /**
         * <your documentation>
         */
        Integer& operator ++ () {
            *this += 1;
            return *this;}

        /**
         * <your documentation>
         */
        Integer operator ++ (int) {
            Integer x = *this;
            ++(*this);
            return x;}

        // -----------
        // operator --
        // -----------

        /**
         * <your documentation>
         */
        Integer& operator -- () {
            *this -= 1;
            return *this;}

        /**
         * <your documentation>
         */
        Integer operator -- (int) {
            Integer x = *this;
            --(*this);
            return x;}

        // -----------
        // operator +=
        // -----------

        /**
         * <your documentation>
         */
        Integer& operator += (const Integer& rhs) {
            // <your code>
            return *this;}

        // -----------
        // operator -=
        // -----------

        /**
         * <your documentation>
         */
        Integer& operator -= (const Integer& rhs) {
            // <your code>
            return *this;}

        // -----------
        // operator *=
        // -----------

        /**
         * <your documentation>
         */
        Integer& operator *= (const Integer& rhs) {
            // <your code>
            return *this;}

        // -----------
        // operator /=
        // -----------

        /**
         * <your documentation>
         * @throws invalid_argument if (rhs == 0)
         */
        Integer& operator /= (const Integer& rhs) {
            // <your code>
            return *this;}

        // -----------
        // operator %=
        // -----------

        /**
         * <your documentation>
         * @throws invalid_argument if (rhs <= 0)
         */
        Integer& operator %= (const Integer& rhs) {
            // <your code>
            return *this;}

        // ------------
        // operator <<=
        // ------------

        /**
         * <your documentation>
         */
        Integer& operator <<= (int n) {
            // <your code>
            return *this;}

        // ------------
        // operator >>=
        // ------------

        /**
         * <your documentation>
         */
        Integer& operator >>= (int n) {
            // <your code>
            return *this;}

        // ---
        // abs
        // ---

        /**
         * absolute value
         * <your documentation>
         */
        Integer& abs () {
            // <your code>
            return *this;}

        // ---
        // pow
        // ---

        /**
         * power
         * <your documentation>
         * @throws invalid_argument if (this == 0) && (e == 0)
         * @throws invalid_argument if (e < 0)
         */
        Integer& pow (int e) {
            // <your code>
            return *this;}};

#endif // Integer_h

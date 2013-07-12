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
#include <stdlib.h>  // atoi, abs
#include <stdio.h>   // isdigit
#include <string>    // string
#include <vector>    // vector
#include <deque>     // deque
#include <typeinfo>  // check data type
#include <ctype.h>   // isdigit
#include <cstddef>   // find_first_of

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
    

    // y.assign(x - count, x);
    // reverse(y.begin(), y.end());
    
    // reverse(y.begin(), y.end());

    vector<int> y(x - count, x);
    x = x - count;
    for(int i = y.size() - 1; i >= 0; --i){
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
            minuend = minuend + 10;
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
    // reverse(y.begin(), y.end());
    x = x - count;
    for(int i = y.size() - 1; i >= 0; --i){
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
    // m1.resize(40000);
    vector<int> m2;
    vector<vector<int> > w;
    // w[0].resize(40000);

    for (int i = 0; i < (e1 - b1); ++i) {
    }
    // if either number is 0, output 0
    if((*b1 == 0 && (e1 - b1) == 1) || (*b2 == 0 && (e2 - b2) == 1)) {
        *x = 0;
        return x + 1;
    }
    
    // order the two numbers, 'biggest' first
    if (s1 >= s2) {
        for(int i = 0; i < s1; ++i) {
            m1.push_back(*(b1 + i));
        }
    
        for(int i = 0; i < s2; ++i) {
            m2.push_back(*(b2 + i));
        }
    }
    else {
        for(int i = 0; i < s1; ++i) {
            m2.push_back(*(b1 + i));
        }
        for(int i = 0; i < s2; ++i) {
            m1.push_back(*(b2 + i));
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
            carry = 0;
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

	int s1 = distance(b1, e1);
    int s2 = distance(b2, e2);
    *x = 0;
    // x++;
    if(s1 < s2)
        return x + 1;

    if(s1 == 1 && *b1 == 0) 
        return x + 1;

    bool less_than = false;
    vector<int> divisor;
    // int* divisor_end = &divisor[0] + 1;
    vector<int> dummy;
    int* dummy_end = &dummy[0] + 1;
    dummy.push_back(0);

    int* x_end = x + 1;

    vector<int> increment;
    increment.push_back(1);

    vector<int> count;
    count.push_back(0);

    for(int i = 0; i < s1; ++i) {
        divisor.push_back(*(b1 + i));
    }

    for(int i = 0; i < s2; ++i) {
    }
    
    while (!less_than) {
        dummy_end = minus_digits(divisor.begin(), divisor.end(), b2, e2, &dummy[0]);

        divisor.clear();
        for (int i = 0; i < (dummy_end - &dummy[0]); ++i) {
            divisor.push_back(dummy[i]);
        }



        x_end = plus_digits(increment.begin(), increment.end(), count.begin(), count.end(), x);
        
        count.clear();
        for(int i = 0; i < (x_end - x); ++i) {
            count.push_back(*(x + i));
        }

        s1 = distance(divisor.begin(), divisor.end());
        
        if(s1 < s2)
            less_than = true;
        if(s1 == s2) {
            for(int i = 0; i < s1; ++i) {
                if(*(divisor.begin() + i) < *(b2 + i)) {
                    less_than = true;
                    break;
                }
                else if (*(divisor.begin() + i) > *(b2 + i))
                    break;
            }
        }

    }
    return x_end;}

// -------
// Integer
// -------

template < typename T, typename C = std::vector<T> >
class Integer  {

    // -----------
    // operator ==
    // -----------

    /**
     *  Check for equivalence between the lhs and rhs
     */
    friend bool operator == (const Integer& lhs, const Integer& rhs) {
        // <your code>
        return lhs.v == rhs.v;
    }

    // -----------
    // operator !=
    // -----------

    /**
     * != works
     */
    friend bool operator != (const Integer& lhs, const Integer& rhs) {
        return !(lhs == rhs);}

    // ----------
    // operator <
    // ----------

    /**
     *  Operator < is a comparison between the two values
     */
    friend bool operator < (const Integer& lhs, const Integer& rhs) {
        // <your code>
        return lhs.v < rhs.v;}

    // -----------
    // operator <=
    // -----------

    /**
     * <= works
     */
    friend bool operator <= (const Integer& lhs, const Integer& rhs) {
        return !(rhs < lhs);}

    // ----------
    // operator >
    // ----------

    /**
     * > works
     */
    friend bool operator > (const Integer& lhs, const Integer& rhs) {
        return (rhs < lhs);}

    // -----------
    // operator >=
    // -----------

    /**
     * >= works
     */
    friend bool operator >= (const Integer& lhs, const Integer& rhs) {
        return !(lhs < rhs);}

    // ----------
    // operator +
    // ----------

    /**
     * addition works
     */
    friend Integer operator + (Integer lhs, const Integer& rhs) {
        return lhs += rhs;}

    // ----------
    // operator -
    // ----------

    /**
     * subtraction works
     */
    friend Integer operator - (Integer lhs, const Integer& rhs) {
        return lhs -= rhs;}

    // ----------
    // operator *
    // ----------

    /**
     * multiplication works
     */
    friend Integer operator * (Integer lhs, const Integer& rhs) {
        return lhs *= rhs;}

    // ----------
    // operator /
    // ----------

    /**
     * division works
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
     * output stream works
     */
    friend std::ostream& operator << (std::ostream& lhs, const Integer& rhs) {
        // <your code>
        // o_interator<int> out_it
        Integer<T, C> x = rhs;
        for (int i = 0; i < x.size(); ++i) {
            lhs << x.value(i);
        }
        return lhs;}

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
        C v;

    private:
        // -----
        // valid
        // -----

        bool valid () const {
            // <your code>

            if(v.empty())
                return false;
            return true;
        }

    public:

        
        // ------------
        // constructors
        // ------------

        /**
         * int constructor
         */
        Integer (int value) {
            // <your code>
            if (value == 0)
                v.push_back(0);
            while(value != 0) {
                int d = value%10;
                if(!(value < 10 && value > -10)){
                    if(d < 0)
                        v.push_back(-d);
                    else v.push_back(d);
                } 
                else v.push_back(d);
                value/=10;
            }
            reverse(v.begin(), v.end());
            assert(valid());}

        /**
         * string constructor
         * @throws invalid_argument if value is not a valid representation of an Integer
         */
        explicit Integer (const std::string& value) {
            // <your code>
            

            bool isInteger = true;
            for (int i = 0; i < value.size(); ++i) {

                if(value.find_first_of('-') == i)
                    continue;
                if(!(isdigit(value[i]))) {
                    isInteger = false;
                    break;   
                }
            }
            if (isInteger) {
                int i = atoi(value.c_str());
                while(i != 0) {
                    int d = i%10;
                    if(!(i < 10 && i > -10)){
                        if(d < 0)
                            v.push_back(-d);
                        else v.push_back(d);
                    } 
                    else v.push_back(d);
                    i/=10;
                }
                reverse(v.begin(),v.end());
                
            }
            
            if (!valid())
                throw std::invalid_argument("Integer::Integer()");
        }

        // typedef Integer<T, C>::iterator iterator;
        // typedef Integer<T, C>::const_iterator const_iterator;


        // iterator begin() {
        //     // cout << typeid(v).name();
        //     return v.begin();
        // }
        // iterator end() {
        //     return v.end();
        // }

        int* array_iter(int i) {
            return &v[i];
        }

        // void push_back(int i) {
        //     v.push_back(i);
        // }

        // void pop_back() {
        //     v.pop_back();
        // }

        // void clear() {
        //     v.clear();
        // }

        // void resize(int n) {
        //     v.resize(n);
        // }


        // Default copy, destructor, and copy assignment.
        // Integer (const Integer&);
        // ~Integer ();
        // Integer& operator = (const Integer&);

        long size() {
            return v.size();
        }

        int value(int i) {
            return v[i];
        }

        void set(int i, int val) {
            v[i] = val;
        }

        // template<class II>
        // void assign (II first, II last){
        //     v.assign(first, last);
        // }

        // void assign (size_type n, const value_type& val) {
        //     v.assign(n, val);
        // }


        // ----------
        // operator -
        // ----------

        /**
         * <your documentation>
         */
        Integer operator - () const {
            // <your code>

            // cout << *((*this).begin()) << endl;
            return *this;}

        // -----------
        // operator ++
        // -----------

        /**                *this

         * <your documentation>
         */
        Integer& operator ++ () {
            *this += 1;
            return *this;}

        /**
         * <your documentation>
         */
        Integer operator ++ (int) {
            Integer<T, C> x = *this;
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
            Integer<T, C> x = *this;
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
            Integer<T, C> x = *this;
            Integer<T, C> y = rhs;
            Integer<T, C> z = Integer<T, C>(0);
            z.v.resize(1332);
            int* begin = (z).array_iter(0);
            int* end = (z).array_iter(0) + 1;
            end = plus_digits(x.v.begin(), x.v.begin() + x.v.size(), y.v.begin(), y.v.begin() + y.v.size(), begin);
            (*this).v.clear();
            (*this).v.resize(1332);
            (*this).v.assign(begin, end);
            return *this;}

        // -----------
        // operator -=
        // -----------

        /**
         * <your documentation>
         */
        Integer& operator -= (const Integer& rhs) {
            // <your code>
            Integer<T, C> x = *this;
            Integer<T, C> y = rhs;
            Integer<T, C> z = Integer<T, C>(0);
            z.v.resize((*this).size());
            int* begin = (z).array_iter(0);
            int* end = (z).array_iter(0) + 1;
            end = minus_digits(x.v.begin(), x.v.begin() + x.v.size(), y.v.begin(), y.v.begin() + y.v.size(), begin);
            (*this).v.clear();
            (*this).v.resize((*this).size());
            (*this).v.assign(begin, end);
            return *this;}

        // -----------
        // operator *=
        // -----------

        /**
         * <your documentation>
         */
        Integer& operator *= (const Integer& rhs) {
            // <your code>
   
            Integer<T, C> x = *this;
            Integer<T, C> y = rhs;
            Integer<T, C> z = Integer<T, C>(0);
            // cout << "+ " << endl;
            z.v.resize((*this).size() + 1);
            int* begin = (z).array_iter(0);
            int* end = (z).array_iter(0) + 1;
            end = multiplies_digits(x.v.begin(), x.v.begin() + x.v.size(), y.v.begin(), y.v.begin() + y.v.size(), begin);
            (*this).v.clear();
            (*this).v.resize((*this).size() + 1);
            (*this).v.assign(begin, end);
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
            Integer<T, C> x = *this;
            Integer<T, C> y = rhs;
            Integer<T, C> z = Integer<T, C>(0);
            z.v.resize((*this).size());
            int* begin = (z).array_iter(0);
            int* end = (z).array_iter(0) + 1;
            end = divides_digits(x.v.begin(), x.v.begin() + x.v.size(), y.v.begin(), y.v.begin() + y.v.size(), begin);
            cout << "divide " << endl;
            for(int i = 0; i < end - begin; ++i) {
                cout << *(begin + 1) ;
            }
            cout << endl;
            (*this).v.clear();
            (*this).v.resize((*this).size());
            (*this).v.assign(begin, end);
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
            Integer<T, C> x = *this;
            Integer<T, C> z = Integer<T, C>(0);
            z.v.resize((*this).size() + n);
            int* begin = (z).array_iter(0);
            int* end = (z).array_iter(0) + 1;
            end = shift_left_digits(x.v.begin(), x.v.begin() + x.v.size(), n, begin);
            // for(int i = 0; i < end - begin; ++i) {
            //     cout << *(begin + 1) ;
            // }
            // cout << endl;
            (*this).v.clear();
            (*this).v.resize((*this).size() + 1);
            (*this).v.assign(begin, end);
            return *this;}

        // ------------
        // operator >>=
        // ------------

        /**
         * <your documentation>
         */
        Integer& operator >>= (int n) {
            // <your code>
            Integer<T, C> x = *this;
            Integer<T, C> z = Integer<T, C>(0);
            z.v.resize((*this).size());
            int* begin = (z).array_iter(0);
            int* end = (z).array_iter(0) + 1;
            end = shift_right_digits(x.v.begin(), x.v.begin() + x.v.size(), n, begin);
            // for(int i = 0; i < end - begin; ++i) {
            //     cout << *(begin + 1) ;
            // }
            // cout << endl;
            (*this).v.clear();
            (*this).v.resize((*this).size());
            (*this).v.assign(begin, end);
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
            // if 
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
            // cout << "+ " << endl;
            Integer<T, C> x = *this;
            // Integer<T, C> y = *this;
            // cout << "+ " << endl;
            for(int i = 1; i < e; ++i) {
                *this *= x;
                // cout << (*this).v.size() << endl;
                // cout << i + 1 << endl;
            }
            
            return *this;
        }
    };

#endif // Integer_h

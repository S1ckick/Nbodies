//
// Created by Максим on 07.01.2021.
//

#ifndef NBODIES_HELPER_H
#define NBODIES_HELPER_H

#ifdef NUMBER_DOUBLE_DOUBLE
#include <qd/dd_real.h>
#include <string>
#include <iostream>
using namespace std;

void round_string(char *s, int precision, int *offset){
    /*
     Input string must be all digits or errors will occur.
     */

    int i;
    int D = precision ;

    /* Round, handle carry */
    if (D>0 && s[D] >= '5') {
        s[D-1]++;

        i = D-1;
        while (i > 0 && s[i] > '9') {
            s[i] -= 10;
            s[--i]++;
        }
    }

    /* If first digit is 10, shift everything. */
    if (s[0] > '9') {
        // e++; // don't modify exponent here
        for (i = precision; i >= 1; i--) s[i+1] = s[i];
        s[0] = '1';
        s[1] = '0';

        (*offset)++ ; // now offset needs to be increased by one
        precision++ ;
    }

    s[precision] = 0; // add terminator for array
}

void append_expn(std::string &str, int expn) {
    int k;

    str += (expn < 0 ? '-' : '+');
    expn = std::abs(expn);

    if (expn >= 100) {
        k = (expn / 100);
        str += '0' + k;
        expn -= 100*k;
    }

    k = (expn / 10);
    str += '0' + k;
    expn -= 10*k;

    str += '0' + expn;
}

string dd_real::to_string(int precision, int width, ios_base::fmtflags fmt,
                          bool showpos, bool uppercase, char fill) const {
    string s;
    bool fixed = (fmt & ios_base::fixed) != 0;
    bool sgn = true;
    int i, e = 0;

    if (isnan()) {
        s = uppercase ? "NAN" : "nan";
        sgn = false;
    } else {
        if (*this < 0.0)
            s += '-';
        else if (showpos)
            s += '+';
        else
            sgn = false;

        if (isinf()) {
            s += uppercase ? "INF" : "inf";
        } else if (*this == 0.0) {
            /* Zero case */
            s += '0';
            if (precision > 0) {
                s += '.';
                s.append(precision, '0');
            }
        } else {
            /* Non-zero case */
            int off = (fixed ? (1 + to_int(floor(log10(abs(*this))))) : 1);
            int d = precision + off;

            int d_with_extra = d;
            if(fixed)
                d_with_extra = std::max(60, d); // longer than the max accuracy for DD

            // highly special case - fixed mode, precision is zero, abs(*this) < 1.0
            // without this trap a number like 0.9 printed fixed with 0 precision prints as 0
            // should be rounded to 1.
            if(fixed && (precision == 0) && (abs(*this) < 1.0)){
                if(abs(*this) >= 0.5)
                    s += '1';
                else
                    s += '0';

                return s;
            }

            // handle near zero to working precision (but not exactly zero)
            if (fixed && d <= 0) {
                s += '0';
                if (precision > 0) {
                    s += '.';
                    s.append(precision, '0');
                }
            } else { // default

                char *t; //  = new char[d+1];
                int j;

                if(fixed){
                    t = new char[d_with_extra+1];
                    to_digits(t, e, d_with_extra);
                }
                else{
                    t = new char[d+1];
                    to_digits(t, e, d);
                }

                off = e + 1;

                if (fixed) {
                    // fix the string if it's been computed incorrectly
                    // round here in the decimal string if required
                    round_string(t, d, &off);

                    if (off > 0) {
                        for (i = 0; i < off; i++) s += t[i];
                        if (precision > 0) {
                            s += '.';
                            for (j = 0; j < precision; j++, i++) s += t[i];
                        }
                    } else {
                        s += "0.";
                        if (off < 0) s.append(-off, '0');
                        for (i = 0; i < d; i++) s += t[i];
                    }
                } else {
                    s += t[0];
                    if (precision > 0) s += '.';

                    for (i = 1; i <= precision; i++)
                        s += t[i];

                }
                delete [] t;
            }
        }

        // trap for improper offset with large values
        // without this trap, output of values of the for 10^j - 1 fail for j > 28
        // and are output with the point in the wrong place, leading to a dramatically off value
        if(fixed && (precision > 0)){
            // make sure that the value isn't dramatically larger
            double from_string = atof(s.c_str());

            // if this ratio is large, then we've got problems
            if( fabs( from_string / this->x[0] ) > 3.0 ){

                int point_position;
                char temp;

                // loop on the string, find the point, move it up one
                // don't act on the first character
                for(i=1; i < s.length(); i++){
                    if(s[i] == '.'){
                        s[i] = s[i-1] ;
                        s[i-1] = '.' ;
                        break;
                    }
                }

                from_string = atof(s.c_str());
                // if this ratio is large, then the string has not been fixed
                if( fabs( from_string / this->x[0] ) > 3.0 ){
                    dd_real::error("Re-rounding unsuccessful in large number fixed point trap.") ;
                }
            }
        }


        if (!fixed && !isinf()) {
            /* Fill in exponent part */
            s += uppercase ? 'E' : 'e';
            append_expn(s, e);
        }
    }

    /* Fill in the blanks */
    int len = s.length();
    if (len < width) {
        int delta = width - len;
        if (fmt & ios_base::internal) {
            if (sgn)
                s.insert(static_cast<string::size_type>(1), delta, fill);
            else
                s.insert(static_cast<string::size_type>(0), delta, fill);
        } else if (fmt & ios_base::left) {
            s.append(delta, fill);
        } else {
            s.insert(static_cast<string::size_type>(0), delta, fill);
        }
    }

    return s;
}

/* Reads in a double-double number from the string s. */
void read(const char *s, dd_real &a) {
    const char *p = s;
    char ch;
    int sign = 0;
    int point = -1;
    int nd = 0;
    int e = 0;
    bool done = false;
    dd_real r = 0.0;
    int nread;

    /* Skip any leading spaces */
    while (*p == ' ')
        p++;

    while (!done && (ch = *p) != '\0') {
        if (ch >= '0' && ch <= '9') {
            int d = ch - '0';
            r *= 10.0;
            r += static_cast<double>(d);
            nd++;
        } else {

            switch (ch) {

                case '.':
                    if (point >= 0)
                        return;
                    point = nd;
                    break;

                case '-':
                case '+':
                    if (sign != 0 || nd > 0)
                        return;
                    sign = (ch == '-') ? -1 : 1;
                    break;

                case 'E':
                case 'e':
                    nread = std::sscanf(p+1, "%d", &e);
                    done = true;
                    if (nread != 1)
                        return;
                    break;

                default:
                    return;
            }
        }

        p++;
    }

    if (point >= 0) {
        e -= (nd - point);
    }

    if (e != 0) {
        r *= (dd_real(10.0) ^ e);
    }

    a = (sign == -1) ? -r : r;
    return;
}



#endif

#endif //NBODIES_HELPER_H

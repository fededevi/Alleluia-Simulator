#pragma once

#include <math.h>
#include <sys/types.h>
#include <vector>

double sum(const std::vector<double> & a) {
    double s = 0;
    for (uint i = 0; i < a.size(); i++) {
        s += a[i];
    }
    return s;
}

double mean(const std::vector<double> &  a) {
    return sum(a) / (double)a.size();
}

double sqsum(const std::vector<double> & a) {
    double s = 0;
    for (uint i = 0; i < a.size(); i++) {
        s += pow(a[i], 2);
    }
    return s;
}

double stdev(const std::vector<double> & nums) {
    double N = (double)nums.size();
    return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

std::vector<double> operator-(const std::vector<double> & a, double b) {
    std::vector<double> retvect;
    retvect.reserve(a.size());
    for (uint i = 0; i < a.size(); i++) {
        retvect.push_back(a[i] - b);
    }
    return retvect;
}

std::vector<double> operator*(const std::vector<double> & a, const std::vector<double> & b) {
    std::vector<double> retvect;
    retvect.reserve(a.size());
    for (uint i = 0; i < a.size() ; i++) {
        retvect.push_back(a[i] * b[i]);
    }
    return retvect;
}

double pearsoncoeff(const std::vector<double> & X, const std::vector<double> & Y) {
    return sum((X - mean(X))*(Y - mean(Y))) / ((double)X.size()*stdev(X)* stdev(Y));
}

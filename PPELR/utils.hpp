//
//  utils.h
//  PPELR
//
//  Created by Dave on 1/6/15.
//  Copyright (c) 2015 Todai. All rights reserved.
//

#ifndef __PPELR__utils__
#define __PPELR__utils__

#include <stdio.h>
#include <vector>
#include <numeric>

template <typename T, typename C>
std::vector<int> sort_permutation(std::vector<T> const& vec, C compare)
{
    std::vector<int> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&](int i, int j){ return compare(vec[i], vec[j]); });
    return p;
}


template <typename T>
std::vector<T> apply_permutation(std::vector<T> const& vec, std::vector<int> const& p)
{
//    std::vector<T> sorted_vec(p.size());
    std::vector<T> sorted_vec(vec); // workardound to be able to create a new vector of ciphers
    std::transform(p.begin(), p.end(), sorted_vec.begin(), [&](int i){ return vec[i]; });
    return sorted_vec;
}


#endif /* defined(__PPELR__utils__) */

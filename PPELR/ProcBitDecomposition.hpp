//
//  ProcBitDecomposition.h
//  PPELR
//
//  Created by Dave dV on 1/8/15.
//  Copyright (c) 2015 Todai. All rights reserved.
//

#ifndef __PPELR__ProcBitDecomposition__
#define __PPELR__ProcBitDecomposition__

#include <stdio.h>

#include "Procedure.hpp"
#include "Cryptosystem.hpp"

namespace ppf {
    
    template<class CryptosystemType>
    class ProcBitDecomposition : public Procedure<CryptosystemType> {
        
    public:
        using LocalRawType = typename CryptosystemType::Cipher::RawType;
        using LocalCipher = typename CryptosystemType::Cipher;
        using RemoteRawType = typename CryptosystemType::PublicKeyCSType::Cipher::RawType;
        using RemoteCipher = typename CryptosystemType::PublicKeyCSType::Cipher;
        
        static const char* procName() { return "BITS"; };
        
        void run(msgpack::rpc::request req) {
            auto outParams = ProcBitDecomposition::template unpack<std::vector<LocalCipher>, int>(req);
            
            std::vector<LocalCipher> innerProdRs = std::get<0>(outParams);
            int nBits = std::get<1>(outParams);
            
            //DEBUG: for verification purposes:
//            std::vector<LocalCipher> innerProds = std::get<2>(outParams);
//            std::vector<LocalCipher> Rs = std::get<3>(outParams);
            
            LocalCipher c(*this->localCS);
            
            SERVER_LOG("got vector of length: " << innerProdRs.size());
            SERVER_LOG("LSB >= " << nBits << " bits");
            
            std::vector<std::vector<LocalCipher>> allBits;
            LocalCipher one(*this->localCS), zero(*this->localCS);
            one = 1;
            zero = 0;

            gmp_randstate_t r_state;
            gmp_randinit_default (r_state);
            Plain r;

            for(int j = 0; j < innerProdRs.size(); j++) {
                int nBitsAdjusted = nBits;
                

                //TODO: compute the exact needed bitLength on the other side
                while(mpz_tstbit(innerProdRs[j].as_mpz_class().get_mpz_t(), nBitsAdjusted) == 0 && nBitsAdjusted < 128)
                    nBitsAdjusted++;
                nBitsAdjusted+=2;
                
//                SERVER_LOG("Using " << nBitsAdjusted << " bits");
                
                std::vector<LocalCipher> bits;
                for(int b = 0; b < nBitsAdjusted; b++) {
                    //Encode from most significant to least
                    if(mpz_tstbit(innerProdRs[j].as_mpz_class().get_mpz_t(), nBitsAdjusted-b-1) == 1) {
                        bits.push_back(one);
                        one.recrypt();
                    }
                    else {
                        bits.push_back(zero);
                        zero.recrypt();
                    }
                }
                
                allBits.push_back(bits);
                
            }
            
            //DEBUG:
//            std::vector<Plain> x(nBits + 64);
//            Plain p, pR;
//            Plain pOne = 1, two = 2, pMOne = -1;
//            
//            for (int i = 0; i < allBits.size(); i++) {
//                std::vector<LocalCipher> y = allBits[i];
//                
//                if((i+1) % 100 == 0)
//                    SERVER_LOG("Running DEBUG Comparison protocol for " << (i+1) << "th value");
//                
//                //Encode from most significant to least:
//                pR = Rs[i];
//                for(int b = 0; b < y.size(); b++) {
//                    p = mpz_tstbit(pR.as_mpz_class().get_mpz_t(), y.size() - b - 1);
//                    x[b] = p;
//                }
//                
//                std::vector<LocalCipher> d(y.size(), *this->localCS);
//                std::vector<LocalCipher> f(y.size(), *this->localCS);
//                std::vector<LocalCipher> gamma(y.size(), *this->localCS);
//                LocalCipher delta(*this->localCS);
//                LocalCipher c(*this->localCS);
//                gamma[0] = 0;
//                
//                bool found = false;
//                
//                for(int j = 0; j < y.size(); j++) {
//                    if(x[j] == 0) {
//                        d[j] = y[j];
//                        f[j] = y[j];
//                    }
//                    else {
//                        d[j] = y[j] - pOne;
//                        f[j] = y[j] * pMOne;
//                        f[j] += pOne;
//                    }
//                    
//                    if(j > 0) {
//                        gamma[j] = gamma[j-1]*two;
//                        gamma[j] += f[j];
//                    }
//                    //Todo: directly assign a random cipher
//                    mpz_urandomb(r.raw().get_mpz_t(), r_state, this->localCS->getBitLength());
//                    
//                    c = gamma[j];
//                    c -= 1;
//                    
//                    delta = c * r + d[j];
//                    
//                    if(delta == 1)
//                        found = true;
//                }
//                
//                if((innerProds[i] >= 0) != found) {
//                    std::cout << "ERROR" << std::endl;
//                    
//                    std::cout << innerProds[i] << " + " << Rs[i] << " = " << innerProdRs[i] << std::endl;
//                    std::cout << "x  |  y  |  d  f  g" << std::endl;
//
//                    for(int j = 0; j < y.size(); j++) {
//                        std::cout << x[j] << " | " << y[j] << "   " << d[j] << "   " << f[j] << "   " << gamma[j] << "  " << (j > 0 ? (gamma[j-1]*two) : -1) << std::endl;
//                    }
//                }
//            }

            
            ProcBitDecomposition::packResult(req, allBits);

        };
        
    };
    
}


#endif /* defined(__PPELR__ProcBitDecomposition__) */

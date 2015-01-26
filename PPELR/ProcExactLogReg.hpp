//
//  ProcExactLogReg.h
//  PPFramework
//
//  Created by Dave on 12/12/14.
//  Copyright (c) 2014 Tsudalab. All rights reserved.
//

#ifndef __PPFramework__ProcExactLogReg__
#define __PPFramework__ProcExactLogReg__

#include <stdio.h>

#pragma GCC diagnostic push
#pragma clang diagnostic ignored "-Wshorten-64-to-32" // Stupid Boost lib is not yet completely 64bit-friendly
#include <boost/thread/thread.hpp>
//#include <boost/thread/shared_mutex.hpp>
#pragma GCC diagnostic pop

#include "Procedure.hpp"
#include "Cryptosystem.hpp"
#include "ProcBitDecomposition.hpp"

namespace ppf {
    
    template<class CryptosystemType>
    class ProcExactLogReg : public Procedure<CryptosystemType> {

    public:
        using LocalRawType = typename CryptosystemType::Cipher::RawType;
        using LocalCipher = typename CryptosystemType::Cipher;
        using RemoteRawType = typename CryptosystemType::PublicKeyCSType::Cipher::RawType;
        using RemoteCipher = typename CryptosystemType::PublicKeyCSType::Cipher;

        static const char* procName() { return "ELR"; };
        
        ProcExactLogReg(std::vector<std::vector<bool>>& _x_1n) : x_1n(_x_1n) {
            boost::lock_guard<boost::shared_mutex> lockX_1(x_1Access);

            if(x_1n.size() == 0)
                PROC_ERROR("Need at least one x1 vector");
            size_t s = x_1n[0].size();
            if(s == 0)
                PROC_ERROR("x1 vectors need to have at least one value");
            for(auto x_1 : x_1n) {
                if(x_1.size() != s)
                    PROC_ERROR("all x1 vectors must have the same length");
            }
        }

        void computeT1s(std::vector<int> useIndices) {
            //Todo make sure that subsequent useIndices are subsets of initial one
            boost::shared_lock<boost::shared_mutex> lock(x_1Access);
            boost::lock_guard<boost::mutex> lockT_1(t_1Access);
            if(t_1n.size() == 0) {
                SERVER_LOG("Computing real t1 vector");

                t_1n.resize(x_1n.size(), **this->remoteCS);
                
                RemoteCipher t_1(**this->remoteCS);
                
                for(int idx : useIndices) {
                    t_1 = 0;
                    for(int j = 0; j < y.size(); j++) {
                        if(x_1n[idx][j])
                            t_1 += y[j];
                    }
                    t_1n[idx] = t_1;
                }
            }
        }
        
        void run(msgpack::rpc::request req) {
            auto outParams = ProcExactLogReg::template unpack<std::vector<RemoteCipher>, std::vector<std::vector<RemoteCipher>>, std::vector<int>>(req);
            

            std::vector<std::vector<RemoteCipher>> allResults;

            std::vector<RemoteCipher> inY = std::get<0>(outParams);
            std::vector<std::vector<RemoteCipher>> y_samples = std::get<1>(outParams);
            std::vector<int> snpIndices = std::get<2>(outParams);

            
            auto totInst = x_1n[0].size();
            
            boost::upgrade_lock< boost::shared_mutex> lock(yAccess);
            if(y.size() == 0 || inY.size() > 0) {
                if(inY.size() != totInst)
                    PROC_ERROR("real y vector must be of same length as x_1");
                boost::upgrade_to_unique_lock<boost::shared_mutex> uniqueLock(lock);
                y = inY;
            }
            
            if(y_samples.size() == 0)
                PROC_ERROR("need at least one y sample");
            for(auto y_sample : y_samples)
                if(y_sample.size() != totInst)
                    PROC_ERROR("sampled y vectors must be of same length as x_1");
            
            SERVER_LOG("Received " << y_samples.size() << " sampled vectors of length: " << totInst);
            
            RemoteCipher t_1(**this->remoteCS), cOne(**this->remoteCS), cMinusOne(**this->remoteCS);
            cOne = 1;
            cMinusOne = -1;
            gmp_randstate_t r_state;
            gmp_randinit_default (r_state);
            Plain r;
            
            
            this->computeT1s(snpIndices);

            for(int snpIdx: snpIndices) {
                SERVER_LOG("Comparing t1 values for SNP #" << snpIdx);

                t_1 = t_1n[snpIdx];

    #ifdef ENABLE_VALIDATION
                std::vector<RemoteCipher> innerProds;
    #endif
                
                std::vector<RemoteCipher> innerProdRs;
                std::vector<int> rBitLengths;
                RemoteCipher innerProd(**this->remoteCS);
                std::vector<mpz_class> allRs;
                
                int minBitLength = ceil(log2(totInst))+1; // Minimum bit length that needs comparing (extend if R[minBitLength] == 1)
                minBitLength++; // account for *2+1 (done to avoid 0)
                
                for(int j = 0; j < y_samples.size(); j++) {
                    auto y_sample = y_samples[j];
                    if((j+1)%100 == 0)
                        SERVER_LOG("Computing t1 " << (j+1) << "th sample");
                    innerProd = 0;
                    for(int i = 0; i < x_1n[snpIdx].size(); i++)
                        if(x_1n[snpIdx][i])
                            innerProd += y_sample[i];
                    
                    innerProd -= t_1;
                    
    #ifdef ENABLE_VALIDATION
                    //DEBUG: for verification purposes:
                    innerProds.push_back(innerProd);
    #endif
                    
                    // Remove zeros (so we can check for strict inequalities):
                    innerProd *= 2;
                    innerProd += cOne;
                    
                    //Todo: directly assign a random cipher
                    mpz_urandomb(r.raw().get_mpz_t(), r_state, (*this->remoteCS)->getBitLength()-1);
                    
                    innerProd += r;
                    innerProdRs.push_back(innerProd);
                    allRs.push_back(r.as_mpz_class());
                    
                    //DEBUG:
    //                DEBUG_cr = r;
    //                DEBUG_cRs.push_back(DEBUG_cr);
                    
    //                int rBitLength = 0;
    //                std::vector<LocalCipher> RBits;
    //                for(int b = 0; b < minBitLength || mpz_tstbit(r.raw().get_mpz_t(), b) == 1; b++) {
    //                    if(mpz_tstbit(r.raw().get_mpz_t(), b) == 1) {
    //                        RBits.push_back(one);
    //                        one.recrypt();
    //                    }
    //                    else {
    //                        RBits.push_back(zero);
    //                        zero.recrypt();
    //                    }
    //                    rBitLength++;
    //                }
    //                
    //                rBitLengths.push_back(rBitLength);
                }

                if(client == NULL)
                    client = new RPCClient<PaillierCryptosystemPrivateKey>(ALICE_SERVER_PORT);
                
                auto allProdBits = client->getRemoteCipherVectorVector(ProcBitDecomposition<PaillierCryptosystemPrivateKey>::procName(), innerProdRs, minBitLength);
                
                
                Plain p;
                RemoteCipher c(**this->remoteCS);
                std::vector<RemoteCipher> x1iResults;
                
                std::vector<Plain> x(minBitLength + 64);

                for (int i = 0; i < allRs.size(); i++) {
                    std::vector<RemoteCipher> y = allProdBits[i];

                    if((i+1) % 100 == 0)
                        SERVER_LOG("Running Comparison protocol for " << (i+1) << "th value");

                    //Encode from most significant to least:
                    for(int b = 0; b < y.size(); b++) {
                        p = mpz_tstbit(allRs[i].get_mpz_t(), y.size() - b - 1);
                        x[b] = p;
                    }
                    
                    std::vector<RemoteCipher> d(y.size(), **this->remoteCS);
                    std::vector<RemoteCipher> f(y.size(), **this->remoteCS);
                    std::vector<RemoteCipher> gamma(y.size(), **this->remoteCS);
                    RemoteCipher delta(**this->remoteCS);
                    RemoteCipher c(**this->remoteCS);
                    gamma[0] = 0;
                    
                    for(int j = 0; j < y.size(); j++) {
                        if(x[j] == 0) {
                            d[j] = y[j];
                            f[j] = y[j];
                        }
                        else {
                            d[j] = y[j] + cMinusOne;
                            f[j] = y[j] * (-1) + cOne;
                        }

                        if(j > 0) {
                            gamma[j] = gamma[j-1]*2;
                            gamma[j] += f[j];
                        }
                        //Todo: directly assign a random cipher
                        mpz_urandomb(r.raw().get_mpz_t(), r_state, (*this->remoteCS)->getBitLength());
                        
                        c = gamma[j];
                        c += cMinusOne;
                        
                        delta = c * r + d[j];
                        
                        x1iResults.push_back(delta);
                    }
                }
                
                allResults.push_back(x1iResults);
                
    #ifdef ENABLE_VALIDATION
                allResults.push_back(innerProds);
    #endif
            
            }
            ProcExactLogReg::packResult(req, allResults);
            SERVER_LOG("Remote procedure returning");

            // NAIVE COMPARISON PROTOCOL:
            /*
            std::vector<RemoteCipher> innerProds;
            RemoteCipher innerProd(**this->remoteCS), innerProdSign(**this->remoteCS);
            gmp_randstate_t r_state;
            gmp_randinit_default (r_state);
            Plain r;

            int j = 0;
            for(auto y_sample : y_samples) {
                innerProd = 0;
                for(int i = 0; i < x_1n[0].size(); i++)
                    if(x_1n[0][i])
                        innerProd += y_sample[i];
                
                innerProd -= t1;
                innerProds.push_back(innerProd);
                for(int i = 0; i <= totInst; i++) {
                    innerProdSign = innerProd - i;
                    
                    mpz_urandomb(r.raw().get_mpz_t(), r_state, (*this->remoteCS)->getBitLength());
                    innerProdSign *= r;
                    
                    innerProds.push_back(innerProdSign);
                }
                
                SERVER_LOG("Packed comparison data for sample #" << ++j);
            }

            ProcExactLogReg::packResult(req, innerProds);
            SERVER_LOG("Remote procedure returning");
            */
            
            
        };
    
    private:
        std::vector<std::vector<bool>>& x_1n;
        std::vector<RemoteCipher> y;
        std::vector<RemoteCipher> t_1n;
        
        boost::shared_mutex x_1Access;
        boost::mutex t_1Access;
        boost::shared_mutex yAccess;
        
        RPCClient<PaillierCryptosystemPrivateKey> *client = NULL;
    };
    
}

#endif /* defined(__PPFramework__ProcExactLogReg__) */

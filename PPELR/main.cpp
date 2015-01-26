//
//  main.cpp
//  PPFramework
//
//  Created by Dave on 11/28/14.
//  Copyright (c) 2014 Tsudalab. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>

#include "RPCServer.hpp"
#include "RPCClient.hpp"

#pragma GCC diagnostic push
#pragma clang diagnostic ignored "-Wshorten-64-to-32" // Stupid Boost lib is not yet completely 64bit-friendly
#include "boost/thread/thread.hpp"
#include "boost/function.hpp"
#include "boost/tokenizer.hpp"
#include "boost/format.hpp"
#include <boost/dynamic_bitset.hpp>
#pragma GCC diagnostic pop


//#define ENABLE_VALIDATION 1
//#define SANITY_CHECKS
//#define SHORTEN_VECTORS 1
#define LOG_DEBUG(_str) /* std::cout << "DEBUG: " << _str << std::endl */
#define LOG_MSG(_str) std::cout << _str << std::endl

#define ALICE_SERVER_PORT 1234
#define BOB_SERVER_PORT 1235

#define PAILLIER_ENCRYPTION_BITS 1024

#define TOT_Y_SAMPLES 1000
#define CHUNK_SIZE 50

#define NUM_X_1 100
#define PICK_RANDOM_X_1

#define TOT_REPEATS 1

#define TOT_THREADS 1

#define SIGNIFICANCE_LEVEL 0.01

#define EARLY_STOPPING 1

#include "PaillierCryptosystem.hpp"
#include "ProcExactLogReg.hpp"
#include "ProcBitDecomposition.hpp"

#include "utils.hpp"


using namespace ppf;

void testCryptosystem(void);
void runAlice(std::string, std::string, std::vector<double>&, std::vector<double>&, std::vector<std::vector<long>>&);
void runBob(char[]);
void runNonOblivious();
void readBoolVectorFromFile(const std::string, std::vector<bool>&, int pos = 0);
void ThrowError(const char*);
void ThrowError(const boost::format);
void TestEncryptSpeed();
void testPPF();

std::mutex bobMtx, elrWrite, threadMtx;
std::condition_variable cv, threadCv;

int activeThreads = 0;
bool bobReady = false;
RPCServer<PaillierCryptosystemPrivateKey>* bobServer = NULL;


int main(int argc, const char * argv[]) {

    boost::thread *bobThread, *aliceThread;

    std::vector<double> totTimesRet;
    std::vector<double> totShuffleTimesRet;
    std::vector<std::vector<long>> pValsRet;
#if TOT_REPEATS > 1
    for(int i = 0; i < TOT_REPEATS; i++) {
        std::cout << std::endl << "----------------------------------------" << std::endl << "BENCHMARK: repeating " << (i+1) << "/" << TOT_REPEATS << std::endl << std::endl;
#endif
        
        bobMtx.unlock();
        elrWrite.unlock();
        threadMtx.unlock();
        
        char x1FileName[] = "../../pplogreg/Data/Yamada-lab/HT/Stage1Data/X.%d.csv";
        bobThread = new boost::thread(runBob, x1FileName);
        
        std::string x2Filename = "../../pplogreg/Data/Yamada-lab/HT/Info/smoke.csv";
        std::string yFilename = "../../pplogreg/Data/Yamada-lab/HT/Info/smoke.csv";

        aliceThread = new boost::thread(runAlice, yFilename, x2Filename, boost::ref(totTimesRet), boost::ref(totShuffleTimesRet), boost::ref(pValsRet));
        aliceThread->join();
        
        if(bobServer) {
            delete bobServer;
            bobServer = NULL;
        }
        delete bobThread;
        delete aliceThread;
#if TOT_REPEATS > 1
    }
    double avgTotTime = std::accumulate(totTimesRet.begin(), totTimesRet.end(), 0.0) / totTimesRet.size();
    double sq_sum = std::inner_product(totTimesRet.begin(), totTimesRet.end(), totTimesRet.begin(), 0.0);
    double sdTotTime = sq_sum ? std::sqrt(sq_sum / totTimesRet.size() - avgTotTime * avgTotTime) : 0;

    double avgShuffleTime = std::accumulate(totShuffleTimesRet.begin(), totShuffleTimesRet.end(), 0.0) / totShuffleTimesRet.size();
    sq_sum = std::inner_product(totShuffleTimesRet.begin(), totShuffleTimesRet.end(), totShuffleTimesRet.begin(), 0.0);
    double sdShuffleTime = sq_sum ? std::sqrt(sq_sum / totShuffleTimesRet.size() - avgShuffleTime * avgShuffleTime) : 0;

    
    std::vector<std::vector<double>> pValsBySNPs(pValsRet[0].size());
    
    for (auto allSNPs : pValsRet)
        for (int i = 0; i < allSNPs.size(); i++)
            pValsBySNPs[i].push_back((double) allSNPs[i] / TOT_Y_SAMPLES);
    
    std::cout << std::endl << "----------------------------------------------" << std::endl;
    std::cout << "Ran " << totTimesRet.size() << " iterations:" << std::endl;
    std::cout << "Avg total time: " << avgTotTime / 1000 << " s. (" << sdTotTime/1000 << ")" << " | avg shuffle time: " << avgShuffleTime/1000 << " s. (" << sdShuffleTime/1000 << ")" << std::endl;
    
    std::cout << std::endl << "P-value averages:" << std::endl;
    for(int i = 0; i < pValsBySNPs.size(); i++) {
        std::cout << "SNP #" << (i+1) << ": ";
        
        double mean = std::accumulate(pValsBySNPs[i].begin(), pValsBySNPs[i].end(), 0.0) / pValsBySNPs[i].size();
        sq_sum = std::inner_product(pValsBySNPs[i].begin(), pValsBySNPs[i].end(), pValsBySNPs[i].begin(), 0.0);
        double sd = sq_sum ? std::sqrt(sq_sum / pValsBySNPs[i].size() - mean * mean) : 0;
        
        std::cout << mean << " (" << sd << ")" << std::endl;

    }
#endif
    
    return 0;
}




void runBob(char namePref[]) {
    try {
        
        bobServer = new RPCServer<PaillierCryptosystemPrivateKey>(NULL, BOB_SERVER_PORT);
        
        std::vector<std::vector<bool>> x_1n;

#ifdef PICK_RANDOM_X_1
        int tot = 690;
#else
        int tot = NUM_X_1;
#endif
        for (int i = 0; i < tot; i++) {
            std::vector<bool> v;
            auto name = boost::format(namePref) % i;
            readBoolVectorFromFile(name.str(), v, 4);
            
#ifdef SHORTEN_VECTORS  //DEBUG: (use shorter vector)
            v = std::vector<bool>(v.begin(), v.begin()+100);
#endif
            x_1n.push_back(v);
        }
        
        auto *elr = new ProcExactLogReg<PaillierCryptosystemPrivateKey>(x_1n);
        
        bobServer->registerProc(elr);
        
        bobReady = true;
        cv.notify_all();
        std::cout << "Starting Bob's server" << std::endl;
        bobServer->run();
        std::cout << "Closing Bob's server" << std::endl;

    }
    catch(const std::exception & ex) {
        std::cerr << "Error starting Server:" << std::endl << ex.what() << std::endl;
    }
}

void runAliceServer(RPCServer<PaillierCryptosystemPrivateKey>* server) {
    std::cout << "Starting Alice's server" << std::endl;
    server->run();
    std::cout << "Closig Alice's server" << std::endl;
}



void runAlice(std::string yFilename, std::string x2Filename, std::vector<double>& totTimesRet, std::vector<double>& totShuffleTimesRet, std::vector<std::vector<long>>& pValsRet) {
    try {
        std::cout << "Using Paillier encryption on " << PAILLIER_ENCRYPTION_BITS << " bits" << std::endl;
        typedef PaillierCryptosystemPrivateKey CS;
        CS localCS(PAILLIER_ENCRYPTION_BITS);
        
        RPCServer<CS> server(&localCS, ALICE_SERVER_PORT);
        server.registerProc(new ProcBitDecomposition<CS>);
        boost::thread serverThread(boost::bind(runAliceServer, &server));

        //*** Load y (outcome)
        std::vector<bool> y;
        readBoolVectorFromFile(yFilename, y);
        
#ifdef SHORTEN_VECTORS  //DEBUG: (use shorter y)
        y = std::vector<bool>(y.begin(), y.begin()+100);
#endif
        //*** Load x2 (covariates)
        std::vector<char> x2;
        std::vector<bool> x2_bool;
        
        readBoolVectorFromFile(x2Filename, x2_bool);
        
#ifdef SHORTEN_VECTORS  //DEBUG: (use shorter vector)
        x2_bool = std::vector<bool>(x2_bool.begin(), x2_bool.begin()+100);
#endif
        
        for (auto b : x2_bool)
            x2.push_back(b);
        
        if(y.size() != x2.size())
            ThrowError("categorical covariate x_2 does not have same length as outcome y");
        
        LOG_DEBUG("Read data: " << std::endl << "y   x2");
        for (int i = 0; i < y.size(); i++)
            LOG_DEBUG((int) y[i] << " | " << (int) x2[i]);
        std::cout << std::endl;
        
        // Reorder by values of categorical variable x_2:
        std::vector<int> p = sort_permutation(x2, [](char const& a, char const& b){ return(a < b); });
        
        x2 = apply_permutation(x2, p);
        
        // vector to reverse the permutation:
        std::vector<int> rev_p(p.size());
        for(int i = 0; i < p.size(); i++)
            rev_p[p[i]] = i;
        
        // Get counts of t_i, n_i for each value of x_2:
        std::vector<int> n_counts;
        int n = 0, prev_x2_val = x2[0];
        for (int i = 0; i < x2_bool.size(); i++) {
            LOG_DEBUG((int) x2[i] << " | " << (int) x2_bool[i]);
            if(x2[i] != prev_x2_val) {
                n_counts.push_back(n);
                n = 0;
                prev_x2_val = x2[i];
            }
            n++;
        }
        n_counts.push_back(n);

        boost::posix_time::ptime begin, end, beginShuffle, endShuffle;
        begin = boost::posix_time::microsec_clock::local_time();
        beginShuffle = boost::posix_time::microsec_clock::local_time();

        std::cout << "CLIENT: Encrypting real vector" << std::endl;

        std::vector<CS::Cipher> c_y;
        CS::Cipher c(localCS);
        for(auto b : y) {
            c = b;
            c_y.push_back(c);
        }
        
        std::cout << "CLIENT: Done encrypting real vector" << std::endl << std::endl;

        // Wait for Bob's server to be started:
        std::cout << "Waiting for Bob's server" << std::endl;
        std::unique_lock<std::mutex> lck(bobMtx);
        while (!bobReady) {  cv.wait(lck); }
        lck.unlock();
        std::cout << "Bob's server started: starting Alice's Client" << std::endl;
        RPCClient<CS> client(&localCS, BOB_SERVER_PORT);


        // Choosing x_1 indices to treat:
        std::vector<bool> useX_1Indices(690, false);
        int skip = 0;
#ifdef PICK_RANDOM_X_1
        //just skipping a random number of SNPs is enough for our needs
        skip = rand();
#endif
        for (int i = 0; i < NUM_X_1; i++)
            useX_1Indices[(i+skip)%useX_1Indices.size()] = true;

        std::vector<int> x_1Indices;

        // Result arrays:
        std::vector<long> pValTots(useX_1Indices.size());
        std::vector<long> pValTotsConfirm(useX_1Indices.size());

        // Reordered and re-encrypted version of y according to x_2 (used for samples):
        auto c_ord_y = apply_permutation(c_y, p);
        for(auto c : c_ord_y)
            c.recrypt();
        
        // Generate shuffled samples:
        auto engine = std::mt19937{};
        engine.seed((unsigned) std::chrono::system_clock::now().time_since_epoch().count());
        std::vector<std::vector<CS::Cipher>> c_y_samples;

        
        endShuffle = boost::posix_time::microsec_clock::local_time();
        auto totalShuffle = endShuffle - beginShuffle;

        for(int iter = 0; iter < (TOT_Y_SAMPLES/CHUNK_SIZE); iter ++) {
            x_1Indices.clear();
            for (int i = 0; i < useX_1Indices.size(); i++)
                if(useX_1Indices[i])
                    x_1Indices.push_back(i);

            std::cout << std::endl << "### Tot SNPs for iter: " << iter << " " << x_1Indices.size() << std::endl << std::endl;
            
            if(x_1Indices.size() == 0)
                break;
            
            std::cout << std::endl << "Iteration #" << (iter+1) << " (/" << (TOT_Y_SAMPLES/CHUNK_SIZE) << "): " << std::endl;
            std::cout << "Generating " << CHUNK_SIZE << " new samples" << std::endl;

            beginShuffle = boost::posix_time::microsec_clock::local_time();

            c_y_samples.clear();
            
            for (int i = 0; i < CHUNK_SIZE; i++) {
                int idx = 0;
                for (auto e : n_counts) {
                    std::shuffle(c_ord_y.begin() + idx, c_ord_y.begin() +idx+e, engine);
                    idx += e;
                }
                
                for (auto c_it : c_ord_y)
                    c_it.recrypt();
                
                c_y_samples.push_back(apply_permutation(c_ord_y, rev_p));
            }

            endShuffle = boost::posix_time::microsec_clock::local_time();
            totalShuffle += endShuffle - beginShuffle;

            if(iter > 0)
                c_y.clear(); // only need to send real Y the first time

            std::vector<std::vector<CS::Cipher>> resV;
            
            std::cout << "Calling ELR procedure" << std::endl;

            resV = client.getLocalCipherVectorVector(ProcExactLogReg<CS>::procName(), c_y, c_y_samples, x_1Indices);

            if(resV.size() > 0)
                std::cout << "Received " << resV.size() << " vectors back (x " << resV[0].size() << ")" << std::endl;
            
            for (int k = 0; k < x_1Indices.size(); k++) {
                auto testPVals = resV[k];
                
                int pValTot = 0;
                Plain t;
                
                for (int i = 0; i < testPVals.size(); i++) {
                    t = testPVals[i];
                    pValTot += (t == 1);
                }
                
                std::cout << "Partial result SNP #" << x_1Indices[k] << " (" << (k+1) << "/" << x_1Indices.size() << ") p-val: " << ((double) pValTot/c_y_samples.size()) << " (" << pValTot << ")" << std::endl;


                int pValIdx;
    #ifdef ENABLE_VALIDATION
                pValIdx = x_1Indices[k/2];
    #else
                pValIdx = x_1Indices[k];
    #endif
                pValTots[pValIdx] += pValTot;

                
    #ifdef EARLY_STOPPING
                if(pValTots[pValIdx] > (double) TOT_Y_SAMPLES*SIGNIFICANCE_LEVEL) {
                    std::cout << "Partial p-val total:" << pValTots[pValIdx] << "/" << TOT_Y_SAMPLES << " -> Non-significant: dropping from pool" << std::endl;
                    useX_1Indices[pValIdx] = false;
                    pValTots[pValIdx] = TOT_Y_SAMPLES;
                }
    #endif
                
    #ifdef ENABLE_VALIDATION
                k++;
                auto confirmPVals = resV[k];
                int pValTotConfirm = 0;
                for (int i = 0; i < confirmPVals.size(); i++) {
                    c = confirmPVals[i];
                    pValTotConfirm += (c >= 0);
                }
                
                if(pValTotConfirm != pValTot) {
                    if(abs(pValTotConfirm - pValTot) > 0.1 * c_y_samples.size())
                        std::cout << "### ERROR: Confirmation failed!" << std::endl;
                    
                    std::cout << "Confirmation pValTot: " << ((double) pValTot/c_y_samples.size()) << " (" << pValTot << ")" << std::endl;
                }
                
                pValTotsConfirm[x_1Indices[(k-1)/2]] += pValTotConfirm;
    #endif
            }
            
            if(x_1Indices.size() == 0)
                break;
        }
        
        
        std::cout << std::endl << "****************************" << std::endl << "FINAL SCORES:" << std::endl;
        for (int i = 0; i < pValTots.size(); i++) {
            std::cout << "SNP #" << i << " (" << (i+1) << "/" << pValTots.size() << ")" << " pVal >= " << ((double) pValTots[i] / TOT_Y_SAMPLES) << " (" << pValTots[i] << ") | validation: " << pValTotsConfirm[i] << std::endl;
        }
        
        end = boost::posix_time::microsec_clock::local_time();
        auto total = end - begin;
        std::cout << "Protocol took: " <<  total.total_milliseconds()/1000.0 << " s (Shuffle generation: " << totalShuffle.total_milliseconds()/1000.0 << " s)" << std::endl;

        pValsRet.push_back(pValTots);
        totShuffleTimesRet.push_back(totalShuffle.total_milliseconds());
        totTimesRet.push_back(total.total_milliseconds());
        
        serverThread.interrupt();
        boost::this_thread::sleep(boost::posix_time::seconds(1)); // let the server finish
        
    }
    catch(const std::exception & ex) {
        std::cerr << "Error running Client:" << std::endl << ex.what() << std::endl;
    }
}


void readBoolVectorFromFile(const std::string filename, std::vector<bool>& v, int pos) {
    std::ifstream input(filename);
    if(input.fail())
        ThrowError(boost::format("Can't open file: \"%s\": %s") % filename % strerror(errno));
        
    std::string line;
    std::string field;
    
    boost::char_separator<char> sep(",");
    
    while(std::getline(input, line)) {
        boost::tokenizer< boost::char_separator<char> > tok(line, sep);
        auto beg = tok.begin();
        for(int i = 0; i < pos && beg != tok.end(); ++beg, i++) {};
        v.push_back(!((*beg)[0] == '0'));
    }
}

void runNonOblivious(std::string yFilename, std::string x2Filename, char x1Filename[]) {
    // Load vector y
    
    std::vector<bool> y_bool;
    
    //Temp (should be reading from file):
    readBoolVectorFromFile(yFilename, y_bool);
    
    boost::dynamic_bitset<> y (y_bool.size());
    for (int i = 0; i < y_bool.size(); i++)
        y[i] = y_bool[i];
    
    // Load vector x2 (categorical covariate)
    std::vector<char> x2;
    std::vector<bool> x2_bool;

    readBoolVectorFromFile(x2Filename, x2_bool);
    for (auto b : x2_bool)
        x2.push_back(b);
    
    if(y_bool.size() != x2.size())
        ThrowError("categorical covariate x_2 does not have same length as outcome y");
    
    LOG_DEBUG("Read data: " << std::endl << "y   x2");
    for (int i = 0; i < y_bool.size(); i++)
        LOG_DEBUG((int) y_bool[i] << " | " << (int) x2[i]);
    std::cout << std::endl;
    
    // Reorder by values of categorical variable x_2:
    std::vector<int> p = sort_permutation(x2, [](char const& a, char const& b){ return(a < b); });
    
    x2 = apply_permutation(x2, p);
    auto ord_y_bool = apply_permutation(y_bool, p);
    
    // Order to reverse the permutation
    std::vector<int> rev_p(p.size());
    for(int i = 0; i < p.size(); i++)
        rev_p[p[i]] = i;
    
    
    // Get counts of t_i, n_i for each value of x_2:
    std::vector<int> n_counts;
    
    int n = 0, prev_x2_val = x2[0];
    for (int i = 0; i < ord_y_bool.size(); i++) {
        LOG_DEBUG((int) x2[i] << " | " << (int) ord_y_bool[i]);
        
        if(x2[i] != prev_x2_val) {
            n_counts.push_back(n);
            n = 0;
            prev_x2_val = x2[i];
        }
        
        n++;
    }
    n_counts.push_back(n);
    
    
//#ifdef SANITY_CHECKS
//    if(tot != y.size())
//        ThrowError("Sum of n_i must be equal to the length of y");
//#endif
    
    // Generate shuffled samples:
    std::vector<boost::dynamic_bitset<>> y_samples;
    
    auto engine = std::mt19937{};
    engine.seed((unsigned) std::chrono::system_clock::now().time_since_epoch().count());
    
    for (int i = 0; i < TOT_Y_SAMPLES; i++) {
        int idx = 0;
        for (auto e : n_counts) {
            std::shuffle(ord_y_bool.begin() + idx, ord_y_bool.begin() +idx+e, engine);
            idx += e;
        }
        
        boost::dynamic_bitset<> yb(ord_y_bool.size());
        auto rev_y_bool = apply_permutation(ord_y_bool, rev_p);
        for(int j = 0; j < ord_y_bool.size(); j++)
            yb[j] = rev_y_bool[j];
        
        y_samples.push_back(yb);
    }
    
    // Load vectors x1
    std::vector<boost::dynamic_bitset<>> x1n;
    int totSNPs = 689;
    
    for (int i = 0; i <= totSNPs; i++) {
        std::vector<bool> v;
        auto name = boost::format(x1Filename) % i;
        readBoolVectorFromFile(name.str(), v, 4);
        
        boost::dynamic_bitset<> vb(v.size());
        for (int i = 0; i < v.size(); i++)
            vb[i] = v[i];
        
        x1n.push_back(vb);
    }
    
    LOG_DEBUG(x1n.size() << " SNPs");

    boost::dynamic_bitset<> x1_y(x1n[0].size());
    std::vector<double> p_values(x1n.size());
    
    
    // run the test for each x1 (e.g. SNP)
    for (int i = 0; i < x1n.size(); i++) {
        auto x1 = x1n[i];
        
        // compute t1 (on the real y and x1)
        
        auto t1 = (x1 & y).count();
        
        LOG_DEBUG("#" << i << ": ");
        LOG_DEBUG("x1  y   x1_y");
        for(int i = 0; i < x1_y.size(); i++)
            LOG_DEBUG(x1[i] << " & " << y[i] << " = " << x1_y[i]);
        LOG_DEBUG("");


        // compare t1 to inner products of x1 and each sample:
        int p_value = 0;
        for (auto y_sample : y_samples) {
            auto tot = (x1 & y_sample).count();
            
            p_value += (tot >= t1);
        }

        if(p_value / (float) y_samples.size() < 0.05) {
            LOG_MSG("SNP #" << i << ": ");
            LOG_MSG("t1 == " << t1 << " (/" << x1_y.size() << ")");
            LOG_MSG("p-val: " << (p_value / (float) y_samples.size()) << " (" << p_value << " / " << y_samples.size() << ")");
            LOG_MSG("");
        }
        p_value /= (float) y_samples.size();
        p_values.push_back(p_value);
        
    }
}

void ThrowError(const boost::format fmt) {
    ThrowError(fmt.str().c_str());
}

void ThrowError(const char* msg) {
    std::cerr << "*** ERROR: " << msg << std::endl;
    exit(-1);
    //    throw msg;
}

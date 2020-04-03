#include "chronobench.h"
#include "codegen_FDFraser.h"
#include "codegen_FDFraserNAddr.h"
#include "cMMVII_Appli.h"
#include "MMVII_Derivatives.h"
#include "MMVII_FormalDerivatives.h"

#include <ceres/jet.h>

using ceres::Jet;
namespace  FD = NS_MMVII_FormalDerivative;

// ========== Define on Jets two optimization as we did on formal

template <typename T, int N>
inline Jet<T, N> square(const Jet<T, N>& f) {
  return Jet<T, N>(FD::square(f.a), 2.0*f.a * f.v);
}

template <typename T, int N>
inline Jet<T, N> cube(const Jet<T, N>& f) {
  return Jet<T, N>(FD::cube(f.a), 3.0*FD::square(f.a) * f.v);
}

using MMVII::cEpsNum;
#include "fraserformulatest.h"


typedef ChronoBench<3> Bench;

class cFraserCamColinearEq : cFraserCamColinearDef
{
    public :
       cFraserCamColinearEq(int aSzBuf, bool Show, int test, Bench &bench);

    private :
       typedef Jet<double,TheNbUk>  tJets;
       typedef cEpsNum<TheNbUk>     tEps;

       FD::cCoordinatorF<double>  mCFD;  /// Coordinator for formal derivative
};




cFraserCamColinearEq::cFraserCamColinearEq(int aSzBuf, bool Show, int test, Bench &bench) :
    // mCFD (aSzBuf,19,11) would have the same effect, but future generated code will be less readable
     mCFD  (aSzBuf,TheVNamesUnknowns,TheVNamesObs)
{
      static bool once=false;
   // In unknown, we set everything to zero exepct focal to 1
   mVUk[11] = 1.0;
   // In obs, we set the current matrix to Id
   mVObs[0] = mVObs[4] = mVObs[8] = 1;


   auto aVFormula = FraserCamColinearEq(mCFD.VUk(),mCFD.VObs());
   mCFD.SetCurFormulasWithDerivative(aVFormula);


   if (Show)
         mCFD.ShowStackFunc();

   int aNbTestTotal =  1e6; ///< Approximative number of Test
   int aNbTestWithBuf = aNbTestTotal/aSzBuf;  ///< Number of time we will fill the buffer
   aNbTestTotal = aNbTestWithBuf * aSzBuf; ///< Number of test with one equation

   // Here we initiate with "peferct" projection, to check something
   const std::vector<double> & aVUk  =  VUk(1.0,2.0,10.0);
   const std::vector<double> & aVObs =  VObs(0.101,0.2);

   // Make the computation with jets
   std::vector<tJets> aJetRes;
   if (test==0)  {
       bench.start();
       bench.start(1);
       std::vector<tJets>  aVJetUk;
       bench.stopStart(1,2);
       for (int aK=0 ; aK<TheNbUk ; aK++)
           aVJetUk.push_back(tJets(aVUk[aK],aK));
       bench.stopStart(2,3);
       for (int aK=0 ; aK<aNbTestTotal ; aK++)  {
           aJetRes = FraserCamColinearEq(aVJetUk,aVObs);
       }
       bench.stop(3);
       bench.stop();
   }

   std::vector<tEps> aEpsRes;
   if (test==1) {
       bench.start();
       bench.start(1);
        std::vector<tEps >  aVEpsUk;
        bench.stopStart(1,2);
        for (int aK=0 ; aK<TheNbUk ; aK++)
            aVEpsUk.push_back(tEps(aVUk[aK],aK));
        bench.stopStart(2,3);
        for (int aK=0 ; aK<aNbTestTotal ; aK++)
        {
            aEpsRes = FraserCamColinearEq(aVEpsUk,aVObs);
        }
        bench.stop(3);
        bench.stop();
   }

   // Make the computation with formal deriv buffered
   if (test == 2)
   {
       bench.start();
       bench.start(1);
         auto aVFormula = FraserCamColinearEq(mCFD.VUk(),mCFD.VObs());
         mCFD.SetCurFormulasWithDerivative(aVFormula);
         bench.stop(1);
       for (int aK=0 ; aK<aNbTestWithBuf ; aK++)
       {
           // Fill the buffers with data
           bench.start(2);
           for (int aKInBuf=0 ; aKInBuf<aSzBuf ; aKInBuf++)
               mCFD.PushNewEvals(aVUk,aVObs);
           bench.stopStart(2,3);
           // Evaluate the derivate once buffer is full
           mCFD.EvalAndClear();
           bench.stop(3);
       }
       bench.stop();
   }

   /////////// Gen Code
   typedef FDFraserNAddr<double> FDFraser;
   FDFraser fdFraser(0);
   if (test == 3) {
       bench.start();
       bench.start(1);
         fdFraser = FDFraser(aSzBuf);
         FDFraser::UkType aVVUk;
         FDFraser::ObsType aVVObs;
         for (size_t i=0; i<aVVUk.size(); i++)
             aVVUk[i] = aVUk[i];
         for (size_t i=0; i<aVVObs.size(); i++)
             aVVObs[i] = aVObs[i];
         bench.stop(1);
         for (int aK=0 ; aK<aNbTestWithBuf ; aK++)
         {
             // Fill the buffers with data
             bench.start(2);
               for (int aKi=0 ; aKi<aSzBuf ; aKi++)
                   fdFraser.pushNewEvals(aVVUk,aVVObs);
               bench.stopStart(2,3);
               fdFraser.evalAndClear();
               bench.stop(3);
         }
         bench.stop();
   }

// Verif result
   if (test==3 && !once) {
       std::cout << "TEST\n";
       once = true;
       std::vector<tJets>  aVJetUk;
       for (int aK=0 ; aK<TheNbUk ; aK++)
           aVJetUk.push_back(tJets(aVUk[aK],aK));
       aJetRes = FraserCamColinearEq(aVJetUk,aVObs);
       for (int aKVal=0 ; aKVal<int(aJetRes.size()) ; aKVal++)
       {
           if (fabs(aJetRes[aKVal].a != fdFraser.val(0,aKVal)) > fabs(aJetRes[aKVal].a) * 1e-15) {
               std::cerr << "CodeGen error: Jet[" << aKVal << "].a != gen.val(0," << aKVal<< ")";
               std::cerr << "  (" <<  aJetRes[aKVal].a - fdFraser.val(0,aKVal) << ")\n" ;

           }
           for (int aKVar=0;  aKVar< TheNbUk ; aKVar++)
           {
               if (fabs(aJetRes[aKVal].v[aKVar] - fdFraser.der(0,aKVal,aKVar)) > fabs(aJetRes[aKVal].v[aKVar]) * 1e-15) {
                   std::cerr << "CodeGen error: Jet[" << aKVal << "].v[" << aKVar << "] != gen.der(0," << aKVal<< "," << aKVar <<")";
                   std::cerr << "  (" << aJetRes[aKVal].v[aKVar] -  fdFraser.der(0,aKVal,aKVar) << ")\n" ;

               }
           }
       }

       auto aGenRes = fdFraser.result();
       for (size_t i=0; i< aGenRes.size(); i++) {
           for (size_t j=0; j<aGenRes[i].size(); j++) {
               if (aGenRes[i][j] != aGenRes[i][0]) {
                   std::cerr << "CodeGen error: aGenRes[" << i << "]["<< j << "] != aGenRes[" << i << "][0]";
                   std::cerr << "  (" << aGenRes[i][j] - aGenRes[i][0] << ")\n" ;
               }
           }
       }
   }
}


static std::ostream &operator<<(std::ostream &os, const Bench::Times &t)
{
    os << "[" ;
    for (size_t i=0; i<t.size(); i++)
        os << t[i] << (i<t.size()-1 ? "," : "");
    os << "]" ;
    return os;
}



static void TestFraserOneShot(size_t szBuf=100)
{
    Bench t;

    cFraserCamColinearEq (szBuf,true,3,t);
    std::cout << t.currents() << "\n";
}

static void TestFraserBenchMark(void)
{
    static const size_t NbTest=20;
    const std::string names[]= {"Jet","Eps","Buf","Gen"};

    for (int test=0; test <4; test++) {
        std::vector<size_t> nb_buffer;
        std::vector<size_t> nb_thread;
        if (test <2) {
            nb_buffer={1};
            nb_thread={1};
        } else {
            nb_buffer={1,10,20,50,100,200,500,700,1000,2000};
            nb_thread={1,2,3,4,5,6,7,8};
        }

        std::ofstream os(names[test] + ".txt");
        os << "0 ";
        for (auto &t : nb_thread) os << t << " ";
        for (auto &t : nb_thread) os << t << " ";
        for (auto &t : nb_thread) os << t << " ";
        os << "\n";
        for (auto &bufs : nb_buffer) {
            std::vector<Bench::Times> allTimes;
            for (auto &threads : nb_thread) {
                omp_set_num_threads(threads);
                Bench bench;
                std::string testName = names[test];
                if (test>1)
                    testName += " B:" + std::to_string(bufs) + " T:" + std::to_string(threads);
                testName += " ";

                while(1) {
                    for (size_t n=0; n<NbTest; n++)  {
                        cFraserCamColinearEq (bufs,false,test,bench);
                        std::cout << testName << bench.currents() << "\n";
                        bench.next();
                    }
                    auto filter=bench.filter();
                    std::cout << testName << "Main Mean:" << filter.mean << " sigma:" << filter.stddev << "\n";
                    if (filter.stddev > filter.mean / 100)
                        std::cout << testName << "Redoing this test ...\n";
                    break;
                }
                std::cout << testName << "Means: " << bench.fltMeans() << "\n\n";
                allTimes.push_back(bench.fltMeans());
            }
            for (auto &t : allTimes) os << t[0] << " ";
            for (auto &t : allTimes) os << t[2] << " ";
            for (auto &t : allTimes) os << t[1] << " ";
            os << "\n";
            os.flush();
        }
    }
}



/* ==================================================== */
/*                                                      */
/*          cAppli_TestCodeGen                          */
/*                                                      */
/* ==================================================== */
namespace MMVII
{


class cAppli_TestCodeGen : public cMMVII_Appli
{
     public :
        cAppli_TestCodeGen(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override {return anArgObl;}
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
     private :
        int mSizeBuf;
        int mThreads;
};


cCollecSpecArg2007 & cAppli_TestCodeGen::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt << AOpt2007(mSizeBuf,"b","Buffer Size",{})
                    << AOpt2007(mThreads,"t","Nb Thread Max",{});
}

cAppli_TestCodeGen::cAppli_TestCodeGen(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec) :
  cMMVII_Appli (aVArgs,aSpec),mSizeBuf(0),mThreads(0)
{
}


int cAppli_TestCodeGen::Exe()
{
    if (mSizeBuf == 0 || mThreads==0) {
        TestFraserBenchMark();
    } else {
        omp_set_num_threads(mThreads);
        TestFraserOneShot(mSizeBuf);
    }
   return EXIT_SUCCESS;
}


tMMVII_UnikPApli Alloc_TestCodeGen(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_TestCodeGen(aVArgs,aSpec));
}

cSpecMMVII_Appli TheSpecTestCodeGen
(
     "TestCodeGen",
      Alloc_TestCodeGen,
      "Code Generator test",
      {eApF::Perso},
      {eApDT::FileSys},
      {eApDT::Xml},
      __FILE__
);

}


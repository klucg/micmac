// File Automatically generated by eLiSe
#include "general/all.h"
#include "private/all.h"


class cEqHomogrSpaceInitX: public cElCompiledFonc
{
   public :

      cEqHomogrSpaceInitX();
      void ComputeVal();
      void ComputeValDeriv();
      void ComputeValDerivHessian();
      double * AdrVarLocFromString(const std::string &);
      void SetXL1(double);
      void SetXL2(double);
      void SetYL1(double);


      static cAutoAddEntry  mTheAuto;
      static cElCompiledFonc *  Alloc();
   private :

      double mLocXL1;
      double mLocXL2;
      double mLocYL1;
};

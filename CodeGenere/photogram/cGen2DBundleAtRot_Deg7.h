// File Automatically generated by eLiSe
#include "StdAfx.h"


class cGen2DBundleAtRot_Deg7: public cElCompiledFonc
{
   public :

      cGen2DBundleAtRot_Deg7();
      void ComputeVal();
      void ComputeValDeriv();
      void ComputeValDerivHessian();
      double * AdrVarLocFromString(const std::string &);
      void SetAmplAttR(double);
      void SetCentrAttR_x(double);
      void SetCentrAttR_y(double);
      void SetDepR1_x(double);
      void SetDepR1_y(double);
      void SetDepR2_x(double);
      void SetDepR2_y(double);
      void SetDepR3_x(double);
      void SetDepR3_y(double);
      void SetRotPt_x(double);
      void SetRotPt_y(double);


      static cAutoAddEntry  mTheAuto;
      static cElCompiledFonc *  Alloc();
   private :

      double mLocAmplAttR;
      double mLocCentrAttR_x;
      double mLocCentrAttR_y;
      double mLocDepR1_x;
      double mLocDepR1_y;
      double mLocDepR2_x;
      double mLocDepR2_y;
      double mLocDepR3_x;
      double mLocDepR3_y;
      double mLocRotPt_x;
      double mLocRotPt_y;
};

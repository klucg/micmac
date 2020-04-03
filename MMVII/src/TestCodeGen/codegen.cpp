#include "MMVII_formulagentpl.h"
#include "MMVII_FormalDerivatives.h"
#include "fraserformulatest.h"


int main(int , char **)
{
    NS_MMVII_FormalDerivative::cCoordinatorF<double>
            mCFD(0,cFraserCamColinearDef::TheVNamesUnknowns,cFraserCamColinearDef::TheVNamesObs);

    auto aVFormula = FraserCamColinearEq(mCFD.VUk(),mCFD.VObs());
    mCFD.SetCurFormulasWithDerivative(aVFormula);
    mCFD.genCodeNAddr("FDFraser");
    mCFD.genCode("FDFraser");
    return EXIT_SUCCESS;
}

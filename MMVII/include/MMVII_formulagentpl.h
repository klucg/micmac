#ifndef FORMULAGENTPL_H
#define FORMULAGENTPL_H

#include <vector>
#include <array>

// TODO: specialized exceptions

template<typename TypeElem,size_t NbUk, size_t NbObs, size_t NbRes, unsigned Interval>
class GenFuncTpl
{
public:
    typedef TypeElem Type;
    typedef std::vector<Type> VType;
    typedef std::vector<VType> VVType;
    typedef std::array<Type,NbUk> UkType;
    typedef std::array<Type,NbObs> ObsType;


    GenFuncTpl(size_t sizeBuf) :
        vvRes(NbRes,VType(sizeBuf)),
        vvUk(NbUk,VType(sizeBuf)),
        vvObs(NbObs,VType(sizeBuf)),
        mSizeBuf(sizeBuf),
        mInBuf(0)
    {}
    void pushNewEvals(const UkType &aVUk,const ObsType &aVObs)
    {
        if (mInBuf >= mSizeBuf)
            throw std::range_error("Codegen buffer overflow");
        for (size_t i=0; i<NbUk; i++)
            vvUk[i][mInBuf] = aVUk[i];
        for (size_t i=0; i<NbObs; i++)
            vvObs[i][mInBuf] = aVObs[i];
        mInBuf++;
    }
    void pushNewEvals(const VType &aVUk,const VType &aVObs)
    {
        if (mInBuf >= mSizeBuf)
            throw std::range_error("Codegen buffer overflow");
        if (aVUk.size() != vvUk.size())
            throw std::range_error("Codegen ukn size error");
        if (aVObs.size() != vvObs.size())
            throw std::range_error("Codegen obs size error");
        for (size_t i=0; i<NbUk; i++)
            vvUk[i][mInBuf] = aVUk[i];
        for (size_t i=0; i<NbObs; i++)
            vvObs[i][mInBuf] = aVObs[i];
        mInBuf++;
    }
    size_t bufferFree() const {return mSizeBuf - mInBuf;}
    size_t bufferSuze() const {return mSizeBuf;}

    void evalAndClear() { mInBuf=0; }
    Type val(size_t numBuf, size_t numVal) const
    {
        return vvRes[numVal * Interval][numBuf];
    }
    Type der(size_t numBuf, size_t numVal,size_t numDer) const
    {
        return vvRes[numVal * Interval + numDer + 1][numBuf];
    }
    const VVType& result() const
    {
        return vvRes;
    }
protected:
    VVType vvRes;
    VVType vvUk;
    VVType vvObs;
    size_t mSizeBuf;
    size_t mInBuf;
};


#endif // FORMULAGENTPL_H

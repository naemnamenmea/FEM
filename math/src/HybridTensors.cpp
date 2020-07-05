
#include "HybridTensors.hpp"

// tFEsSetOfTensor2& tFEsSetOfTensor2::Assign0()
//{
// for_each(pData.begin(),pData.end(),fDelete_object<Tensor2s>());
/* for (vector<Tensor2s*>::iterator ppcurrent=pData.begin(); ppcurrent<pData.end(); ++ppcurrent)
   {
	delete *ppcurrent;     *ppcurrent = nullptr;
   }*/
// return *this;
//}

/*tFEsSetOfTensor2& tFEsSetOfTensor2::operator+= (const tFEsSetOfTensor2& addend_)
{
#ifdef STRONGCHECK
 Assert(pData.size()==addend_.pData.size(), "different dimensions in tFEsSetOfTensor2::operator+=");
#endif//def STRONGCHECK
 vector<Tensor2s*>::iterator       ppcurrent =         pData.begin();
 vector<Tensor2s*>::const_iterator ppaddend  = addend_.pData.begin();
 for (; ppcurrent<pData.end(); ++ppcurrent,++ppaddend)
	if (*ppaddend != nullptr)
	  {
	   if (*ppcurrent == nullptr)
			  (*ppcurrent) = new Tensor2s(**ppaddend);
		else
			  (*ppcurrent)->operator+=(**ppaddend);
	  }
 return *this;
}

tFEsSetOfTensor2& tFEsSetOfTensor2::operator*= (real_t factor_)
{
 for (vector<Tensor2s*>::iterator ppcurrent=pData.begin(); ppcurrent<pData.end(); ++ppcurrent)
	if (*ppcurrent != nullptr)
	   (*ppcurrent)->operator*=(factor_);
 return *this;
}

tFEsSetOfTensor2& tFEsSetOfTensor2::operator/= (real_t divisor_)
{
 for (vector<Tensor2s*>::iterator ppcurrent=pData.begin(); ppcurrent<pData.end(); ++ppcurrent)
	if (*ppcurrent != nullptr)
	   (*ppcurrent)->operator/=(divisor_);
 return *this;
}*/

//==============================================================================
//!
//! \file CurlEquation.C
//!
//! \date Nov 24 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for the curl equations.
//!
//==============================================================================

#include "CurlEquation.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "TensorFunction.h"
#include "Vec3Oper.h"
#include "Utilities.h"


CurlEquation::CurlEquation(unsigned short int n, TimeIntegration::Method)
  : IntegrandBase(n)
{
  primsol.clear();
}


LocalIntegral* CurlEquation::getLocalIntegral (const std::vector<size_t>& nen,
                                               size_t, bool neumann) const
{
  return new MixedElmMats(nsd, nen);
}


bool CurlEquation::initElement (const std::vector<int>& MNPC,
                                const std::vector<size_t>& elem_sizes,
                                const std::vector<size_t>& basis_sizes,
                                LocalIntegral& elmInt)
{
  int ierr = 0;
  Vectors& vec = elmInt.vec;
  if (primsol.empty())
    return true;

  using IntVec = std::vector<int>;

  size_t ofs = 0;
  auto fstart = MNPC.begin();
  for (size_t i = 0; i < nsd; ++i) {
    ierr += utl::gather(IntVec(fstart,fstart+elem_sizes[i]),0,1,
                        primsol[0],vec[i],ofs,ofs);
    fstart += elem_sizes[i];
    ofs += basis_sizes[i];
  }

  if (ierr != 0)
    std::cerr << " *** CurlEquation::initElement: Detected " << ierr
              << " node numbers out of range." << std::endl;

  return ierr;
}


bool CurlEquation::evalIntMx (LocalIntegral& elmInt,
                              const MxFiniteElement& fe,
                              const TimeDomain& time,
                              const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  Vec3 f;
   if (source)
     f = (*source)(X);

  WeakOps::Curl(elMat.A, fe, 1.0);
  WeakOps::Mass(elMat.A, fe, 1.0);
  WeakOps::Source(elMat.b, fe, f);

  return true;
}


std::string CurlEquation::getField1Name (size_t i, const char* prefix) const
{
  const char* s[] = {"E_x", "E_y", "E_z"};

  if (i > 10)
    i -= 11;

  return prefix ? prefix + std::string(" ") + s[i] : s[i];
}


void CurlEquation::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;

  if (mode == SIM::RECOVERY)
    primsol.resize(1);
  else
    primsol.clear();
}


NormBase* CurlEquation::getNormIntegrand (AnaSol* asol) const
{
  return new CurlEquationNorm(*const_cast<CurlEquation*>(this), asol);
}


CurlEquationNorm::CurlEquationNorm (CurlEquation& p, AnaSol* a) :
  NormBase(p), anasol(a)
{
  nrcmp = p.getNoSpaceDim();
}


bool CurlEquationNorm::evalIntMx (LocalIntegral& elmInt,
                                  const MxFiniteElement& fe,
                                  const Vec3& X) const
{
  CurlEquation& problem = static_cast<CurlEquation&>(myProblem);
  ElmNorm& pNorm = static_cast<ElmNorm&>(elmInt);
  size_t nsd = problem.getNoSpaceDim();

  Vec3 Eh;
  Tensor grad(nsd);
  for (size_t i = 0; i < nsd; ++i) {
    Eh[i] = elmInt.vec[i].dot(fe.basis(1+i));
    for (size_t j = 1; j <= nsd; ++j)
      grad(i,j) = elmInt.vec[i].dot(fe.grad(i+1).getColumn(j));
  }

  auto&& curl = [nsd](const Tensor& grad)
                {
                  Vec3 c;
                  if (nsd == 3) {
                    c[0] = grad(3,2)-grad(2,3);
                    c[1] = grad(1,3)-grad(3,1);
                    c[2] = grad(2,1)-grad(1,2);
                  } else
                    c[2] = grad(2,1)-grad(1,2);

                  return c;
                };

  pNorm[0] += Eh*Eh*fe.detJxW;
  Vec3 cEh = curl(grad);
  pNorm[1] += cEh*cEh*fe.detJxW;

  if (anasol) {
    Vec3 E = (*anasol->getVectorSol())(X);
    pNorm[2] += E*E*fe.detJxW;
    E -= Eh;
    pNorm[3] += E*E*fe.detJxW;
    Vec3 cE = curl((*anasol->getVectorSecSol())(X));
    pNorm[4] += cE*cE*fe.detJxW;
    pNorm[5] += (cE-cEh)*(cE-cEh)*fe.detJxW;
  }

  return true;
}


size_t CurlEquationNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else if (group > 1)
    return 0;
  else
    return anasol ? 6 : 2;
}


std::string CurlEquationNorm::getName (size_t i, size_t j,
                                       const char* prefix) const
{
  static const char* s[6] = {
    "(E^h,E^h)^0.5",
    "(curl E^h,curl E^h)^0.5",
    "(E,E)^0.5",
    "(e^h, e^h)^0.5, e^h = E^h-E"
    "(curl E,curl E)^0.5",
    "(curl e^h, curl e^h)^0.5, e^h = E^h-E"
  };

  return prefix + std::string(" ") + s[j];
}

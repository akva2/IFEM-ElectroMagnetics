//==============================================================================
//!
//! \file Maxwell.C
//!
//! \date Sep 20 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Maxwells' equations.
//!
//==============================================================================

#include "Maxwell.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "Utilities.h"


Maxwell::Maxwell(unsigned short int n,
                TimeIntegration::Method method, bool m)
  : IntegrandBase(n), bdf(TimeIntegration::Order(method)), mixed(m)
{
  primsol.resize(1+bdf.getActualOrder());
}


LocalIntegral* Maxwell::getLocalIntegral(size_t nen, size_t, bool neumann) const
{
  return new StdElmMats(nsd, nen);
}


LocalIntegral* Maxwell::getLocalIntegral (const std::vector<size_t>& nen,
                                          size_t, bool neumann) const
{
  return new MixedElmMats(nsd, nen[0], nen[1]);
}


bool Maxwell::initElement (const std::vector<int>& MNPC,
                           LocalIntegral& elmInt)
{
  int ierr = 0;
  elmInt.vec.resize(2*primsol.size());
  for (size_t k = 0; k < primsol.size(); ++k) {
    Vector tmp(2*nsd*MNPC.size());
    ierr  += utl::gather(MNPC, 2*nsd, primsol[k], tmp);
    elmInt.vec[2*k].resize(nsd*MNPC.size());
    elmInt.vec[2*k+1].resize(nsd*MNPC.size());
    for (size_t i = 1; i <= MNPC.size(); ++i)
      for (size_t c = 1; c <= nsd; ++c) {
        elmInt.vec[2*k]((i-1)*nsd+c) = tmp((i-1)*2*nsd+c);
        elmInt.vec[2*k+1]((i-1)*nsd+c) = tmp((i-1)*2*nsd+nsd+c);
      }
  }

  if (ierr != 0)
    std::cerr << " *** Maxwell::initElement: Detected " << ierr/3
              << " node numbers out of range." << std::endl;

  return ierr == 0;
}


bool Maxwell::initElement (const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes,
                           LocalIntegral& elmInt)
{
  int ierr = 0;
  auto fstart = MNPC.begin() + elem_sizes[0];
  auto fend = fstart + elem_sizes[1];
  elmInt.vec.resize(2*primsol.size());
  for (size_t k = 0; k < primsol.size(); ++k) {
    ierr  += utl::gather(std::vector<int>(MNPC.begin(), fstart), nsd,
                         primsol[k], elmInt.vec[k*2]);
    ierr  += utl::gather(std::vector<int>(fstart, fend), nsd,
                         primsol[k], elmInt.vec[k*2+1]);
  }

  if (ierr != 0)
    std::cerr << " *** Maxwell::initElement: Detected " << ierr/3
              << " node numbers out of range." << std::endl;

  return ierr == 0;
}


bool Maxwell::evalInt (LocalIntegral& elmInt,
                      const FiniteElement& fe,
                      const TimeDomain& time,
                      const Vec3& X) const
{
  if (!props) {
    std::cerr << "No electro-magnetic properties specified." << std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  double eps = props->getDielectricConstant(X);
  double mu = props->getMagneticPermeability(X);
  double sigma = props->getConductivity(X);
  Vec3 Js;
   if (source)
     Js = (*source)(X);

  Vec3 E;
  for (int t = 1; t <= bdf.getOrder(); ++t)
    for (size_t k = 0; k < nsd; ++k)
      E[k] += eps*bdf[t]/time.dt * elmInt.vec[2*t].dot(fe.basis(1), k, nsd);

  Vec3 H;
  for (int t = 1; t <= bdf.getOrder(); ++t)
    for (size_t k = 0; k < nsd; ++k)
      H[k] += mu*bdf[t]/time.dt * elmInt.vec[2*t+1].dot(fe.basis(2), k, nsd);

//  WeakOps::Curl(elMat.A[3], fe, -1.0, 2, 1);
 // WeakOps::Curl(elMat.A[4], fe, 1.0, 1, 2);

  double fac = sigma;
  if (bdf.getOrder() > 0) {
    fac += eps*bdf[0]/time.dt;
    WeakOps::Mass(elMat.A[2], fe, mu*bdf[0]/time.dt, 2);
    WeakOps::Source(elMat.b[2], fe, H, 1.0, 2);
    E -= Js;
  }

  WeakOps::Mass(elMat.A[1], fe, fac, 1);
  WeakOps::Source(elMat.b[1], fe, E);

  return true;
}


std::string Maxwell::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11) {
    if (mixed)
      return nsd == 2 ? "E_x&&E_y" : "E_x&&E_y&&E_z";
    else
      return nsd == 2 ? "E_x&&E_y&&H_x&&H_y" : "E_x&&E_y&&E_z&&H_x&&H_y&&H_z";
  } else if (i == 12)
    return nsd == 2 ? "H_x&&H_y" : "H_x&&H_y&&H_z";

  const char* s[] = {"E_x", "E_y", "E_z", "H_x", "H_y", "H_z"};
  if (nsd == 2 && i >= 2)
    ++i;

  return prefix ? prefix + std::string(" ") + s[i] : s[i];
}


void Maxwell::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;

  if (mode == SIM::RECOVERY)
    primsol.resize(1);
  else if (mode == SIM::DYNAMIC)
    primsol.resize(1+bdf.getActualOrder());
  else
    primsol.clear();
}


NormBase* Maxwell::getNormIntegrand (AnaSol* asol) const
{
  return new MaxwellNorm(*const_cast<Maxwell*>(this),
                         asol ? asol->getVectorSol() : nullptr);
}


MaxwellNorm::MaxwellNorm (Maxwell& p, VecFunc* a) :
  NormBase(p), anasol(a)
{
  nrcmp = 6;
}


bool MaxwellNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                           const Vec3& X) const
{
  Maxwell& problem = static_cast<Maxwell&>(myProblem);
  ElmNorm& pNorm = static_cast<ElmNorm&>(elmInt);
  size_t nsd = problem.getNoSpaceDim();

  Vec3 Eh;
  for (size_t i = 0; i < nsd; ++i)
    Eh[i] = elmInt.vec[0].dot(fe.basis(1), i, nsd);

  Vec3 Hh;
  for (size_t i = 0; i < nsd; ++i)
    Hh[i] = elmInt.vec[1].dot(fe.basis(2), i, nsd);

  pNorm[0] += Eh*Eh*fe.detJxW;
  pNorm[1] += Hh*Hh*fe.detJxW;

  if (anasol) {
    std::vector<double> EH = anasol->getValue(X);
    Vec3 E(EH.data(),nsd);
    Vec3 H(EH.data()+nsd,nsd);
    pNorm[2] += E*E*fe.detJxW;
    E -= Eh;
    pNorm[3] += E*E*fe.detJxW;
    pNorm[4] += H*H*fe.detJxW;
    H -= Hh;
    pNorm[5] += H*H*fe.detJxW;
  }

  return true;
}


size_t MaxwellNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else if (group > 1)
    return 0;
  else
    return anasol ? 6 : 2;
}


std::string MaxwellNorm::getName (size_t i, size_t j, const char* prefix) const
{
  static const char* s[6] = {
    "(E^h,E^h)^0.5",
    "(H^h,H^H)^0.5",
    "(E,E)^0.5",
    "(e^h, e^h)^0.5, e^h = E^h-E"
    "(H,H)^0.5",
    "(e^h, e^h)^0.5, e^h = H^h-H"
  };

  return prefix + std::string(" ") + s[j];
}

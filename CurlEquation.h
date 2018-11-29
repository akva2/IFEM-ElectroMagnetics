// $Id$
//==============================================================================
//!
//! \file CurlEquation.h
//!
//! \date Sep 20 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for the curl equations.
//!
//==============================================================================

#ifndef _CURL_EQUATIONS_H
#define _CURL_EQUATIONSL_H

#include "BlockElmMats.h"
#include "ElectroMagneticProperties.h"
#include "IntegrandBase.h"
#include "CompatibleOperators.h"
#include "TimeIntUtils.h"
#include <array>
#include <cassert>
#include <memory>


/*!
  \brief Class representing the integrand of the curl equations.
*/

class CurlEquation : public IntegrandBase
{
public:
  using WeakOps = CompatibleOperators::Weak;     //!< Convenience renaming

  //! \brief Class representing the local block system of the curl equations.
  class MixedElmMats : public BlockElmMats {
  public:
    //! \brief Constructor.
    //! \param nsd Number of spatial dimensions
    //! \param ndof Number of basis functions
    MixedElmMats(size_t nsd,
                const std::vector<size_t>& ndof) : BlockElmMats(nsd, nsd)
    {
      this->resize(1+nsd*nsd, 1+nsd);
      for (size_t i = 1; i <= nsd; ++i)
        this->redim(i, ndof[i-1], 1, i);
      if (nsd == 3)
        assert(0);
      else {
        this->redimOffDiag(3,0);
        this->redimOffDiag(4,0);
      }
      this->redimNewtonMat();
    }

    //! \brief Empty destructor.
    virtual ~MixedElmMats() {}
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] method The time integration method to use
  //! \param[in] m \e true to use a mixed formulation
  CurlEquation(unsigned short int n, TimeIntegration::Method);

  //! \brief Empty destructor.
  virtual ~CurlEquation() {}

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                  size_t, bool neumann) const override;

  using IntegrandBase::initElement;
  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC Nodal point correspondance for the bases
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC,
                   const std::vector<size_t>& elem_sizes,
                   const std::vector<size_t>& basis_sizes,
                   LocalIntegral& elmInt) override;

  using IntegrandBase::evalIntMx;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const TimeDomain& time, const Vec3& X) const override;

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix = 0) const override;

  //! \brief Set current material.
  void setMaterial(ElectroMagneticProperties* mat) { props = mat; }

  //! \brief Set source function.
  void setSource(VecFunc* s) { source = s; }

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld = 2) const override { return fld > 1 ? 0 : nsd; }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol = nullptr) const override;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  void setMode(SIM::SolutionMode mode) override;

  void advanceStep() {}

protected:
  ElectroMagneticProperties* props = nullptr; //!< Electro-magnetic properties

  VecFunc* source = nullptr; //!< Source electric density
};


/*!
  \brief Class representing the integrand of curl equation norms.
*/

class CurlEquationNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] integrandType Integrand type flag
  //! \param[in] a The analytical solution (optional)
  CurlEquationNorm(CurlEquation& p, AnaSol* a = nullptr);
  //! \brief Empty destructor.
  virtual ~CurlEquationNorm() {}

  using NormBase::evalIntMx;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const Vec3& X) const override;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group = 0) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t i, size_t j, const char* prefix) const override;

private:
  AnaSol* anasol; //!< Analytical fields
};

#endif

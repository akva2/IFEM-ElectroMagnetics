// $Id$
//==============================================================================
//!
//! \file Maxwell.h
//!
//! \date Sep 20 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Maxwells' equations.
//!
//==============================================================================

#ifndef _MAXWELL_H
#define _MAXWELL_H

#include "BDF.h"
#include "BlockElmMats.h"
#include "ElectroMagneticProperties.h"
#include "IntegrandBase.h"
#include "TimeIntUtils.h"
#include "EqualOrderOperators.h"
#include <array>
#include <memory>


/*!
  \brief Class representing the integrand of Maxwells' equations.

  \details Time stepping is done using BDF1/BDF2
*/

class Maxwell : public IntegrandBase
{
public:
  using WeakOps = EqualOrderOperators::Weak;     //!< Convenience renaming

  //! \brief Class representing the local block system of Maxwells equations.
  class StdElmMats : public BlockElmMats {
  public:
    //! \brief Constructor.
    //! \param nsd Number of spatial dimensions
    //! \param ndof Number of basis functions
    StdElmMats(size_t nsd, size_t ndof) : BlockElmMats(2, 1)
    {
      this->resize(5, 3);
      this->redim(1, ndof, nsd, 1);
      this->redim(2, ndof, nsd, 1);
      this->redimOffDiag(3,0);
      this->redimOffDiag(4,0);
      this->redimNewtonMat();
    }

    //! \brief Empty destructor.
    virtual ~StdElmMats() = default;
  };

  //! \brief Class representing the local block system of Maxwells equations.
  class MixedElmMats : public BlockElmMats {
  public:
    //! \brief Constructor.
    //! \param nsd Number of spatial dimensions
    //! \param ndof Number of basis functions
    MixedElmMats(size_t nsd, size_t ndof1, size_t ndof2) : BlockElmMats(2, 2)
    {
      this->resize(5, 3);
      this->redim(1, ndof1, nsd, 1);
      this->redim(2, ndof2, nsd, 2);
      this->redimOffDiag(3,0);
      this->redimOffDiag(4,0);
      this->redimNewtonMat();
    }

    //! \brief Empty destructor.
    virtual ~MixedElmMats() = default;
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] method The time integration method to use
  //! \param[in] m \e true to use a mixed formulation
  Maxwell(unsigned short int n = 3,
          TimeIntegration::Method method = TimeIntegration::BE, bool m = false);

  //! \brief Empty destructor.
  virtual ~Maxwell() = default;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen,
                                  size_t, bool neumann) const override;

  //! \brief Returns a local integral container for the given element
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                  size_t, bool neumann) const override;

  using IntegrandBase::initElement;
    //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt) override;

  //! \brief Initializes current element for numerical integration (mixed).
  //! \param[in] MNPC Nodal point correspondance for the bases
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC,
                   const std::vector<size_t>& elem_sizes,
                   const std::vector<size_t>& basis_sizes,
                   LocalIntegral& elmInt) override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

  using IntegrandBase::evalIntMx;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const TimeDomain& time, const Vec3& X) const override
  {
    return this->evalInt(elmInt,fe,time,X);
  }

  //! \brief Advances the time stepping scheme.
  void advanceStep() { bdf.advanceStep(); }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix = 0) const override;

  //! \brief Set current material.
  void setMaterial(ElectroMagneticProperties* mat) { props = mat; }

  //! \brief Set source function.
  void setSource(VecFunc* s) { source = s; }

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld = 2) const override { return fld > 1 ? 0 : 2*nsd; }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol = nullptr) const override;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  void setMode(SIM::SolutionMode mode) override;

protected:
  TimeIntegration::BDF bdf; //!< BDF helper class
  ElectroMagneticProperties* props = nullptr; //!< Electro-magnetic properties

  VecFunc* source = nullptr; //!< Source electric density

  double eps = 1.0; //!< Dielectric constant
  double sigma = 1.0; //!< Conductivity of medium
  bool mixed = false; //!< \e true for mixed formulation
};


/*!
  \brief Class representing the integrand of Maxwells' norms.
*/

class MaxwellNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] integrandType Integrand type flag
  //! \param[in] a The analytical solution (optional)
  MaxwellNorm(Maxwell& p, VecFunc* a = nullptr);
  //! \brief Empty destructor.
  virtual ~MaxwellNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  using NormBase::evalIntMx;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const TimeDomain& time, const Vec3& X) const override
  { return this->evalInt(elmInt,fe,time,X); }

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
  VecFunc* anasol; //!< Analytical fields
};

#endif

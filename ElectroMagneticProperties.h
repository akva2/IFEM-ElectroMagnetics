// $Id$
//==============================================================================
//!
//! \file ElectroMagneticProperties.h
//!
//! \date Sep 20 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for electro-magnetic properties.
//!
//==============================================================================

#ifndef _ELECTRO_MAGNETIC_PROPERTIES_H
#define _ELECTRO_MAGNETIC_PROPERTIES_H

class TiXmlElement;

#include "Function.h"


/*!
  \brief Class representing electro-magnetic properties for Maxwells problems.
*/

class ElectroMagneticProperties
{
public:
  //! \brief Parses material parementers from an XML element.
  void parse(const TiXmlElement*);

  //! \brief Prints out fluid properties to the log stream.
  void printLog() const;

  //! \brief Returns the mass density.
  double getConductivity(const Vec3& X) const { return sigma.evaluate(X); }

  //! \brief Returns thermal diffusivity.
  double getDielectricConstant(const Vec3& X) const { return eps.evaluate(X); }

  double getMagneticPermeability(const Vec3& X) const { return mu.evaluate(X); }

protected:
  //! \brief Helper template for wrapping a constant/function pair.
  template<class Function> struct FuncConstPair
  {
    Function* function;                 //!< Function definition
    typename Function::Output constant; //!< Constant

    //! \brief Default constructor.
    FuncConstPair() { function = nullptr; constant = 0.0; }

    //! \brief Parses an XML element. Specialized per function type.
    Function* parse(const char*, const std::string&) { return nullptr; }

    //! \brief Evaluates the function at the given point \b X.
    typename Function::Output evaluate(const typename Function::Input& X) const
    {
      return function ? (*function)(X) : constant;
    }
  };

  FuncConstPair<RealFunc> eps; //!< Dielectric constant
  FuncConstPair<RealFunc> mu; //!< Magnetic permeability of medium
  FuncConstPair<RealFunc> sigma; //!< Conductivity of medium
};

#endif

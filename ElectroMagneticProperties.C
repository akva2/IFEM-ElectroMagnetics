// $Id$
//==============================================================================
//!
//! \file ElectroMagneticProperties.C
//!
//! \date Sep 20 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for electro-magnetic properties.
//!
//==============================================================================


#include "ElectroMagneticProperties.h"

#include "IFEM.h"
#include "Functions.h"
#include "tinyxml.h"
#include "Utilities.h"
#include "Vec3.h"

/*!
  \brief Parses an XML element. Specialization for RealFunc.
*/

template<>
RealFunc* ElectroMagneticProperties::FuncConstPair<RealFunc>::parse (const char* val,
                                                                     const std::string& type)
{
  return utl::parseRealFunc(val, type);
}



/*!
  \brief Template function to parse a property value from an XML-element.
*/

template<class T>
static bool propertyParse (ElectroMagneticProperties::FuncConstPair<T>& data,
                           const TiXmlElement* elem,
                           const char* attr, const char* tag)
{
  if (utl::getAttribute(elem,attr,data.constant))
    return true;

  const TiXmlElement* child = elem->FirstChildElement(tag);
  const TiXmlNode* aval = child ? child->FirstChild() : nullptr;
  if (!aval) return false;

  IFEM::cout <<" ";
  std::string type;
  utl::getAttribute(child,"type",type,true);
  data.function = data.parse(aval->Value(),type);

  return data.function != nullptr;
}


void ElectroMagneticProperties::parse(const TiXmlElement* elem)
{
  propertyParse(eps, elem, "eps", "dielectric");
  propertyParse(mu, elem, "mu", "magneticpermeability");
  propertyParse(sigma, elem, "sigma", "conductivity");
}


void ElectroMagneticProperties::printLog() const
{
  IFEM::cout << "\tElectro-magnetic properties:";
  IFEM::cout << "\n\t\t\t Dielectric constant, eps = " << eps.constant;
  IFEM::cout << "\n\t\t\tMagnetic permeability, mu = " << mu.constant;
  IFEM::cout << "\n\t\t\t      Conductivity, sigma = " << sigma.constant;
}

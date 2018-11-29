// $Id$
//==============================================================================
//!
//! \file SIMElectroMagnetics.h
//!
//! \date Sep 20 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for Isogeometric FE analysis of electromagnetics.
//!
//==============================================================================

#ifndef _SIM_ELECTRO_MAGNETICS_H_
#define _SIM_ELECTRO_MAGNETICS_H_

#include "IFEM.h"
#include "AnaSol.h"
#include "DataExporter.h"
#include "ElectroMagneticProperties.h"
#include "Functions.h"
#include "SIMenums.h"
#include "SIMsolution.h"
#include "TimeStep.h"
#include "tinyxml.h"
#include "TimeIntUtils.h"
#include "Utilities.h"


/*!
  \brief Driver class for isogeometric FE analysis of electromagnetics problems.
*/

template<class Dim, class Integrand>
class SIMElectroMagnetics : public Dim, public SIMsolution
{
public:
  //! \brief Mixed constructor.
  SIMElectroMagnetics(TimeIntegration::Method method) :
    Dim(std::vector<unsigned char>(Dim::dimension, 1)), em(Dim::dimension, method)
  {
    Dim::myProblem = &em;
  }

  //! \brief Destructor.
  virtual ~SIMElectroMagnetics()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "ElectroMagnetics"; }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    // Write solution fields
    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (!this->writeGlvS(this->solution.front(),iDump,nBlock,tp.time.t))
      return false;

    return this->writeGlvStep(iDump,tp.time.t);
  }

  //! \brief Initializes for time-dependent simulation.
  bool init()
  {
    this->initSolution(this->getNoDOFs(), em.getNoSolutions());
    this->registerField("solution1",this->SIMsolution::getSolution(0));
    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    this->setMode(tp.multiSteps() ? SIM::DYNAMIC : SIM::STATIC);
    if (tp.multiSteps() && Dim::msgLevel >= 0)
      IFEM::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    Vector dummy;
    if (!this->updateDirichlet(tp.time.t, &dummy))
      return false;

    if (!this->assembleSystem(tp.time, this->solution))
      return false;

    if (!this->solveSystem(this->solution.front(),1))
      return false;

    return true;
  }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep& tp)
  {
    em.advanceStep();
    this->pushSolution();
    return true;
  }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override
  {
    if (strcasecmp(elem->Value(),"electromagnetics"))
      return this->Dim::parse(elem);

    bool result = true;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"emproperties")) {
        IFEM::cout <<"  Parsing <"<< child->Value() <<">"<< std::endl;
        int code = this->parseMaterialSet(child,mVec.size());
        IFEM::cout <<"\tMaterial code "<< code <<":";
        mVec.push_back(std::make_unique<ElectroMagneticProperties>());
        mVec.back()->parse(child);
        if (code == 0)
          em.setMaterial(mVec.back().get());
      } else if (!strcasecmp(child->Value(),"source")) {
        std::string type;
        utl::getAttribute(child, "type", type, true);

        if (type == "expression" && child->FirstChild()) {
          std::string variables;
          utl::getAttribute(child,"variables",variables);
          IFEM::cout << "\n\tSource function:";
          VecFunc *func = utl::parseVecFunc(child->FirstChild()->Value(), type, variables);
          if (!variables.empty())
            IFEM::cout << "\n\t\tVariables: "
                       << variables;
          IFEM::cout << std::endl;
          em.setSource(func);
        }
      } else if (!strcasecmp(child->Value(),"anasol")) {
        std::string type;
        utl::getAttribute(child,"type",type,true);
        if (type == "expression") {
          IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
          this->mySol = new AnaSol(child,false);
        }
        else
          std::cerr <<"  ** SIMElectroMagnetics::parse: Invalid analytical solution "
            << type <<" (ignored)"<< std::endl;
      } else
        result &= this->Dim::parse(child);

    return true;
  }

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(SIMsolution::SerializeMap& data) const override
  {
    return this->saveSolution(data,this->getName());
  }

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const SIMsolution::SerializeMap& data) override
  {
    if (!this->restoreSolution(data,this->getName()))
      return false;

    em.advanceStep();
    return true;
  }

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  bool initMaterial(size_t propInd) override
  {
    if (propInd >= mVec.size())
      propInd = mVec.size()-1;

    em.setMaterial(mVec[propInd].get());
    return true;
  }

  //! \brief Prints a summary of the calculated solution to std::cout.
  //! \param[in] solvec The solution vector
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  void printSolutionSummary(const Vector& solvec, int, const char*,
                            std::streamsize outPrec) override
  {
    const size_t nsd = this->getNoSpaceDim();
    size_t iMax[2*nsd];
    double dMax[2*nsd];
    double dNorm;
    if (this->getNoBasis() > 1) {
      double dNormE = this->solutionNorms(solvec, dMax, iMax, nsd, 'D');
      double dNormH = this->solutionNorms(solvec, dMax+nsd, iMax+nsd, nsd, 'P');
      dNorm = sqrt(dNormE*dNormE+dNormH*dNormH);
    } else
      dNorm = this->solutionNorms(solvec, dMax, iMax, 2*nsd, 'D');

    std::stringstream str;
    if (Dim::adm.getProcId() == 0)
    {
      if (outPrec > 0) str.precision(outPrec);

      str <<"  Primary solution summary: L2-norm              : "<< utl::trunc(dNorm);

      char D = 'X';
      for (size_t d = 0; d < nsd; d++, D++)
        if (utl::trunc(dMax[d]) != 0.0)
          str <<"\n                            Max "<< char('X'+d)
              <<" electric field : "<< dMax[d] <<" node "<< iMax[d];
      for (size_t d = nsd; d < 2*nsd; d++, D++)
        if (utl::trunc(dMax[d]) != 0.0)
          str <<"\n                            Max "<< char('X'+(d-nsd))
              <<" magnetic field : "<< dMax[d] <<" node "<< iMax[d];
    }

    IFEM::cout << str.str() << std::endl;
  }

  //! \brief Evaluates and prints out solution norms.
  void printFinalNorms(const TimeStep& tp)
  {
    Vectors gNorm;
    this->setMode(SIM::RECOVERY);
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,this->solution,gNorm))
      return;
    else if (gNorm.empty())
      return;

    IFEM::cout << ">>> Norm summary for Maxwells <<<" << std::endl;
    IFEM::cout <<"  L2 norm |E^h| = (E^h,E^h)^0.5       : "<< gNorm[0](1);
    IFEM::cout <<"\n  L2 norm |H^h| = (H^h,H^h)^0.5       : "<< gNorm[0](2);
    if (this->haveAnaSol() && gNorm[0].size() >= 6)
      IFEM::cout <<"\n  L2 norm |E|   = (E,E)^0.5           : "<< gNorm[0](3)
                 <<"\n  L2 norm |e|   = (e,e)^0,5, e=E-E^h  : "<< gNorm[0](4)
                 <<"\n  L2 norm |H|   = (H,H)^0.5           : "<< gNorm[0](5)
                 <<"\n  L2 norm |e|   = (e,e)^0,5, e=H-H^h  : "<< gNorm[0](6);
    IFEM::cout << std::endl;
  }

  //! \brief Register fields for simulation result export.
  void registerFields(DataExporter& exporter)
  {
    int results = DataExporter::PRIMARY;

    if (!Dim::opt.pSolOnly && em.getNoFields(2) > 0)
      results |= DataExporter::SECONDARY;

    if (Dim::opt.saveNorms)
      results |= DataExporter::NORMS;

    exporter.registerField("u","solution",DataExporter::SIM,results);
    exporter.setFieldValue("u", this, &this->solution.front());
  }

private:
  Integrand em; //!< Problem integrand
  std::vector<std::unique_ptr<ElectroMagneticProperties>> mVec; //!< Material properties
};

#endif

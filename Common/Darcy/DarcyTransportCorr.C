// $Id$
//==============================================================================
//!
//! \file DarcyTransportCorr.C
//!
//! \date Jan 14 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Darcy transport correction problem.
//!
//==============================================================================

#include "DarcyTransportCorr.h"

#include "BlockElmMats.h"
#include "EqualOrderOperators.h"
#include "FiniteElement.h"
#include "Function.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "SIMbase.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <array>
#include <cmath>
#include <vector>


DarcyTransportCorr::DarcyTransportCorr (unsigned short int n, int torder) :
  IntegrandBase(n)
{
  npv = n;
  primsol.resize(1);
}


DarcyTransportCorr::~DarcyTransportCorr() = default;


LocalIntegral*
DarcyTransportCorr::getLocalIntegral (const std::vector<size_t>& nen,
                                      size_t, bool neumann) const
{
  BlockElmMats* result = new BlockElmMats(3, true ? 3 : 2);

  result->resize(NMAT, NVEC);
  result->redim(1, nen[0], nsd);
  result->redim(2, nen[1], 1, -2);
  result->redim(3, 1, true ? 1 : 0, -3);
  result->redimOffDiag(ql, 1);
  result->redimOffDiag(lg, 1);
  result->finalize();

  return result;
}

LocalIntegral* DarcyTransportCorr::getLocalIntegral (size_t nen, size_t,
                                                     bool neumann) const
{
  ElmMats* result = new ElmMats();
  result->resize(3,3);
  result->redim(nen*nsd);

  return result;
}

bool DarcyTransportCorr::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                  const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  const double C = (*observed_C)(X);
  const Vec3 dC = observed_C->gradient(X);
  const double dCdt = observed_C->timeDerivative(X);
  const Vec3 q = (*input_q)(X);
  const double f = (*this->input_source)(X);

  EqualOrderOperators::Weak::Mass(elMat.A[0], fe);
  for (size_t i = 1; i <= fe.N.size(); ++i) {
    for (size_t j = 1; j <= fe.N.size(); ++j)
      for (int d = 1; d <= nsd; ++d) {
          for (int d2 = 1; d2 <= nsd; ++d2) {
            const double mass_term = fe.dNdX(i,d) * fe.dNdX(j,d2);
            const double transport_term = (dC[d2-1]*fe.N(j) + C * fe.dNdX(j,d2)) *
                                          (dC[d-1]*fe.N(i) + C * fe.dNdX(i,d));
            elMat.A[1]((i-1)*nsd + d, (j-1)*nsd + d2) += mass_term * fe.detJxW;
            elMat.A[2]((i-1)*nsd + d, (j-1)*nsd + d2) += transport_term * fe.detJxW;
          }
      }
  }

  for (size_t i = 1; i <= fe.N.size(); ++i) {
    for (int d = 1; d <= nsd; ++d) {
      double transport_term = 0.0;
      for (int d2 = 1; d2 <= nsd; ++d2)
        transport_term += (f - dCdt) * (dC[d2-1]*fe.N(i) + C * fe.dNdX(i,d2));
      elMat.b[2]((i-1)*nsd + d) += transport_term * fe.detJxW;
    }
  }

  EqualOrderOperators::Weak::Source(elMat.b[0], fe, q);

  return true;
}


bool DarcyTransportCorr::evalIntMx (LocalIntegral& elmInt, const MxFiniteElement& fe,
                                  const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  const double C = (*observed_C)(X);
  const Vec3 dC = observed_C->gradient(X);
  const double dCdt = observed_C->timeDerivative(X);
  const Vec3 q = (*input_q)(X);
  const double R = (*this->input_source)(X);//this->residual(X);

  EqualOrderOperators::Weak::Mass(elMat.A[qq], fe);
  for (size_t i = 1; i <= fe.basis(2).size(); ++i)
    for (size_t j = 1; j <= fe.N.size(); ++j)
      for (int d = 1; d <= nsd; ++d) {
        elMat.A[ql]((j-1)*nsd+d, i) += (dC[d-1]*fe.N(j) + C * fe.dNdX(j,d)) * fe.basis(2)(i) * fe.basis(1)(j) * fe.detJxW;
      }

  EqualOrderOperators::Weak::Source(elMat.b[Fq], fe, q);
  EqualOrderOperators::Weak::Source(elMat.b[Fl], fe, R - dCdt);

  EqualOrderOperators::Weak::ItgConstraint(elMat.A[lg], fe, 1.0, 2);

  return true;
}


bool DarcyTransportCorr::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                              const Vec3& X, const Vec3& normal) const
{
  return true;
}


bool DarcyTransportCorr::finalizeElement (LocalIntegral& A)
{
  ElmMats& E = static_cast<ElmMats&>(A);

  E.A[0].add(E.A[1], alpha);
  E.A[0].add(E.A[2], beta);

  E.b[0].add(E.b[1], alpha);
  E.b[0].add(E.b[2], beta);

  return true;
}


bool DarcyTransportCorr::evalSol (Vector& s, const FiniteElement& fe,
                                  const Vec3& X,
                                  const std::vector<int>& MNPC) const
{
  s.resize(this->getNoFields(2));
  s[0] = (*observed_C)(X);
  const auto q = (*input_q)(X);
  const double f = (*input_source)(X);
  s[1] = this->residual(X);
  s[2] = f;
  for (size_t i = 0; i < nsd; ++i)
    s[3+i] = q[i];
  return true;
}


bool DarcyTransportCorr::suppressOutput (size_t i, ASM::ResultClass type) const
{
  return (i == 2 && type == ASM::PRIMARY);
}


std::string DarcyTransportCorr::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11 && nsd == 2)
    return "q_x&&q_y";
  if (i == 11 && nsd == 3)
    return "q_x&&q_y&&q_z";

  const char postfix = 'x' + i;
  std::string name = "q_";
  name += postfix;

  return prefix ? prefix + std::string(" ") + name : name;
}

std::string DarcyTransportCorr::getField2Name (size_t i, const char* prefix) const
{
    const auto names2 = std::array {
        "observed_C",
        "residual",
        "source",
        "input q_x",
        "input q_y",
    };
    const auto names3 = std::array {
        "observed_C",
        "residual",
        "source",
        "input q_x",
        "input q_y",
        "input q_z",
    };

    const std::string res = nsd == 2 ? names2[i] : names3[i];
    return prefix ? prefix + std::string(" ") + res : res;
}


void DarcyTransportCorr::setObservedConcentration (std::unique_ptr<RealFunc> obs_c)
{
  observed_C = std::move(obs_c);
}


void DarcyTransportCorr::setInputSource (std::unique_ptr<RealFunc> f)
{
  input_source = std::move(f);
}


void DarcyTransportCorr::setInputVelocity (std::unique_ptr<VecFunc> inp_q)
{
  input_q = std::move(inp_q);
}


double DarcyTransportCorr::residual (const Vec3& X) const
{
  const double C = (*observed_C)(X);
  const Vec3 q = (*input_q)(X);
  const Vec3 grad_C = observed_C->gradient(X);
  const auto grad_q = input_q->gradient(X);
  const double f = (*input_source)(X);
  const double eq = -(grad_C[0]*q[0] + C*grad_q(1,1) +
                      grad_C[1]*q[1] + C*grad_q(2,2));
  return f + eq;
}

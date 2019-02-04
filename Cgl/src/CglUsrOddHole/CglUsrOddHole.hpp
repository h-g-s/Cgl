#ifndef _CglUsrOddHole_h_
#define _CglUsrOddHole_h_

#include "CglCutGenerator.hpp"

class CglUsrOddHole : public CglCutGenerator {
public:
  static int sepCuts;
  static double sepTime;

  CglUsrOddHole();
  CglUsrOddHole(const CglUsrOddHole &rhs);
  virtual CglCutGenerator *clone() const;
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info = CglTreeInfo());
  virtual ~CglUsrOddHole();

private:
  void fillColSolution(const OsiSolverInterface &si, double colSol[]);
  void fillReducedCost(const OsiSolverInterface &si, double rCost[]);
};

#endif

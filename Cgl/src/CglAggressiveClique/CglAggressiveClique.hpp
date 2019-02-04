#ifndef _CglAggressiveClique_h_
#define _CglAggressiveClique_h_

#include "CglCutGenerator.hpp"

class CglAggressiveClique : public CglCutGenerator {
public:
  static int sepCuts;
  static double sepTime;

  CglAggressiveClique();
  CglAggressiveClique(const CglAggressiveClique &rhs);
  virtual CglCutGenerator *clone() const;
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info = CglTreeInfo());
  virtual ~CglAggressiveClique();

  void setMaxItBK(int _maxItBK);
  int getMaxItBK() { return maxItBK; }
  void setMaxItBKExt(int _maxItBKExt);
  int getMaxItBKExt() { return maxItBKExt; }
  void setExtendingMethod(int _extMethod);
  int getExtendingMethod() { return extMethod; }

private:
  int maxItBK, maxItBKExt, extMethod;
};

#endif // CGLAGGRESSIVECLIQUE_HPP

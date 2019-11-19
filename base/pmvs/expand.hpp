#ifndef PMVS3_EXPAND_H
#define PMVS3_EXPAND_H

#include <vector>
#include <queue>
#include <list>
#include "patchOrganizerS.hpp"

namespace PMVS3 {
class CFindMatch;
  
class CExpand {
 public:
  CExpand(CFindMatch& findMatch);
  ~CExpand() {};

  void init(void);
  void run(void);

  float computeRadius(const Patch::CPatch& patch);
  
 protected:
  int expandSub(const Patch::PPatch& orgppatch, const int id,
                const Vec4f& canCoord);
  
  int updateCounts(const Patch::CPatch& patch);
  
  int checkCounts(Patch::CPatch& patch);
  
  void findEmptyBlocks(const Patch::PPatch& ppatch,
                       std::vector<std::vector<Vec4f> >& canCoords);
 protected:

  std::priority_queue<Patch::PPatch, std::vector<Patch::PPatch>, P_compare>
    _queue;
  
  CFindMatch& _fm;
  
  //-----------------------------------------------------------------
  // thread related
  //-----------------------------------------------------------------  
  void expandThread(void);
  static int expandThreadTmp(void* arg);

  // Number of trials
  std::vector<int> _ecounts;
  // Number of failures in the prep
  std::vector<int> _fcounts0;
  // Number of failures in the post processing
  std::vector<int> _fcounts1;
  // Number passes
  std::vector<int> _pcounts;
};
};

#endif // PMVS3_EXPAND_H

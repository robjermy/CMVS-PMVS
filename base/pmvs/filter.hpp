#ifndef PMVS3_FILTER_H
#define PMVS3_FILTER_H

#include "patch.hpp"
#include <list>
#include "../numeric/vec2.hpp"

namespace PMVS3 {
  
class CFindMatch;
  
class CFilter {
 public:
  CFilter(CFindMatch& findMatch);

  void init(void);
  void run(void);

  float computeGain(const Patch::CPatch& patch, const int lock);

  int filterQuad(const Patch::CPatch& patch,
                 const std::vector<Patch::PPatch>& neighbors) const;
  
  
 protected:
  void filterOutside(void);
  void filterOutsideThread(void);
  static int filterOutsideThreadTmp(void* arg);

  void filterExact(void);
  void filterExactThread(void);
  static int filterExactThreadTmp(void* arg);
  
  void filterNeighbor(const int time);
  void filterSmallGroups(void);
  void filterSmallGroupsSub(const int pid, const int id,
                            std::vector<int>& label,
                            std::list<int>& ltmp) const;
  void setDepthMaps(void);
  void setDepthMapsVGridsVPGridsAddPatchV(const int additive);
  
  void setConf(const int image);

  std::vector<float> _gains;

  std::vector<std::vector<int> > _newimages, _removeimages;
  std::vector<std::vector<TVec2<int> > > _newgrids, _removegrids;

  int _time;
  std::vector<int> _rejects;
  
  //----------------------------------------------------------------------
  // Thread related
  //----------------------------------------------------------------------
  void setDepthMapsThread(void);
  static int setDepthMapsThreadTmp(void* arg);
  
  void addPatchVThread(void);
  static int addPatchVThreadTmp(void* arg);
  
  void setVGridsVPGridsThread(void);
  static int setVGridsVPGridsThreadTmp(void* arg);

  void filterNeighborThread(void);
  static int filterNeighborThreadTmp(void* arg);
  
  CFindMatch& _fm;
  
};
};

#endif // PMVS3_FILTER_H

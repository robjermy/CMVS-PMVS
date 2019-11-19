#pragma once

#include "patch.hpp"
#include <list>
#include "../numeric/vec2.hpp"

namespace PMVS3 {
  class CFindMatch;

  class CFilter {
  public:
    CFilter(CFindMatch& findMatch);

    void init();
    void run();

    float computeGain(const Patch::CPatch& patch, const int lock);
    int filterQuad(const Patch::CPatch& patch, const std::vector<Patch::PPatch>& neighbors) const;

  protected:
    void filterOutside();
    void filterOutsideThread();
    static int filterOutsideThreadTmp(void* arg);

    void filterExact();
    void filterExactThread();
    static int filterExactThreadTmp(void* arg);

    void filterNeighbor(const int time);
    void filterSmallGroups();
    void filterSmallGroupsSub(const int pid, const int id, std::vector<int>& label, std::list<int>& ltmp) const;
    void setDepthMaps();
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
    void setDepthMapsThread();
    static int setDepthMapsThreadTmp(void* arg);

    void addPatchVThread();
    static int addPatchVThreadTmp(void* arg);

    void setVGridsVPGridsThread();
    static int setVGridsVPGridsThreadTmp(void* arg);

    void filterNeighborThread();
    static int filterNeighborThreadTmp(void* arg);

    CFindMatch& _fm;
  };
}

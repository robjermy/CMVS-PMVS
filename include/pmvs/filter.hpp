#pragma once

#include "patch.hpp"
#include <list>
#include <thread>
#include "numeric/vec2.hpp"

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

    void filterExact();
    void filterExactThread();

    void filterNeighbor(const int time);
    void filterSmallGroups();
    void filterSmallGroupsSub(const int pid, const int id, std::vector<int>& label, std::list<int>& ltmp) const;
    void setDepthMaps();
    void setDepthMapsVGridsVPGridsAddPatchV(const int additive);

    void setConf(const int image);

    std::vector<float> _gains;

    std::vector<std::vector<int>> _newimages, _removeimages;
    std::vector<std::vector<TVec2<int>>> _newgrids, _removegrids;

    int _time;
    std::vector<int> _rejects;

    //----------------------------------------------------------------------
    // Thread related
    //----------------------------------------------------------------------
    void setDepthMapsThread();
    void addPatchVThread();
    void setVGridsVPGridsThread();
    void filterNeighborThread();

    CFindMatch& _fm;
  };
}

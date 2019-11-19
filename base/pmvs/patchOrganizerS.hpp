#ifndef PMVS3_PATCHORGANIZERS_H
#define PMVS3_PATCHORGANIZERS_H

#include "patch.hpp"
#include <queue>

namespace PMVS3 {

class CFindMatch;

class P_compare {
public:
  bool operator()(const Patch::PPatch& lhs, const Patch::PPatch& rhs) const {
    return lhs->m_tmp < rhs->m_tmp;
  }
};
 
class CPatchOrganizerS {
 public:
  CPatchOrganizerS(CFindMatch& findMatch);

  void init(void);
  void collectPatches(const int target = 0);
  void collectPatches(std::priority_queue<Patch::PPatch, std::vector<Patch::PPatch>,
                      P_compare>& pqpatches);
  
  void collectPatches(const int index,
                      std::priority_queue<Patch::PPatch, std::vector<Patch::PPatch>,
                      P_compare>& pqpatches);
  void collectNonFixPatches(const int index, std::vector<Patch::PPatch>& ppatches);
  
  void writePatches2(const std::string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet);
  
  void writePLY(const std::vector<Patch::PPatch>& patches,
                const std::string filename);
  void writePLY(const std::vector<Patch::PPatch>& patches,
                const std::string filename,
                const std::vector<Vec3i>& colors);
  
  void readPatches(void);
  
  void clearCounts(void);
  void clearFlags(void);

  void setGridsImages(Patch::CPatch& patch,
                      const std::vector<int>& images) const;
  void addPatch(Patch::PPatch& ppatch);
  void removePatch(const Patch::PPatch& ppatch);
  void setGrids(Patch::PPatch& ppatch) const;
  void setGrids(Patch::CPatch& patch) const;
  void setVImagesVGrids(Patch::PPatch& ppatch);
  void setVImagesVGrids(Patch::CPatch& patch);
  void updateDepthMaps(Patch::PPatch& ppatch);
  
  int isVisible(const Patch::CPatch& patch, const int image,
                const int& ix, const int& iy,
                const float strict, const int lock);
  int isVisible0(const Patch::CPatch& patch, const int image,
                 int& ix, int& iy,
                 const float strict, const int lock);

  void findNeighbors(const Patch::CPatch& patch,
                     std::vector<Patch::PPatch>& neighbors,
                     const int lock,
                     const float scale = 1.0f,
                     const int margin = 1,
                     const int skipvis = 0);
  
  void setScales(Patch::CPatch& patch) const;
  
  float computeUnit(const Patch::CPatch& patch) const;

  // change the contents of m_images from images to indexes
  void image2index(Patch::CPatch& patch);
  // change the contents of m_images from indexes to images
  void index2image(Patch::CPatch& patch);
  
  //----------------------------------------------------------------------
  // Widths of grids
  std::vector<int> m_gwidths;
  std::vector<int> m_gheights;
  //----------------------------------------------------------------------
  // image, grid
  std::vector<std::vector<std::vector<Patch::PPatch> > > m_pgrids;  
  // image, grid
  std::vector<std::vector<std::vector<Patch::PPatch> > > m_vpgrids;
  // Closest patch
  std::vector<std::vector<Patch::PPatch> > m_dpgrids;

  // all the patches in the current level of m_pgrids 
  std::vector<Patch::PPatch> m_ppatches;

  // Check how many times patch optimization was performed for expansion
  std::vector<std::vector<unsigned char> > m_counts;

  static Patch::PPatch m_MAXDEPTH;
  static Patch::PPatch m_BACKGROUND;
  
 protected:
  CFindMatch& m_fm;
};
};

#endif //PMVS3_PATCHORGANIZERS_H

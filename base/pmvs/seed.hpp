#ifndef PMVS3_SEED_H
#define PMVS3_SEED_H

#include <boost/shared_ptr.hpp>
#include <vector>
#include "patch.hpp"
#include "point.hpp"

#include "tinycthread.h"

namespace PMVS3 {
class CFindMatch;
typedef boost::shared_ptr<CPoint> PPoint;
  
class CSeed {
 public:
  CSeed(CFindMatch& findMatch);
  virtual ~CSeed() {};

  void init(const std::vector<std::vector<CPoint> >& points);
  void run(void);
  void clear(void);

 protected:
  void readPoints(const std::vector<std::vector<CPoint> >& points);
  int canAdd(const int index, const int x, const int y);  

  void initialMatch(const int index, const int id);
  void collectCells(const int index0, const int index1,
                    const CPoint& p0, std::vector<Vec2i>& cells);
  
  void collectCandidates(const int index, const std::vector<int>& indexes,
                         const CPoint& point, std::vector<PPoint>& vcp);

  int initialMatchSub(const int index0, const int index1,
                      const int id, Patch::CPatch& patch);
  
  void unproject(const int index0, const int index1,
                 const CPoint& p0, const CPoint& p1,
                 Vec4f& coord) const;

  //----------------------------------------------------------------------
  CFindMatch& m_fm;
  // points in a grid. For each index, grid
  std::vector<std::vector<std::vector<PPoint> > > m_ppoints;

  //----------------------------------------------------------------------
  // thread related
  //----------------------------------------------------------------------
  void initialMatchThread(void);
  static int initialMatchThreadTmp(void* arg);

  // Number of trials
  std::vector<int> m_scounts;
  // Number of failures in the prep
  std::vector<int> m_fcounts0;
  // Number of failures in the post processing
  std::vector<int> m_fcounts1;
  // Number passes
  std::vector<int> m_pcounts;
};
};

#endif // PMVS3_SEED_H

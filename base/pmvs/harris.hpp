#ifndef PMVS3_HARRIS_H
#define PMVS3_HARRIS_H

#include <vector>
#include <set>
#include "../numeric/vec3.hpp"
#include "detector.hpp"
#include "point.hpp"

namespace PMVS3 {

class CHarris: public CDetector {
 public:
  void run(const std::vector<unsigned char>& image,
	   const std::vector<unsigned char>& mask,
	   const std::vector<unsigned char>& edge,           
	   const int width, const int height,	   
	   const int gspeedup, const float sigma,
	   std::multiset<CPoint> & result);

  virtual ~CHarris() {
  }
  
 protected:
  float m_sigmaD;
  float m_sigmaI;
  
  std::vector<float> m_gaussD;
  std::vector<float> m_gaussI;
  
  std::vector<std::vector<Vec3f> > m_dIdx;
  std::vector<std::vector<Vec3f> > m_dIdy;
  
  std::vector<std::vector<float> > m_dIdxdIdx;
  std::vector<std::vector<float> > m_dIdydIdy;
  std::vector<std::vector<float> > m_dIdxdIdy;

  std::vector<std::vector<float> > m_response;
  
  void init(const std::vector<unsigned char>& image,
            const std::vector<unsigned char>& mask,
	    const std::vector<unsigned char>& edge);
  
  void setDerivatives(void);
  void preprocess(void);
  void preprocess2(void);
  void setResponse(void);
};
};

#endif // PMVS3_HARRIS_H

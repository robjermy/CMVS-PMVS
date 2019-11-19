#include <iostream>
#include <fstream>
#include "../image/image.hpp"
#include "detectFeatures.hpp"
#include "harris.hpp"
#include "dog.hpp"
#include "point.hpp"

using namespace PMVS3;
using namespace std;
using namespace Image;

CDetectFeatures::CDetectFeatures(void) {
  mtx_init(&m_rwlock, mtx_plain | mtx_recursive);
}

CDetectFeatures::~CDetectFeatures() {
  mtx_destroy(&m_rwlock);
}

void CDetectFeatures::run(const CPhotoSetS& pss, const int num,
                          const int csize, const int level,
                          const int CPU) {
  m_ppss = &pss;
  m_csize = csize;
  m_level = level;
  m_CPU = CPU;

  m_points.clear();
  m_points.resize(num);
  
  //----------------------------------------------------------------------
  for (int index = 0; index < num; ++index)
    m_jobs.push_back(index);
  
  vector<thrd_t> threads(m_CPU);
  for (int i = 0; i < m_CPU; ++i)
    thrd_create(&threads[i], &runThreadTmp, (void*)this);
  for (int i = 0; i < m_CPU; ++i)
    thrd_join(threads[i], NULL);
  //----------------------------------------------------------------------
  cerr << "done" << endl;
}

int CDetectFeatures::runThreadTmp(void* arg) {
  CDetectFeatures* detectFeatures = (CDetectFeatures*)arg;  
  detectFeatures->runThread();
  return 0;
}

void CDetectFeatures::runThread(void) {
  while (1) {
    int index = -1;
    mtx_lock(&m_rwlock);
    if (!m_jobs.empty()) {
      index = m_jobs.front();
      m_jobs.pop_front();
    }
    mtx_unlock(&m_rwlock);
    if (index == -1)
      break;
    
    const int image = m_ppss->m_images[index];
    cerr << image << ' ' << flush;

    //?????????????  May need file lock, because targetting images
    //should not overlap among multiple processors.    
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.affin%d", m_ppss->m_prefix.c_str(), image, m_level);
    ifstream ifstr;
    ifstr.open(buffer);
    if (ifstr.is_open()) {
      ifstr.close();
      continue;
    }
    ifstr.close();
    
    //----------------------------------------------------------------------
    // parameters
    // for harris
    const float sigma = 4.0f;
    // for DoG
    const float firstScale = 1.0f;    const float lastScale = 3.0f;

    //----------------------------------------------------------------------
    // Harris
    {
      CHarris harris;
      multiset<CPoint> result;
      harris.run(m_ppss->m_photos[index].getImage(m_level),
                 m_ppss->m_photos[index].CImage::getMask(m_level),
                 m_ppss->m_photos[index].CImage::getEdge(m_level),
                 m_ppss->m_photos[index].getWidth(m_level),
                 m_ppss->m_photos[index].getHeight(m_level), m_csize, sigma, result);
      
      multiset<CPoint>::reverse_iterator rbegin = result.rbegin();
      while (rbegin != result.rend()) {
        m_points[index].push_back(*rbegin);
        rbegin++;
      }
    }

    //----------------------------------------------------------------------
    // DoG
    {
      CDifferenceOfGaussians dog;
      multiset<CPoint> result;
      dog.run(m_ppss->m_photos[index].getImage(m_level),
              m_ppss->m_photos[index].CImage::getMask(m_level),
              m_ppss->m_photos[index].CImage::getEdge(m_level),
              m_ppss->m_photos[index].getWidth(m_level),
              m_ppss->m_photos[index].getHeight(m_level),
              m_csize, firstScale, lastScale, result);
      
      multiset<CPoint>::reverse_iterator rbegin = result.rbegin();      
      while (rbegin != result.rend()) {
        m_points[index].push_back(*rbegin);
        rbegin++;
      }
    }
  }
}

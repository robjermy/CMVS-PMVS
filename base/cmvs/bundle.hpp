#pragma once

#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <ctime>
#include <time.h> // PM

#include "disjoint.hpp"
#include "../stann/sfcnn.hpp"
#include "../numeric/mat3.hpp"
#include "../image/photoSetS.hpp"

namespace CMVS {
  struct Sadd {
    Sadd(const int image, const float gain) : _image(image), _gain(gain) {};

    int _image;
    float _gain;
  };

  struct Ssfm2 {
    Ssfm2(void) : _cluster(-1), _score(-2.0f), _scoreThreshold(-1.0f), _satisfied(1) {}
    // which cluster it belongs to currently
    int _cluster;

    // current score
    float _score;

    // score threshold
    float _scoreThreshold;

    // If SFM is satisfied or not.
    // 1: satisfied, 0: not satisfied
    //
    // In adding images,
    // 2: currently not satisfied, 1: satisfied,
    // 0: not satisfied from the begining and no hope
    char _satisfied;

    // For SfM point that has not bee satisfied, compute several number
    // of images that can be added to gain more info
    std::vector<Sadd> _adds;

    // best images
    std::vector<int> _uimages;
  };

  class CBundle {
  public:
    CBundle();
    virtual ~CBundle();

    void run(const std::string prefix, const int imageThreshold, const int tau, const float scoreRatioThreshold, const float coverageThreshold, const int pnumThreshold, const int CPU);
    // root dir
    std::string _prefix;

    // # of cameras
    int _cnum;
    // # of points
    int _pnum;

    // Point params
    std::vector<Vec4f> _coords;
    std::vector<std::vector<int> > _visibles;

    std::vector<Vec3f> _colors;

    // A set of point ids visible in each camera
    std::vector<std::vector<int> > _vpoints;

    std::vector<int> _pweights;

    //----------------------------------------------------------------------
    // Generated data
    //----------------------------------------------------------------------
    // A list of connected images, for each camera.
    std::vector<std::vector<int> > _neighbors;

    // Width and height of depth map
    std::vector<int> _widths;
    std::vector<int> _heights;
    // scale
    std::vector<int> _levels;

    //----------------------------------------------------------------------
    // Output
    //----------------------------------------------------------------------
    // clusters
    // _timages need to be sorted for addImages.
    std::vector<std::vector<int> > _timages;
    std::vector<std::vector<int> > _oimages;

  protected:
    void prep(const std::string prefix, const int imageThreshold, const int tau, const float scoreRatioThreshold, const float coverageThreshold, const int pnumThreshold, const int CPU);

    void prep2(void);

    void readBundle(const std::string file);
    void setWidthsHeightsLevels(void);
    void setNeighbors(void);

    int totalNum(void) const;

    // set _scoreThresholds
    void setScoreThresholds(void);

    void resetVisibles(void);

    // set new images without image while taking into account _removed
    void setNewImages(const int pid, const int rimage, std::vector<int>& newimages);

    void sRemoveImages(void);
    void checkImage(const int image);

    void setCluster(const int p);

    void setScoresClusters(void);

    // For unsatisfied sfm points, update cluster
    void setClusters(void);

    void slimNeighborsSetLinks(void);

    float computeLink(const int image0, const int image1);

    void addImagesP(void);
    int addImages(void);
    int addImagesSub(const std::vector<std::map<int, float> >& cands);

    // angle score
    static float angleScore(const Vec4f& ray0, const Vec4f& ray1);

    void mergeSfM(void);
    void mergeSfMP(void);
    void mergeSfMPThread(void);

    std::vector<char> _merged;

    void findPNeighbors(sfcnn<const float*, 3, float>& tree, const int pid, std::vector<int>& pneighbors);

    void resetPoints(void);

    static void mymerge(const std::vector<int>& lhs, const std::vector<int>& rhs, std::vector<int>& output);

    static int my_isIntersect(const std::vector<int>& lhs, const std::vector<int>& rhs);

    // Cluster images
    void setTimages(void);
    void divideImages(const std::vector<int>& lhs, std::vector<std::vector<int> >& rhs);

    float computeScore2(const Vec4f& coord, const std::vector<int>& images) const;
    float computeScore2(const Vec4f& coord, const std::vector<int>& images, std::vector<int>& uimages) const;
    // Enforce the specified image to be inside
    float computeScore2(const Vec4f& coord, const std::vector<int>& images, const int index) const;

    void writeCameraCenters(void);
    void writeVis(void);
    void writeGroups(void);
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Link info
    std::vector<std::vector<float> > _links;

    // Removed or not for images
    std::vector<int> _removed;
    // Number of SFM points above threshold for each image. We can
    // remove images until _allows is non-negative for all the images.
    // Used in removing images
    std::vector<int> _allows;
    // Used in adding images
    std::vector<int> _lacks;

    // For an image, how sfm point (_vpoints) changes if the image is
    // removed.
    //  0: unsatisfy
    //  1: satisfy->satisfy
    //  2: satisfy->unsatisfy
    std::vector<char> _statsT;
    // image under consideration
    int _imageT;
    // The value of lacks
    int _lacksT;

    // sfm information used in addimages
    std::vector<Ssfm2> _sfms2;

    // Number of images used in computeScore2
    int _tau;
    // union find operations to be executed
    std::vector<std::vector<std::vector<int> > > _ufsT;
    // Smallest scale
    std::vector<float> _minScales;

    // add nums
    std::vector<int> _addnums;

    // scaling factor for depth
    float _dscale;
    // scaling for kdtree version
    float _dscale2;

    //----------------------------------------------------------------------
    Image::CPhotoSetS _pss;

    // depth level
    int _dlevel;
    // maxLevel in _pss.
    int _maxLevel;

    int _imageThreshold;
    // Num of points for images to be connected
    int _pnumThreshold;

    // link threshold for neighbor
    float _linkThreshold;

    // Score ratio threshold. Optimal score using all the visible images
    // times this threshold is the mimimum possible score to be
    // satisfied.
    float _scoreRatioThreshold;
    // How much SFM must be satisfied in each image.
    float _coverageThreshold;

    // union find for sfm points
    DisjointSetForest<int>* _puf2;

    sfcnn<const float*, 3, float>* _ptree;

    //----------------------------------------------------------------------
    // Threads
    int _CPU;
    mtx_t _lock;
    std::list<int> _jobs;
    int _junit;
    int _thread;
    int _count;

    int _debug;

    void startTimer(void);
    time_t curTimer(void);

    time_t _tv; //PM
    time_t _curtime;
  };
}

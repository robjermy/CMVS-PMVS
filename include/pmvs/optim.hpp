#pragma once

#include <vector>

#include "patch.hpp"

namespace PMVS3 {
  class CFindMatch;

  class COptim {
  public:
    COptim(CFindMatch& findMatch);

    void init(void);

    //-----------------------------------------------------------------
    // Image manipulation
    //-----------------------------------------------------------------
    void collectImages(const int index, std::vector<int>& indexes) const;
    void addImages(Patch::CPatch& patch) const;
    void removeImagesEdge(Patch::CPatch& patch) const;

    float getUnit(const int index, const Vec4f& coord) const;

    void computeUnits(const Patch::CPatch& patch, std::vector<int>& indexes, std::vector<float>& fineness, std::vector<Vec4f>& rays) const;
    void computeUnits(const Patch::CPatch& patch, std::vector<float>& fineness) const;

    //-----------------------------------------------------------------
    // Optimization
    //-----------------------------------------------------------------
    int preProcess(Patch::CPatch& patch, const int id, const int seed);
    void refinePatch(Patch::CPatch& patch, const int id, const int time);

    bool refinePatchBFGS(Patch::CPatch& patch, const int id, const int time, const int ncc);

    int postProcess(Patch::CPatch& patch, const int id, const int seed);

    void setRefImage(Patch::CPatch& patch, const int id);

    int check(Patch::CPatch& patch);

    std::vector<int> _status;

  protected:
    void filterImagesByAngle(Patch::CPatch& patch);

    void sortImages(Patch::CPatch& patch) const;
    void constraintImages(Patch::CPatch& patch, const float nccThreshold, const int id);
    void setRefConstraintImages(Patch::CPatch& patch, const float nccThreshold, const int id);

    void setINCCs(const Patch::CPatch& patch, std::vector<float> & nccs, const std::vector<int>& indexes, const int id, const int robust);

    void setINCCs(const Patch::CPatch& patch, std::vector<std::vector<float> >& nccs, const std::vector<int>& indexes, const int id, const int robust);

    int grabTex(const Vec4f& coord, const Vec4f& pxaxis, const Vec4f& pyaxis, const Vec4f& pzaxis, const int index, const int size, std::vector<float>& tex) const;

    int grabSafe(const int index, const int size, const Vec3f& center, const Vec3f& dx, const Vec3f& dy, const int level) const;

    double computeINCC(const Vec4f& coord, const Vec4f& normal, const std::vector<int>& indexes, const Vec4f& pxaxis, const Vec4f& pyaxis, const int id, const int robust);

  public:
    static void normalize(std::vector<float>& tex);
    static void normalize(std::vector<std::vector<float> >& texs, const int size);

    float dot(const std::vector<float>& tex0, const std::vector<float>& tex1) const;
    float ssd(const std::vector<float>& tex0, const std::vector<float>& tex1) const;

  protected:
    static void lfunc(double* p, double* hx, int m, int n, void* adata);
    void func(int m, int n, double* x, double* fvec, int* iflag, void* arg);

    //BFGS
    static double my_f(unsigned n, const double *x, double *grad, void *my_func_data);

    void encode(const Vec4f& coord, double* const vect, const int id) const;
    void encode(const Vec4f& coord, const Vec4f& normal, double* const vect, const int id) const;
    void decode(Vec4f& coord, Vec4f& normal, const double* const vect, const int id) const;
    void decode(Vec4f& coord, const double* const vect, const int id) const;

  public:
    void setWeightsT(const Patch::CPatch& patch, const int id);

    double computeINCC(const Vec4f& coord, const Vec4f& normal, const std::vector<int>& indexes, const int id, const int robust);
    void getPAxes(const int index, const Vec4f& coord, const Vec4f& normal, Vec4f& pxaxis, Vec4f& pyaxis) const;

    static inline float robustincc(const float rhs) {
      return rhs / (1 + 3 * rhs);
    }

    static inline float unrobustincc(const float rhs) {
      return rhs / (1 - 3 * rhs);
    }

  protected:
    void setAxesScales(void);

    static COptim* _one;
    CFindMatch& _fm;

    //-----------------------------------------------------------------
    // Axes
    std::vector<Vec3f> _xaxes;
    std::vector<Vec3f> _yaxes;
    std::vector<Vec3f> _zaxes;
    // Scales
    std::vector<float> _ipscales;

    //-----------------------------------------------------------------
    // For threads
    std::vector<float> _vect0T;
    std::vector<Vec4f> _centersT;
    std::vector<Vec4f> _raysT;
    std::vector<std::vector<int> > _indexesT;
    std::vector<float> _dscalesT;
    std::vector<float> _ascalesT;

    // stores current parameters for derivative computation
    std::vector<Vec3f> _paramsT;

    // Grabbed texture
    std::vector<std::vector<std::vector<float> > > _texsT; // last is 7x7x3 patch
    // weights for refineDepthOrientationWeighed
    std::vector<std::vector<float> > _weightsT;
    // Working array for levmar
    std::vector<std::vector<double> > _worksT;
  };
}

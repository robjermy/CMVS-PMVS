#pragma once

#include <iomanip>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace PMVS3 {
  struct SOption{
    public:
    int _level;
    int _csize;
    float _threshold;
    int _wsize;
    int _minImageNum;
    int _CPU;
    float _setEdge;
    int _useBound;
    int _useVisData;
    int _sequence;

    float _maxAngleThreshold;
    float _quadThreshold;

    std::string _prefix;
    std::string _option;
    int _tflag;
    std::vector<int> _timages;
    int _oflag;
    std::vector<int> _oimages;

    std::map<int, int> _dict;

    std::vector<int> _bindexes;
    std::vector<std::vector<int> > _visdata;
    std::vector<std::vector<int> > _visdata2;

    SOption(const std::string& config);
    SOption(const std::string& prefix, const std::string& option);

    void init(const std::string& prefix, const std::string& option);

    friend std::ostream& operator<<(std::ostream& ostr, const SOption& options);
    protected:
    SOption();

    void initOimages(void);
    void initBindexes(const std::string& sbimages);
    void initVisdata(void);
    void initVisdata2(void);

    void getImagesFromConfig(const nlohmann::json& useImages, std::vector<int>& images);
  };

  inline std::ostream &operator<<(std::ostream& out, const SOption& options) {
    const int ow = 20;
    return out << "Options" << std::endl
      << std::setw(ow) << "level: " << options._level << std::endl
      << std::setw(ow) << "cellSize: " << options._csize << std::endl
      << std::setw(ow) << "windowSize: " << options._wsize << std::endl
      << std::setw(ow) << "minImages: " << options._minImageNum << std::endl
      << std::setw(ow) << "threshold: " << options._threshold << std::endl
      << std::setw(ow) << "threads: " << options._CPU << std::endl
      << std::setw(ow) << "setEdge: " << options._setEdge << std::endl
      << std::setw(ow) << "useBound: " << options._useBound << std::endl
      << std::setw(ow) << "useVisData: " << options._useVisData << std::endl
      << std::setw(ow) << "maxAngleThreshold: " << options._maxAngleThreshold << std::endl
      << std::setw(ow) << "quadThreshold: " << options._quadThreshold << std::endl;
  }
}

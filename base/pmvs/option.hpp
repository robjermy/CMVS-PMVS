#ifndef PMVS3_OPTION_H
#define PMVS3_OPTION_H

#include <string>
#include <vector>
#include <map>

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
  
  SOption(void);
  
  void init(const std::string prefix, const std::string option);
  
  protected:
  void initOimages(void);
  void initBindexes(const std::string sbimages);
  void initVisdata(void);
  void initVisdata2(void);
};
};

#endif // PMVS3_OPTION_H

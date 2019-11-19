#include <algorithm>
#include "harris.hpp"

using namespace PMVS3;
using namespace std;

void CHarris::init(const std::vector<unsigned char>& image,
		   const std::vector<unsigned char>& mask,
                   const std::vector<unsigned char>& edge) {
  _image.clear();
  _image.resize(_height);
  int count = 0;
  for (int y = 0; y < _height; ++y) {
    _image[y].resize(_width);
    for (int x = 0; x < _width; ++x) {
      _image[y][x][0] = ((int)image[count++]) / 255.0f;
      _image[y][x][1] = ((int)image[count++]) / 255.0f;
      _image[y][x][2] = ((int)image[count++]) / 255.0f;
    }
  }

  _mask.clear();
  if (!mask.empty() || !edge.empty()) {
    _mask.resize(_height);
    count = 0;
    for (int y = 0; y < _height; ++y) {
      _mask[y].resize(_width);
      for (int x = 0; x < _width; ++x) {
        if (mask.empty())
          _mask[y][x] = edge[count++];
        else if (edge.empty())
          _mask[y][x] = mask[count++];
        else {
          if (mask[count] && edge[count])
            _mask[y][x] = (unsigned char)255;
          else
            _mask[y][x] = 0;
          count++;
        }
      }
    }
  }
  
  setGaussD(_sigmaD, _gaussD);
  setGaussI(_sigmaI, _gaussI); 
}

void CHarris::setDerivatives(void) {
  // set _dIdx, _dIdy
  preprocess();
  
  // now set _dIdxdIdx, _dIdydIdy, _dIdxDIdy
  preprocess2();
}

void CHarris::preprocess2(void) {
  _dIdxdIdx.clear();  _dIdydIdy.clear();  _dIdxdIdy.clear();
  _dIdxdIdx.resize(_height);
  _dIdydIdy.resize(_height);
  _dIdxdIdy.resize(_height);
  for (int y = 0; y < _height; ++y) {
    _dIdxdIdx[y].resize(_width);
    _dIdydIdy[y].resize(_width);
    _dIdxdIdy[y].resize(_width);
    for (int x = 0; x < _width; ++x) {
      _dIdxdIdx[y][x] = _dIdydIdy[y][x] = _dIdxdIdy[y][x] = 0.0;
      if (!_mask.empty() && _mask[y][x] == 0)	continue;
      
      _dIdxdIdx[y][x] += _dIdx[y][x] * _dIdx[y][x];
      _dIdydIdy[y][x] += _dIdy[y][x] * _dIdy[y][x];
      _dIdxdIdy[y][x] += _dIdx[y][x] * _dIdy[y][x];
    }
  }

  {
    vector<vector<Vec3f> >().swap(_dIdx);
    vector<vector<Vec3f> >().swap(_dIdy);
  }
  
  //----------------------------------------------------------------------
  // blur
  vector<vector<float> > vvftmp;
  vvftmp.resize(_height);
  for (int y = 0; y < _height; ++y) {
    vvftmp[y].resize(_width);
    for (int x = 0; x < _width; ++x)
      vvftmp[y][x] = 0.0;
  }
  
  //----------------------------------------------------------------------
  // _dIdxdIdx
  convolveX(_dIdxdIdx, _mask, _gaussI, vvftmp);
  convolveY(_dIdxdIdx, _mask, _gaussI, vvftmp);
  
  //----------------------------------------------------------------------
  // _dIdydIdy
  convolveX(_dIdydIdy, _mask, _gaussI, vvftmp);
  convolveY(_dIdydIdy, _mask, _gaussI, vvftmp);
  
  //----------------------------------------------------------------------
  // _dIdxdIdy
  convolveX(_dIdxdIdy, _mask, _gaussI, vvftmp);
  convolveY(_dIdxdIdy, _mask, _gaussI, vvftmp);
}

void CHarris::preprocess(void) {
  vector<vector<Vec3f> > vvvftmp;
  vvvftmp.resize(_height);
  for (int y = 0; y < _height; ++y) {
    vvvftmp[y].resize(_width);
    for (int x = 0; x < _width; ++x)
      vvvftmp[y][x] = Vec3f();
  }
  
  _dIdx = _image;

  vector<float> dfilter, ifilter;
  dfilter.resize(3);
  dfilter[0] = -0.5;  dfilter[1] = 0;  dfilter[2] = 0.5;
  ifilter.resize(3);
  ifilter[0] = 1.0 / 3.0;  ifilter[1] = 1.0 / 3.0;  ifilter[2] = 1.0 / 3.0;

  convolveX(_dIdx, _mask, dfilter, vvvftmp);
  convolveY(_dIdx, _mask, ifilter, vvvftmp);
  
  _dIdy = _image;
  convolveX(_dIdy, _mask, ifilter, vvvftmp);
  convolveY(_dIdy, _mask, dfilter, vvvftmp);
}

void CHarris::setResponse(void) {
  _response.clear();
  _response.resize(_height);
  for (int y = 0; y < _height; ++y) {
    _response[y].resize(_width);
    for (int x = 0; x < _width; ++x) {
      _response[y][x] = 0.0;
      if (!_mask.empty() && _mask[y][x] == 0)	continue;
      
      const float D = _dIdxdIdx[y][x] * _dIdydIdy[y][x] - _dIdxdIdy[y][x] * _dIdxdIdy[y][x];
      const float tr = _dIdxdIdx[y][x] + _dIdydIdy[y][x];
      _response[y][x] = D - 0.06 * tr * tr;
    }
  }
  
  //----------------------------------------------------------------------
  // suppress non local max
  vector<vector<float> > vvftmp = _response;
  for (int y = 1; y < _height - 1; ++y) {
    for (int x = 1; x < _width - 1; ++x) {
      if (_response[y][x] < _response[y][x+1] ||
	  _response[y][x] < _response[y][x-1] ||
	  _response[y][x] < _response[y+1][x] ||
	  _response[y][x] < _response[y-1][x])
	vvftmp[y][x] = 0.0;
    }
  }
  
  vvftmp.swap(_response);
}

void CHarris::run(const std::vector<unsigned char>& image,
                  const std::vector<unsigned char>& mask,
                  const std::vector<unsigned char>& edge,
		  const int width, const int height,
		  const int gspeedup, const float sigma,
		  std::multiset<CPoint> & result) {

  cerr << "Harris running ..." << flush;
  _width = width;       _height = height;
  _sigmaD = sigma;      _sigmaI = sigma;
  init(image, mask, edge);
  setDerivatives();  setResponse();

  const int factor = 2;
  const int maxPointsGrid = factor * factor;
  const int gridsize = gspeedup * factor;

  const int w = (_width + gridsize - 1) / gridsize;
  const int h = (_height + gridsize - 1) / gridsize;
  
  vector<vector<multiset<CPoint> > > resultgrids;
  resultgrids.resize(h);
  for (int y = 0; y < h; ++y)
    resultgrids[y].resize(w);

  const int margin = (int)_gaussD.size() / 2;
  for (int y = margin; y < _height - margin; ++y) {
    for (int x = margin; x < _width - margin; ++x) {
      if (_response[y][x] == 0.0)
	continue;

      const int x0 = min(x / gridsize, w - 1);
      const int y0 = min(y / gridsize, h - 1);
      
      if ((int)resultgrids[y0][x0].size() < maxPointsGrid ||
	  resultgrids[y0][x0].begin()->_response < _response[y][x]) {
	CPoint p;
	p._icoord = Vec3f(x, y, 1.0f);
	p._response = _response[y][x];
	p._type = 0;
	
	resultgrids[y0][x0].insert(p);
	if (maxPointsGrid < (int)resultgrids[y0][x0].size())
	  resultgrids[y0][x0].erase(resultgrids[y0][x0].begin ());
      }
    }
  }  
  
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x) {
      //const float threshold = setThreshold(resultgrids[y][x]);      
      multiset<CPoint>::iterator begin = resultgrids[y][x].begin();
      multiset<CPoint>::iterator end = resultgrids[y][x].end();
      while (begin != end) {
	//if (threshold <= begin->_response)
	  result.insert(*begin);
	begin++;
      }
    }

  cerr << (int)result.size() << " harris done" << endl;
}

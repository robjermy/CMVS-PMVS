#include <algorithm>
#include "dog.hpp"
#include "point.hpp"

using namespace PMVS3;
using namespace std;

int CDifferenceOfGaussians::notOnEdge(const std::vector<std::vector<float> >& dog, int x, int y) {
  return 1;
  //const float thresholdEdge = 0.06f;
  const float thresholdEdge = 0.06f;
  
  const float H00 = dog[y][x-1] - 2.0f * dog[y][x] + dog[y][x+1];
  const float H11 = dog[y-1][x] - 2.0f * dog[y][x] + dog[y+1][x];
  const float H01 = ((dog[y+1][x+1] - dog[y-1][x+1]) - (dog[y+1][x-1] - dog[y-1][x-1])) / 4.0f;
  const float det = H00 * H11 - H01 * H01;
  const float trace = H00 + H11;
  return det > thresholdEdge * trace * trace;
}

float CDifferenceOfGaussians::getResponse(const std::vector<std::vector<float> >& pdog,
		  const std::vector<std::vector<float> >& cdog,
		  const std::vector<std::vector<float> >& ndog,
		  const int x, const int y) {
  return fabs(getResponse(pdog, x, y) + getResponse(cdog, x, y) + getResponse(ndog, x, y) +
	      (cdog[y][x] - pdog[y][x]) + (cdog[y][x] - ndog[y][x]));
}

float CDifferenceOfGaussians::getResponse(const std::vector<std::vector<float> >& dog,
		  const int x, const int y) {
  const float sum =
    dog[y-1][x-1] + dog[y  ][x-1] +
    dog[y+1][x-1] + dog[y-1][x  ] +
    dog[y+1][x  ] + dog[y-1][x+1] +
    dog[y  ][x+1] + dog[y+1][x+1];

  return 8 * dog[y][x] - sum;
}  

int CDifferenceOfGaussians::isLocalMax(const std::vector<std::vector<float> >& dog,
		     const int x, const int y) {
  const float value = dog[y][x];

  if (0.0 < value) {
    if (dog[y-1][x-1] < value && dog[y  ][x-1] < value &&
	dog[y+1][x-1] < value && dog[y-1][x  ] < value &&
	dog[y+1][x  ] < value && dog[y-1][x+1] < value &&
	dog[y  ][x+1] < value && dog[y+1][x+1] < value)
      return 1;
  }
  else {
    if (dog[y-1][x-1] > value && dog[y  ][x-1] > value &&
	dog[y+1][x-1] > value && dog[y-1][x  ] > value &&
	dog[y+1][x  ] > value && dog[y-1][x+1] > value &&
	dog[y  ][x+1] > value && dog[y+1][x+1] > value)
      return -1;
  }
  return 0;
}

int CDifferenceOfGaussians::isLocalMax(const std::vector<std::vector<float> >& pdog,
		     const std::vector<std::vector<float> >& cdog,
		     const std::vector<std::vector<float> >& ndog,
		     const int x, const int y) {
  /*
  const int flag0 = isLocalMax(pdog, x, y);
  const int flag1 = isLocalMax(cdog, x, y);
  const int flag2 = isLocalMax(ndog, x, y);
  
  if (flag0 == 1 && flag1 == 1 && flag2 == 1) {
    if (pdog[y][x] < cdog[y][x] && ndog[y][x] < cdog[y][x])
      return 1;
    else
      return 0;
  }
  else if (flag0 == -1 && flag1 == -1 && flag2 == -1) {
    if (cdog[y][x] < pdog[y][x] && cdog[y][x] < ndog[y][x])
      return -1;
    else
      return 0;
  }
  return 0;
  */

  const int flag = isLocalMax(cdog, x, y);
  
  if (flag == 1) {
    if (pdog[y][x] < cdog[y][x] && ndog[y][x] < cdog[y][x])
      return 1;
    else
      return 0;
  }
  else if (flag == -1) {
    if (cdog[y][x] < pdog[y][x] && cdog[y][x] < ndog[y][x])
      return -1;
    else
      return 0;
  }
  return 0;
}

void CDifferenceOfGaussians::setDOG(const std::vector<std::vector<float> >& cres,
		  const std::vector<std::vector<float> >& nres,
		  std::vector<std::vector<float> >& dog) {
  const int height = (int)cres.size();
  const int width = (int)cres[0].size();

  dog = nres;
  for (int y = 0; y < height; ++y)
    for (int x = 0; x < width; ++x)
      dog[y][x] -= cres[y][x];
}

void CDifferenceOfGaussians::run(const std::vector<unsigned char>& image,
               const std::vector<unsigned char>& mask,
               const std::vector<unsigned char>& edge,
	       const int width, const int height,
	       const int gspeedup,
	       const float firstScale,   // 1.4f
	       const float lastScale,    // 4.0f
	       std::multiset<CPoint> & result) {
  cerr << "DoG running..." << flush;
  _width = width;
  _height = height;
  
  _firstScale = firstScale;
  _lastScale = lastScale;
  
  init(image, mask, edge);

  const int factor = 2;
  const int maxPointsGrid = factor * factor;
  const int gridsize = gspeedup * factor;
  
  const int w = (_width + gridsize - 1) / gridsize;
  const int h = (_height + gridsize - 1) / gridsize;
  
  /*
  const int gridsize = 50;
  const int w = (int)ceil(_width / (float)gridsize);
  const int h = (int)ceil(_height / (float)gridsize);
  const int maxPointsGrid = (int)(_width * _height * 0.0025 / (w * h));
  */
  
  vector<vector<multiset<CPoint> > > resultgrids;
  resultgrids.resize(h);
  for (int y = 0; y < h; ++y)
    resultgrids[y].resize(w);

  const float scalestep = pow(2.0f, 1 / 2.0f);
  //const float scalestep = pow(2.0f, 1.0f);
  const int steps = max(4, (int)ceil(log(_lastScale / _firstScale) / log(scalestep)));
  
  vector<vector<float> > pdog, cdog, ndog, cres, nres;

  setRes(_firstScale, cres);
  setRes(_firstScale * scalestep, nres);
  setDOG(cres, nres, cdog);
  cres.swap(nres);
  setRes(_firstScale * scalestep * scalestep, nres);
  setDOG(cres, nres, ndog);

  vector<vector<unsigned char> > alreadydetected;
  alreadydetected.resize(_height);
  for (int y = 0; y < _height; ++y) {
    alreadydetected[y].resize(_width);
    for (int x = 0; x < _width; ++x) {
      alreadydetected[y][x] = (unsigned char)0;
    }
  }

  for (int i = 2; i <= steps - 1; ++i) {
    const float cscale = _firstScale * pow(scalestep, i + 1);
    cres.swap(nres);
    setRes(cscale, nres);
    
    pdog.swap(cdog);
    cdog.swap(ndog);
    setDOG(cres, nres, ndog);
    
    const int margin = (int)ceil(2 * cscale);
    // now 3 response maps are ready
    for (int y = margin; y < _height - margin; ++y) {
      for (int x = margin; x < _width - margin; ++x) {
	if (alreadydetected[y][x])
	  continue;
	if (cdog[y][x] == 0.0)
	  continue;
	//if (isCloseBoundary(x, y, margin))
        //continue;
	
	// check local maximum
	if (isLocalMax(pdog, cdog, ndog, x, y) && notOnEdge(cdog, x, y)) {
	  const int x0 = min(x / gridsize, w - 1);
	  const int y0 = min(y / gridsize, h - 1);

	  alreadydetected[y][x] = 1;
	  CPoint p;
	  p._icoord = Vec3f(x, y, 1.0f);
	  p._response = fabs(cdog[y][x]);
	  p._type = 1;
	  
	  resultgrids[y0][x0].insert(p);

	  if (maxPointsGrid < (int)resultgrids[y0][x0].size())
	    resultgrids[y0][x0].erase(resultgrids[y0][x0].begin());
	}
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
  
  cerr << (int)result.size() << " dog done" << endl;  
}

void CDifferenceOfGaussians::setRes(const float sigma,
		  std::vector<std::vector<float> >& res) {
  vector<float> gauss;
  setGaussI(sigma, gauss);
  
  vector<vector<Vec3f> > vvftmp;
  vvftmp.resize((int)_image.size());
  for (int y = 0; y < (int)_image.size(); ++y)
    vvftmp[y].resize((int)_image[y].size());

  vector<vector<Vec3f> > restmp = _image;
  convolveX(restmp, _mask, gauss, vvftmp);
  convolveY(restmp, _mask, gauss, vvftmp);

  res.resize((int)_image.size());
  for (int y = 0; y < (int)_image.size(); ++y) {
    res[y].resize((int)_image[y].size());
    for (int x = 0; x < (int)_image[y].size(); ++x)
      res[y][x] = norm(restmp[y][x]);
  }
}		  

void CDifferenceOfGaussians::init(const std::vector<unsigned char>& image,
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
}

#include <fstream>
#include <iterator>
#include <numeric>
#include <random>

#include "cmvs/bundle.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


#define _USE_MATH_DEFINES
#include <math.h>

extern "C"
{
  int boundary_points;
  int spectral_initialization = 0;
  int cutType;
  int memory_saving;
};

CMVS::CBundle::CBundle() {
  _CPU = 8;
  _junit = 100;
  _debug = 0;
  // _puf = NULL;
  _puf2 = nullptr;
  _ptree = NULL;
}

CMVS::CBundle::~CBundle() {}

void CMVS::CBundle::prep(const std::string prefix, const int imageThreshold, const int tau, const float scoreRatioThreshold, const float coverageThreshold, const int pnumThreshold, const int CPU) {
  if (pnumThreshold != 0) {
    std::cerr << "Should use pnumThreshold = 0" << std::endl;
    exit (1);
  }

  _prefix = prefix;
  _imageThreshold = imageThreshold;
  _tau = tau;
  _scoreRatioThreshold = scoreRatioThreshold;
  _coverageThreshold = coverageThreshold;
  _pnumThreshold = pnumThreshold;
  _CPU = CPU;

  _linkThreshold = 2.0f;

  _dscale = 1 / 100.0f;
  _dscale2 = 1.0f;

  char buffer[1024];
  sprintf(buffer, "%sbundle.rd.out", prefix.c_str());
  std::cerr << "Reading bundle..." << flush;
  readBundle(buffer);
  std::cerr << std::endl;

  std::vector<int> images;
  for (int c = 0; c < _cnum; ++c) {
    images.push_back(c);
  }

  _dlevel = 7;
  _maxLevel = 12;
  _pss.init(images, prefix, _maxLevel + 1, 5, 0);

  std::cerr << "Set widths/heights..." << flush;
  setWidthsHeightsLevels();
  std::cerr << "done" << flush;
}

void CMVS::CBundle::prep2() {
  { // Used in mergeSfMP now.
    _pweights.resize((int)_coords.size(), 1);
    _sfms2.resize((int)_coords.size());
    startTimer();
    setScoreThresholds();
    std::cerr << "done\t" << curTimer()/CLOCKS_PER_SEC << " secs" << std::endl;
    startTimer();
    std::cerr << "slimNeighborsSetLinks..." << flush;
    slimNeighborsSetLinks();
    std::cerr << "done\t" << curTimer()/CLOCKS_PER_SEC << " secs" << std::endl;
  }

  // Improve visibility by using texture analysis
  startTimer();
  std::cerr << "mergeSFM..." << flush;
  mergeSfMP();
  std::cerr << '\t' << curTimer()/CLOCKS_PER_SEC << " secs" << std::endl;

  _sfms2.clear();
  _sfms2.resize((int)_coords.size());
  std::cerr << "setScoreThresholds..." << flush;
  startTimer();
  setScoreThresholds();
  std::cerr << "done\t" << curTimer()/CLOCKS_PER_SEC << " secs" << std::endl;

  // Remove redundant images first
  std::cerr << "sRemoveImages... " << flush;
  startTimer();

  sRemoveImages();
  std::cerr << '\t' << curTimer()/CLOCKS_PER_SEC << " secs" << std::endl;

  // use _removed to change _visibles and update _neighbors
  startTimer();
  resetVisibles();
  setNeighbors();
  std::cerr << "slimNeighborsSetLinks..." << flush;
  slimNeighborsSetLinks();
  std::cerr << "done\t" << curTimer()/CLOCKS_PER_SEC << " secs" << std::endl;

  // Init _timages by mutually exclusive clustering
  setTimages();
  _oimages.resize((int)_timages.size());
}

void CMVS::CBundle::run(const std::string prefix, const int imageThreshold, const int tau, const float scoreRatioThreshold, const float coverageThreshold, const int pnumThreshold, const int CPU) {
  startTimer();

  prep(prefix, imageThreshold, tau, scoreRatioThreshold, coverageThreshold, pnumThreshold, CPU);

  std::cerr << '\t' << curTimer()/CLOCKS_PER_SEC << " secs" << std::endl;

  prep2();

  // Assumed variables that must be set properly here
  std::cerr << "Adding images: " << std::endl;
  startTimer();
  // Add images
  // Repeat until all the clusters become at most _imageThreshold.
  while (1) {
    addImagesP();

    int change = 0;
    std::vector<std::vector<int> > newtimages;
    cout << "Divide: " << flush;
    for (int i = 0; i < (int)_timages.size(); ++i) {
      if ((int)_timages[i].size() <= _imageThreshold) {
        newtimages.push_back(_timages[i]);
        continue;
      } else {
        cout << i << ' ';

        change = 1;
        // divide
        std::vector<std::vector<int> > vvi;
        divideImages(_timages[i], vvi);
        for (int j = 0; j < (int)vvi.size(); ++j) {
          newtimages.push_back(vvi[j]);
        }
      }
    }

    cout << std::endl;

    _timages.swap(newtimages);
    if (change == 0) break;
  }
  std::cerr << "done\t" << curTimer()/CLOCKS_PER_SEC << " secs" << std::endl;

  _oimages.resize((int)_timages.size());

  writeCameraCenters();

  // Output results
  writeVis();
  writeGroups();
}

float CMVS::CBundle::computeLink(const int image0, const int image1) {
  std::vector<int> common;
  set_intersection(_vpoints[image0].begin(), _vpoints[image0].end(), _vpoints[image1].begin(), _vpoints[image1].end(), back_inserter(common));

  float score = 0.0f;
  for (int i = 0; i < (int)common.size(); ++i) {
    const int pid = common[i];
    std::vector<int> vtmp;
    vtmp.push_back(image0);
    vtmp.push_back(image1);
    const float ftmp = computeScore2(_coords[pid], vtmp);
    if (_sfms2[pid]._score != 0.0f) {
      score += _pweights[pid] * ftmp / (_sfms2[pid]._scoreThreshold / _scoreRatioThreshold);
    }
  }

  return score;
}

void CMVS::CBundle::slimNeighborsSetLinks() {
  const int maxneighbor = 30;
  _links.clear();
  _links.resize(_cnum);

#pragma omp parallel for
  for (int image = 0; image < _cnum; ++image) {
    _links[image].resize((int)_neighbors[image].size());

    for (int i = 0; i < (int)_neighbors[image].size(); ++i) {
      _links[image][i] = computeLink(image, _neighbors[image][i]);
    }

    if ((int)_neighbors[image].size() < 2) continue;

    std::vector<int> newneighbors;
    std::vector<float> newlinks;

    std::vector<Vec2f> vv;
    for (int i = 0; i < (int)_neighbors[image].size(); ++i) {
      vv.push_back(Vec2(-_links[image][i], _neighbors[image][i]));
    }
    std::sort(vv.begin(), vv.end(), Svec2cmp<float>());

    const int itmp = std::min(maxneighbor, (int)vv.size());
    for (int i = 0; i < itmp; ++i) {
      newneighbors.push_back((int)vv[i][1]);
      newlinks.push_back(-vv[i][0]);
    }

    _neighbors[image].swap(newneighbors);
    _links[image].swap(newlinks);
  }
}

void CMVS::CBundle::setScoreThresholds() {
#pragma omp parallel for
  for (int p = 0; p < (int)_coords.size(); ++p) {
    _sfms2[p]._scoreThreshold = computeScore2(_coords[p], _visibles[p], _sfms2[p]._uimages) * _scoreRatioThreshold;
  }
}

void CMVS::CBundle::sRemoveImages() {
  _removed.clear();
  _removed.resize(_cnum, 0);

  _allows.resize(_cnum);
  for (int c = 0; c < _cnum; ++c) {
    _allows[c] = (int)ceil((int)_vpoints[c].size() * (1.0f - _coverageThreshold));
  }

  // Sort all the images in an increasing order of resolution
  std::vector<Vec2i> vvi;
  for (int c = 0; c < _cnum; ++c) {
    const int res = _pss._photos[c].getWidth(0) * _pss._photos[c].getHeight(0);
    vvi.push_back(Vec2i(res, c));
  }
  std::sort(vvi.begin(), vvi.end(), Svec2cmp<int>());

  const int tenth = std::max(1, _cnum / 10);
  for (int i = 0; i < (int)vvi.size(); ++i) {
    if (i % tenth == 0) {
      std::cerr << '*' << flush;
    }
    const int image = vvi[i][1];
    checkImage(image);
  }
  std::cerr << std::endl;

  std::cerr << "Kept: ";
  int kept = 0;
  for (int c = 0; c < _cnum; ++c) {
    if (_removed[c] == 0) {
      ++kept;
      std::cerr << c << ' ';
    }
  }
  std::cerr << std::endl << std::endl;

  std::cerr << "Removed: ";
  for (int c = 0; c < _cnum; ++c) {
    if (_removed[c]) {
      std::cerr << c << ' ';
    }
  }
  std::cerr << std::endl;
  std::cerr << "sRemoveImages: " << _cnum << " -> " << kept << flush;
}

void CMVS::CBundle::resetVisibles(void) {
  // reset _visibles. remove "removed images" from the list.
  for (int p = 0; p < (int)_visibles.size(); ++p) {
    std::vector<int> newimages;

    setNewImages(p, -1, newimages);
    _visibles[p].swap(newimages);
  }
}

void CMVS::CBundle::setNewImages(const int pid, const int rimage, std::vector<int>& newimages) {
  for (int i = 0; i < (int)_visibles[pid].size(); ++i) {
    const int image = _visibles[pid][i];
    if (_removed[image] || image == rimage) continue;

    newimages.push_back(image);
  }
}

void CMVS::CBundle::checkImage(const int image) {
  // For each SfM, check if it will be unsatisfied by removing image.
  //  0: unsatisfy
  //  1: satisfy->satisfy
  //  2: satisfy->unsatisfy
  _statsT.resize((int)_vpoints[image].size());
  _jobs.clear();

#pragma omp parallel for
  for (int p = 0; p < (int)_statsT.size(); ++p) {
    const int pid = _vpoints[image][p];
    if (_sfms2[pid]._satisfied == 0) {
      _statsT[p] = 0;
    } else {
      _statsT[p] = 1;

      // if _uimages of p is not removed, if image is not in
      // _uimages, the point is still satisfied.
      int valid = 1;
      int inside = 0;

      for (int i = 0; i < (int)_sfms2[pid]._uimages.size(); ++i) {
        const int itmp = _sfms2[pid]._uimages[i];
        if (itmp == image) {
          inside = 1;
        }

        if (_removed[itmp]) {
          valid = 0;
          break;
        }
      }

      if (valid == 1 && inside == 0) continue;

      std::vector<int> newimages;
      setNewImages(pid, image, newimages);
      const float cscore = computeScore2(_coords[pid], newimages);
      if (cscore < _sfms2[pid]._scoreThreshold) {
        _statsT[p] = 2;
      }
    }
  }

  // For each image, how many SFM points are removed if you remove
  // "image"
  std::vector<int> decrements;
  decrements.resize(_cnum, 0);

  // Look at points with _stastT[p] = 2, to see if _allows are still
  // above thresholds.
  for (int p = 0; p < (int)_statsT.size(); ++p) {
    if (_statsT[p] == 2) {
      const int pid = _vpoints[image][p];
      for (int i = 0; i < (int)_visibles[pid].size(); ++i) {
        const int itmp = _visibles[pid][i];
        ++decrements[itmp];
      }
    }
  }

  // If _allows can cover decrements, go ahead.
  int rflag = 1;
  for (int c = 0; c < _cnum; ++c) {
    if (_allows[c] < decrements[c]) {
      rflag = 0;
      break;
    }
  }

  // remove
  if (rflag) {
    _removed[image] = 1;
    for (int p = 0; p < (int)_statsT.size(); ++p) {
      if (_statsT[p] == 2) {
        _sfms2[_vpoints[image][p]]._satisfied = 0;
      }
    }

    for (int c = 0; c < _cnum; ++c) {
      _allows[c] -= decrements[c];
    }

    // Update uimages for points that are still satisfied and still
    // contains image in _uimages
#pragma omp parallel for
    for (int p = 0; p < (int)_statsT.size(); ++p) {
      const int pid = _vpoints[image][p];
      if (_statsT[p] == 1) {
        int contain = 0;
        for (int i = 0; i < (int)_sfms2[pid]._uimages.size(); ++i) {
          if (_sfms2[pid]._uimages[i] == image) {
            contain = 1;
            break;
          }
        }

        if (contain) {
          std::vector<int> newimages;
          setNewImages(pid, -1, newimages);
          const float cscore = computeScore2(_coords[pid], newimages, _sfms2[pid]._uimages);
          if (cscore < _sfms2[pid]._scoreThreshold) {
            _sfms2[pid]._satisfied = 0;
          }
        }
      }
    }
  }
}

void CMVS::CBundle::setNeighbors(void) {
  _neighbors.clear();
  _neighbors.resize(_cnum);
#pragma omp parallel for
  for (int image = 0; image < _cnum; ++image) {
    std::vector<int> narray;
    narray.resize(_cnum, 0);

    for (int p = 0; p < (int)_visibles.size(); ++p) {
      if (!binary_search(_visibles[p].begin(), _visibles[p].end(), image)) continue;

      for (int i = 0; i < (int)_visibles[p].size(); ++i) {
        narray[_visibles[p][i]] = 1;
      }
    }

    for (int i = 0; i < _cnum; ++i) {
      if (narray[i] && i != image) {
        _neighbors[image].push_back(i);
      }
    }
  }
}

void CMVS::CBundle::setTimages(void) {
  std::vector<int> lhs;
  for (int c = 0; c < _cnum; ++c) {
    if (_removed[c] == 0) {
      lhs.push_back(c);
    }
  }

  _timages.clear();

  if ((int)lhs.size() <= _imageThreshold) {
    _timages.push_back(lhs);
  } else {
    divideImages(lhs, _timages);
  }

  std::cerr << std::endl << "Cluster sizes: " << std::endl;
  for (int i = 0; i < (int)_timages.size(); ++i) {
    std::cerr << (int)_timages[i].size() << ' ';
  }
  std::cerr << std::endl;
}

void CMVS::CBundle::divideImages(const std::vector<int>& lhs, std::vector<std::vector<int> >& rhs) {
  const float iratio = 125.0f / 150.0f;
  // Candidates
  list<std::vector<int> > candidates;
  // initialize first cluster
  candidates.push_back(std::vector<int>());
  std::vector<int>& vitmp = candidates.back();
  for (int i = 0; i < (int)lhs.size(); ++i) {
    vitmp.push_back(lhs[i]);
  }

  // Iterate
  while (1) {
    if (candidates.empty()) break;

    std::vector<int> cand = candidates.front();
    candidates.pop_front();
    if ((int)cand.size() <= _imageThreshold * iratio) {
      rhs.push_back(cand);
      continue;
    }

    // Divide into 2 groups
    std::vector<idxtype> xadj, adjncy, adjwgt, part;
    const int nparts = 2;
    const int cutType = 0;
    // Define graphs
    for (int i = 0; i < (int)cand.size(); ++i) {
      xadj.push_back((int)adjncy.size());

      for (int j = 0; j < (int)cand.size(); ++j) {
        if (i == j) continue;

        // Check if cand[i] and cand[j] are connected
        const int cid0 = cand[i];
        const int cid1 = cand[j];
        std::vector<int>::const_iterator nite = find(_neighbors[cid0].begin(), _neighbors[cid0].end(), cid1);

        if (nite != _neighbors[cid0].end()) {
          adjncy.push_back(j);
          const int offset = nite - _neighbors[cid0].begin();
          const int weight = std::min(5000, (int)floor(10.0f * _links[cid0][offset] + 0.5f));
          adjwgt.push_back(weight);
        }
      }
    }
    xadj.push_back((int)adjncy.size());

    CGraclus::runE(xadj, adjncy, adjwgt, nparts, cutType, part);

    // Divide into 2 groups
    std::vector<int> cand1, cand2;
    for (int i = 0; i < (int)part.size(); ++i) {
      if (part[i] == 0) {
        cand1.push_back(cand[i]);
      } else {
        cand2.push_back(cand[i]);
      }
    }

    if (cand1.empty() || cand2.empty()) {
      std::cerr << "Error. Normalized cuts produced an empty cluster: "
           << (int)part.size() << " -> "
           << (int)cand1.size() << ' '
           << (int)cand2.size() << std::endl;
      exit (1);
    }

    if ((int)cand1.size() <= _imageThreshold * iratio) {
      rhs.push_back(cand1);
    } else {
      candidates.push_back(std::vector<int>());
      (candidates.back()).swap(cand1);
    }

    if ((int)cand2.size() <= _imageThreshold * iratio) {
      rhs.push_back(cand2);
    } else {
      candidates.push_back(std::vector<int>());
      (candidates.back()).swap(cand2);
    }
  }
}

void CMVS::CBundle::readBundle(const std::string file) {
  // For each valid image, a list of connected images, and the second
  // value is the number of common points.
  std::ifstream ifstr;  ifstr.open(file.c_str());
  if (!ifstr.is_open()) {
    std::cerr << "Bundle file not found: " << file << std::endl;
    exit (1);
  }

  while (1) {
    unsigned char uctmp;
    ifstr.read((char*)&uctmp, sizeof(unsigned char));
    ifstr.putback(uctmp);
    if (uctmp == '#') {
      char buffer[1024];      ifstr.getline(buffer, 1024);
    } else {
      break;
    }
  }
  int cnum, pnum;
  ifstr >> cnum >> pnum;
  std::vector<int> ids;
  ids.resize(cnum);

  std::cerr << cnum << " cameras -- " << pnum << " points in bundle file" << std::endl;
  _cnum = 0;
  for (int c = 0; c < cnum; ++c) {
    ids[c] = -1;
    float params[15];
    for (int i = 0; i < 15; ++i) {
      ifstr >> params[i];
    }
    if (params[0] != 0.0f) {
      ids[c] = _cnum++;
    }
  }

  _vpoints.resize(_cnum);

  _coords.clear();
  _visibles.clear();
  _colors.clear();
  _pnum = 0;

  const int tenth = std::max(1, pnum / 10);
  for (int p = 0; p < pnum; ++p) {
    if (p % tenth == 0) {
      std::cerr << '*' << flush;
    }
    int num;
    Vec3f color;
    Vec4f coord;
    ifstr >> coord[0] >> coord[1] >> coord[2] >> color >> num;
    coord[3] = 1.0f;

    std::vector<int> visibles;
    for (int i = 0; i < num; ++i) {
      int itmp;
      ifstr >> itmp;
      if (cnum <= itmp) {
        double dtmp;
        ifstr >> itmp >> dtmp >> dtmp;
        continue;
      }

      if (ids[itmp] == -1) {
        std::cerr << "impossible " << itmp << ' ' << ids[itmp] << std::endl;
        exit (1);
      }
      visibles.push_back(ids[itmp]);

      // Based on the bundler version, the number of parameters here
      // are either 1 or 3. Currently, it seems to be 3.
      double dtmp;
      ifstr >> itmp >> dtmp >> dtmp;
      //ifstr >> dtmp;
    }

    if ((int)visibles.size() < 2) continue;

    std::sort(visibles.begin(), visibles.end());

    for (int i = 0; i < (int)visibles.size(); ++i) {
       _vpoints[visibles[i]].push_back(_pnum);
    }

    _visibles.push_back(visibles);
    _coords.push_back(coord);
    _colors.push_back(color);
    ++_pnum;
  }
  ifstr.close();
  setNeighbors();

  std::cerr << std::endl << _cnum << " cameras -- " << _pnum << " points" << flush;
}

void CMVS::CBundle::findPNeighbors(sfcnn<const float*, 3, float>& tree, const int pid, std::vector<int>& pneighbors) {
  const float thresh = _dscale2 * _dscale2 * _minScales[pid] * _minScales[pid];

  std::vector<long unsigned int> ids;
  std::vector<double> dists;
  int k = std::min((int)_coords.size(), 400);

  while (1) {
    ids.clear();
    dists.clear();
    tree.ksearch(&_coords[pid][0], k, ids, dists);

    if (thresh < dists[k - 1]) break;
    if (k == (int)_coords.size()) break;

    k = std::min((int)_coords.size(), k * 2);
  }

  for (int i = 0; i < k; ++i) {
    if (thresh < dists[i]) break;

    // If distance from pid2 is also within threshold ok
    const float thresh2 = _dscale2 * _dscale2 * _minScales[ids[i]] * _minScales[ids[i]];
    if (thresh2 < dists[i]) continue;

    if (ids[i] != (unsigned int)pid) {
      pneighbors.push_back(ids[i]);
    }
  }
}

void CMVS::CBundle::mergeSfMPThread(void) {
  const int tenth = std::max(1, (int)_coords.size() / 10);

  while (1) {
    int pid = -1;
    _lock.lock();
    if (!_jobs.empty()) {
      pid = _jobs.front();
      _jobs.pop_front();
    }
    if (pid != -1 && _merged[pid]) {
      pid = -2;
    }

    if (_count % tenth == 0) {
      std::cerr << '*' << flush;
    }
    ++_count;
    _lock.unlock();

    if (pid == -2) continue;
    if (pid == -1) break;

    // Process, and check later if merged flag is on
    std::vector<int> pneighbors;
    findPNeighbors(*_ptree, pid, pneighbors);
    const int psize = (int)pneighbors.size();
    // visible images and its neighbors
    std::vector<int> visibles = _visibles[pid];
    std::vector<int> vitmp;
    for (int i = 0; i < (int)_visibles[pid].size(); ++i) {
      const int imagetmp = _visibles[pid][i];
      vitmp.clear();
      mymerge(visibles, _neighbors[imagetmp], vitmp);
      vitmp.swap(visibles);
    }

    std::vector<char> visflag(psize, 0);
    // test visibility
    for (int i = 0; i < psize; ++i) {
      const int& pid2 = pneighbors[i];
      if (my_isIntersect(visibles, _visibles[pid2])) {
        visflag[i] = 1;
      }
    }

    // Now lock and try to register
    _lock.lock();
    // If the main one is removed, over... waste.
    if (_merged[pid] == 0) {
      _merged[pid] = 1;
      for (int i = 0; i < psize; ++i) {
        const int pid2 = pneighbors[i];
        if (visflag[i] && _merged[pid2] == 0) {
          _merged[pid2] = 1;
          // _puf->union_set(pid, pid2);
          _puf2->union_sets(pid, pid2);
        }
      }
    }
    _lock.unlock();
  }
}

// Arbitrary seed for deterministic pseudorandomness.
static const unsigned int RANDOM_SEED = 42;

void CMVS::CBundle::mergeSfMP(void) {
  // Repeat certain times until no change
  const int cpnum = (int)_coords.size();
  _minScales.resize(cpnum);
  for (int p = 0; p < cpnum; ++p) {
    _minScales[p] = INT_MAX/2;
    for (int i = 0; i < (int)_visibles[p].size(); ++i) {
      const int image = _visibles[p][i];
      const float stmp = _pss._photos[image].getScale(_coords[p], _levels[image]);
      _minScales[p] = std::min(_minScales[p], stmp);
    }
  }

  // _puf = new boost::disjoint_sets_with_storage<>(cpnum);
  _puf2 = new CDisjointSet<int>();
  for (int p = 0; p < cpnum; ++p) {
    // _puf->make_set(p);
    _puf2->add_element(p);
  }

  {
    // kdtree
    std::vector<const float*> ppoints;
    ppoints.resize((int)_coords.size());
    for (int i = 0; i < (int)_coords.size(); ++i) {
      ppoints[i] = &(_coords[i][0]);
    }
    _ptree = new sfcnn<const float*, 3, float>(&ppoints[0], (int)ppoints.size());

    _merged.resize((int)_coords.size(), 0);
    _jobs.clear();
    std::vector<int> order;
    order.resize(cpnum);
    for (int p = 0; p < cpnum; ++p) {
      order[p] = p;
    }
    std::mt19937 gen(RANDOM_SEED);
    shuffle(order.begin(), order.end(), gen);

    for (int p = 0; p < cpnum; ++p) {
      _jobs.push_back(order[p]);
    }

    _count = 0;
    std::vector<std::thread> threads(_CPU);
    for (int c = 0; c < _CPU; ++c) {
      threads[c] = std::thread(&CBundle::mergeSfMPThread, this);
    }

    for (int c = 0; c < _CPU; ++c) {
      if (threads[c].joinable()) {
        threads[c].join();
      }
    }

    int newpnum = 0;
    // Mapping from _pnum to new id for reps
    std::vector<int> dict;
    std::vector<int> reps;
    dict.resize((int)_coords.size(), -1);
    for (int p = 0; p < (int)_coords.size(); ++p) {
      // const int ptmp = _puf->find_set(p);
      const int ptmp = _puf2->find_set(p);
      if (p == ptmp) {
        dict[p] = newpnum;
        reps.push_back(p);
        ++newpnum;
      }
    }
  }

  // Based on _puf, reset _coords, _coords, _visibles, _vpoints
  std::cerr << "resetPoints..." << flush;
  resetPoints();
  std::cerr << "done" << std::endl;

  // delete _puf;
  delete _puf2;
  // _puf = NULL;
  _puf2 = nullptr;
  delete _ptree;
  _ptree = NULL;

  const int npnum = (int)_coords.size();
  std::cerr << "Rep counts: " << cpnum << " -> " << npnum << "  " << flush;
}

// Based on _puf, reset _coords, _coords, _visibles, _vpoints
void CMVS::CBundle::resetPoints(void) {
  std::vector<int> counts;
  std::vector<int> smallestids;
  counts.resize((int)_coords.size(), 0);
  smallestids.resize((int)_coords.size(), INT_MAX/2);
  for (int p = 0; p < (int)_coords.size(); ++p) {
    // const int ptmp = _puf->find_set(p);
    const int ptmp = _puf2->find_set(p);
    smallestids[ptmp] = std::min(smallestids[ptmp], p);
    ++counts[ptmp];
  }
  const int mthreshold = 2;
  std::vector<Vec2i> vv2;
  for (int p = 0; p < (int)_coords.size(); ++p) {
    if (mthreshold <= counts[p]) {
      vv2.push_back(Vec2i(smallestids[p], p));
    }
  }
  std::sort(vv2.begin(), vv2.end(), Svec2cmp<int>());

  std::vector<int> newpids;
  newpids.resize((int)_coords.size(), -1);
  int newpnum = (int)vv2.size();
  for (int i = 0; i < newpnum; ++i) {
    newpids[vv2[i][1]] = i;
  }

  std::vector<Vec4f> newcoords;  newcoords.resize(newpnum);
  std::vector<std::vector<int> > newvisibles;  newvisibles.resize(newpnum);

  for (int p = 0; p < (int)_coords.size(); ++p) {
    // const int ptmp = _puf->find_set(p);
    const int ptmp = _puf2->find_set(p);
    if (counts[ptmp] < mthreshold) continue;

    const int newpid = newpids[ptmp];
    newcoords[newpid] += _coords[p];
    std::vector<int> vitmp;
    mymerge(newvisibles[newpid], _visibles[p], vitmp);
    vitmp.swap(newvisibles[newpid]);
  }

  for (int i = 0; i < newpnum; ++i) {
    newcoords[i] = newcoords[i] / newcoords[i][3];
  }

  // Update _vpoints
  _coords.swap(newcoords);
  _visibles.swap(newvisibles);

  for (int c = 0; c < _cnum; ++c) {
    _vpoints[c].clear();
  }

  for (int p = 0; p < (int)_coords.size(); ++p) {
    for (int i = 0; i < (int)_visibles[p].size(); ++i) {
      const int itmp = _visibles[p][i];
      _vpoints[itmp].push_back(p);
    }
  }

  _pweights.clear();
  for (int i = 0; i < newpnum; ++i) {
    _pweights.push_back(counts[vv2[i][1]]);
  }
}

void CMVS::CBundle::setScoresClusters(void) {
#pragma omp parallel for
  for (int p = 0; p < (int)_coords.size(); ++p) {
    // if _satisfied is 0, no hope (even all the images are in the
    // cluster, not satisfied
    if (_sfms2[p]._satisfied == 0) continue;

    _sfms2[p]._satisfied = 2;
    setCluster(p);
  }
}

void CMVS::CBundle::setCluster(const int p) {
  _sfms2[p]._score = -1.0f;
  _sfms2[p]._cluster = -1;
  for (int c = 0; c < (int)_timages.size(); ++c) {
    // Select images in cluster
    std::vector<int> vitmp;
    set_intersection(_timages[c].begin(), _timages[c].end(), _visibles[p].begin(), _visibles[p].end(), back_inserter(vitmp));

    const float stmp = computeScore2(_coords[p], vitmp);
    if (_sfms2[p]._score < stmp) {
      _sfms2[p]._cluster = c;
      _sfms2[p]._score = stmp;
    }
  }

  // If no cluster contains 2 images
  if (_sfms2[p]._cluster == -1) {
    int find = 0;
    for (int j = 0; j < (int)_visibles[p].size(); ++j) {
      if (find) break;

      for (int c = 0; c < (int)_timages.size(); ++c)
        if (binary_search(_timages[c].begin(), _timages[c].end(), _visibles[p][j])) {
          _sfms2[p]._cluster = c;
          _sfms2[p]._score = 0.0f;
          find = 1;
          break;
        }
    }

    // If none of the visibles images are
    if (find == 0) {
      std::cerr << "Impossible in setscoresclustersthread" << std::endl
           << (int)_visibles[p].size() << std::endl
           << (int)_sfms2[p]._satisfied << std::endl;
      exit (1);
    }
  }

  if (_sfms2[p]._scoreThreshold <= _sfms2[p]._score) {
    // SATISFIED
    _sfms2[p]._satisfied = 1;

    //  update _lacks
    _lock.lock();
    for (int i = 0; i < (int)_visibles[p].size(); ++i) {
      const int image = _visibles[p][i];
      --_lacks[image];
    }
    _lock.unlock();
  }
}

float CMVS::CBundle::angleScore(const Vec4f& ray0, const Vec4f& ray1) {
  const static float lsigma = 5.0f * M_PI / 180.f;
  const static float rsigma = 15.0f * M_PI / 180.0f;
  const static float lsigma2 = 2.0f * lsigma * lsigma;
  const static float rsigma2 = 2.0f * rsigma * rsigma;
  const static float pivot = 20.0f * M_PI / 180.0f;

  const float angle = acos(std::min(1.0f, ray0 * ray1));
  const float diff = angle - pivot;

  if (angle < pivot) {
    return exp(- diff * diff / lsigma2);
  } else {
    return exp(- diff * diff / rsigma2);
  }
}

void CMVS::CBundle::setClusters(void) {
#pragma omp parallel for
  for (int p = 0; p < (int)_sfms2.size(); ++p) {
    if (_sfms2[p]._satisfied != 2) continue;

    setCluster(p);
  }
}

void CMVS::CBundle::addImagesP(void) {
  // set _lacks (how many more sfm points need to be satisfied)
  _lacks.resize(_cnum);
  for (int c = 0; c < _cnum; ++c) {
    if (_removed[c]) {
      _lacks[c] = 0;
    } else {
      _lacks[c] = (int)floor((int)_vpoints[c].size() * _coverageThreshold);
    }
  }

  // Compute best score given all the images. Upper bound on the score
  setScoresClusters();

  // Add images to clusters to make _lacks at most 0 for all the
  // images. In practice, for each cluster, identify corresponding sfm
  // points that have not been satisfied. These points compute gains.

  // set _adds for each sfm point, and add images

  _addnums.clear();
  _addnums.resize((int)_timages.size(), 0);

  while (1) {
    const int totalnum = addImages();
    // Very end
    if (totalnum == 0) break;

    // If one of the cluster gets more than _imageThreshold, divide
    int divide = 0;
    for (int i = 0; i < (int)_timages.size(); ++i) {
      if (_imageThreshold < (int)_timages[i].size()) {
        divide = 1;
        break;
      }
    }
    if (divide) break;

    setClusters();
  }

  for (int i = 0; i < (int)_addnums.size(); ++i) {
    std::cerr << _addnums[i] << ' ';
  }
  std::cerr << std::endl;

  int totalnum = 0;
  for (int c = 0; c < (int)_timages.size(); ++c) {
    totalnum += (int)_timages[c].size();
  }
  int beforenum = 0;
  for (int c = 0; c < _cnum; ++c) {
    if (_removed[c] == 0) {
      ++beforenum;
    }
  }

  std::cerr << "Image nums: "
       << _cnum << " -> " << beforenum <<  " -> " << totalnum << std::endl;
}

int CMVS::CBundle::addImages(void) {
  _thread = 0;
  _jobs.clear();

  // we think about sfm points that belong to images that have not been satisfied
  std::vector<char> flags;
  flags.resize((int)_coords.size(), 0);
  for (int c = 0; c < _cnum; ++c) {
    if (_lacks[c] <= 0) continue;

    for (int i = 0; i < (int)_vpoints[c].size(); ++i) {
      const int pid = _vpoints[c][i];
      if (_sfms2[pid]._satisfied == 2) {
        flags[pid] = 1;
      }
    }
  }

#pragma omp parallel for
  for (int p = 0; p < (int)_sfms2.size(); ++p) {
    if (flags[p] == 0) continue;

    _sfms2[p]._adds.clear();

    const int cluster = _sfms2[p]._cluster;
    // Try to add an image to _timages[cluster]
    std::vector<int> cimages;
    set_intersection(_timages[cluster].begin(), _timages[cluster].end(), _visibles[p].begin(), _visibles[p].end(), back_inserter(cimages));
    std::sort(cimages.begin(), cimages.end());

    for (int i = 0; i < (int)_visibles[p].size(); ++i) {
      const int image = _visibles[p][i];

      if (binary_search(cimages.begin(), cimages.end(), image)) continue;

      std::vector<int> vitmp = cimages;
      vitmp.push_back(image);
      const float newscore = computeScore2(_coords[p], vitmp);
      if (newscore <= _sfms2[p]._score) continue;

      _sfms2[p]._adds.push_back(Sadd(image, (newscore - _sfms2[p]._score) / _sfms2[p]._scoreThreshold));
    }
  }

  // Accumulate information from SfM points. For each cluster.
  // For each cluster, for each image, sum of gains.
  std::vector<std::map<int, float> > cands;
  cands.resize((int)_timages.size());

  for (int p = 0; p < (int)_sfms2.size(); ++p) {
    if (flags[p] == 0) continue;

    const int cluster = _sfms2[p]._cluster;

    for (int i = 0; i < (int)_sfms2[p]._adds.size(); ++i) {
      Sadd& stmp = _sfms2[p]._adds[i];
      cands[cluster][stmp._image] += stmp._gain;
    }
  }

  return addImagesSub(cands);
}

int CMVS::CBundle::addImagesSub(const std::vector<std::map<int, float>>& cands) {
  // Vec3f (gain, cluster, image). Start adding starting from the good
  // one, block neighboring images.
  std::vector<Vec3f> cands2;
  for (int i = 0; i < (int)_timages.size(); ++i) {
    std::map<int, float>::const_iterator mbegin = cands[i].begin();
    std::map<int, float>::const_iterator mend = cands[i].end();

    while (mbegin != mend) {
      cands2.push_back(Vec3f(-mbegin->second, i, mbegin->first));
      ++mbegin;
    }
  }

  if (cands2.empty()) return 0;

  std::sort(cands2.begin(), cands2.end(), Svec3cmp<float>());

  std::vector<char> blocked;
  blocked.resize(_cnum, 0);
  std::vector<int> addnum;
  addnum.resize((int)_timages.size(), 0);

  // A bit of tuning is possible here. For the paper, we used 0.7f,
  //but 0.9f seems to produce better results in general.  const float
  //gainThreshold = -cands2[0][0] * 0.7f;
  const float gainThreshold = -cands2[0][0] * 0.9f;

  // use rato threshold for blocked one. if not blocked, just keep on
  // going.
  for (int i = 0; i < (int)cands2.size(); ++i) {
    const float gain = -cands2[i][0];
    if (gain < gainThreshold) break;

    const int cluster = (int)cands2[i][1];
    const int image = (int)cands2[i][2];

    if (blocked[image]) continue;

    // Add
    ++addnum[cluster];
    ++_addnums[cluster];
    blocked[image] = 1;

    for (int i = 0; i < (int)_neighbors[image].size(); ++i) {
      blocked[_neighbors[image][i]] = 1;
    }

    _timages[cluster].push_back(image);
    // _score, _cluster, _satisfied, _lacks are updated in
    // setClusters.  So, no need to update anything.
  }

  for (int i = 0; i < (int)_timages.size(); ++i) {
    std::sort(_timages[i].begin(), _timages[i].end());
  }

  return std::accumulate(addnum.begin(), addnum.end(), 0);
}

int CMVS::CBundle::totalNum(void) const {
  int totalnum = 0;
  for (int c = 0; c < (int)_timages.size(); ++c) {
    totalnum += (int)_timages[c].size();
  }
  return totalnum;
}

int CMVS::CBundle::my_isIntersect(const std::vector<int>& lhs, const std::vector<int>& rhs) {
  std::vector<int>::const_iterator b0 = lhs.begin();
  std::vector<int>::const_iterator e0 = lhs.end();

  std::vector<int>::const_iterator b1 = rhs.begin();
  std::vector<int>::const_iterator e1 = rhs.end();

  while (1) {
    if (b0 == e0) return 0;

    if (b1 == e1) return 0;

    if (*b0 == *b1) {
      return 1;
    } else if (*b0 < *b1) {
      ++b0;
    } else {
      ++b1;
    }
  }
}

void CMVS::CBundle::mymerge(const std::vector<int>& lhs, const std::vector<int>& rhs, std::vector<int>& output) {
  std::vector<int>::const_iterator b0 = lhs.begin();
  std::vector<int>::const_iterator e0 = lhs.end();

  std::vector<int>::const_iterator b1 = rhs.begin();
  std::vector<int>::const_iterator e1 = rhs.end();

  while (1) {
    if (b0 == e0) {
      output.insert(output.end(), b1, e1);
      break;
    }
    if (b1 == e1) {
      output.insert(output.end(), b0, e0);
      break;
    }

    if (*b0 == *b1) {
      output.push_back(*b0);
      ++b0;
      ++b1;
    } else if (*b0 < *b1) {
      output.push_back(*b0);
      ++b0;
    } else {
      output.push_back(*b1);
      ++b1;
    }
  }
}

// Calculates log2 of number.
template <typename T>
T Log2(T n) {
  return log(n) / log(T(2));
}

void CMVS::CBundle::setWidthsHeightsLevels(void) {
  _widths.resize(_cnum);
  _heights.resize(_cnum);
  _levels.resize(_cnum);

  // Determine level. SfM was done on 2M pixels.
  for (int c = 0; c < _cnum; ++c) {
    _levels[c] = _dlevel;

    _widths[c] = _pss._photos[c].getWidth(_levels[c]);
    _heights[c] = _pss._photos[c].getHeight(_levels[c]);
  }
}


float CMVS::CBundle::computeScore2(const Vec4f& coord, const std::vector<int>& images) const {
  std::vector<int> uimages;
  return computeScore2(coord, images, uimages);
}

float CMVS::CBundle::computeScore2(const Vec4f& coord, const std::vector<int>& images, std::vector<int>& uimages) const {
  const int inum = (int)images.size();
  uimages.clear();
  if (inum < 2) return -1.0f;

  std::vector<Vec4f> rays;
  rays.resize(inum);

  std::vector<float> scales;
  scales.resize(inum);

  for (int r = 0; r < inum; ++r) {
    rays[r] = _pss._photos[images[r]].OpticalCenter() - coord;
    unitize(rays[r]);

    scales[r] = 1.0f / _pss._photos[images[r]].getScale(coord, 0);
  }

  // Find the best pair of images
  Vec2i bestpair;
  float bestscore = -INT_MAX/2;
  for (int i = 0; i < inum; ++i) {
    for (int j = i+1; j < inum; ++j) {
      const float ftmp = angleScore(rays[i], rays[j]) * scales[i] * scales[j];

      if (bestscore < ftmp) {
        bestscore = ftmp;
        bestpair = Vec2i(i, j);
      }
    }
  }

  std::vector<int> inout;
  inout.resize(inum, 1);
  inout[bestpair[0]] = 0;
  inout[bestpair[1]] = 0;
  uimages.push_back(bestpair[0]);
  uimages.push_back(bestpair[1]);

  for (int t = 2; t < std::min(_tau, inum); ++t) {
    // Add the best new image
    int ansid = -1;
    float ansscore = -INT_MAX/2;
    for (int j = 0; j < inum; ++j) {
      if (inout[j] == 0) continue;

      float score = 0.0f;
      for (int k = 0; k < (int)uimages.size(); ++k) {
        const int iid = uimages[k];
        score += angleScore(rays[j], rays[iid]) * scales[j] * scales[iid];
      }

      if (ansscore < score) {
        ansscore = score;
        ansid = j;
      }
    }

    if (ansid == -1) {
      std::cerr << "Impossible 2 in compureScore" << std::endl;      exit (1);
    }

    inout[ansid] = 0;
    uimages.push_back(ansid);
    bestscore += ansscore;
  }

  for (int i = 0; i < (int)uimages.size(); ++i) {
    uimages[i] = images[uimages[i]];
  }

  return bestscore;
}

float CMVS::CBundle::computeScore2(const Vec4f& coord, const std::vector<int>& images, const int index) const {
  std::vector<int> uimages;
  const int inum = (int)images.size();
  if (inum < 2) return -1.0f;

  std::vector<Vec4f> rays;    rays.resize(inum);
  std::vector<float> scales;  scales.resize(inum);

  for (int r = 0; r < inum; ++r) {
    rays[r] = _pss._photos[images[r]].OpticalCenter() - coord;
    unitize(rays[r]);

    scales[r] = 1.0f / _pss._photos[images[r]].getScale(coord, 0);
  }

  float bestscore = 0.0f;
  std::vector<int> inout;
  inout.resize(inum, 1);
  inout[index] = 0;
  uimages.push_back(index);

  for (int t = 1; t < std::min(_tau, inum); ++t) {
    // Add the best new image
    int ansid = -1;
    float ansscore = -INT_MAX/2;
    for (int j = 0; j < inum; ++j) {
      if (inout[j] == 0) continue;

      float score = 0.0f;
      for (int k = 0; k < (int)uimages.size(); ++k) {
        const int iid = uimages[k];
        score += angleScore(rays[j], rays[iid]) * scales[j] * scales[iid];
      }

      if (ansscore < score) {
        ansscore = score;
        ansid = j;
      }
    }

    if (ansid == -1) {
      std::cerr << "Impossible 2 in compureScore" << std::endl;      exit (1);
    }

    inout[ansid] = 0;
    uimages.push_back(ansid);
    bestscore += ansscore;
  }
  return bestscore;
}

void CMVS::CBundle::writeVis(void) {
  std::ofstream ofstr;
  char buffer[1024];
  sprintf(buffer, "%svis.dat", _prefix.c_str());

  ofstr.open(buffer);
  ofstr << "VISDATA" << std::endl;
  ofstr << _cnum << std::endl;

  int numer = 0;
  int denom = 0;

  for (int c = 0; c < _cnum; ++c) {
    if (_removed[c]) {
      ofstr << c << ' ' << 0 << std::endl;
    } else {
      ofstr << c << ' ' << (int)_neighbors[c].size() << "  ";
      for (int i = 0; i < (int)_neighbors[c].size(); ++i) {
        ofstr << _neighbors[c][i] << ' ';
      }
      ofstr << std::endl;

      numer += (int)_neighbors[c].size();
      ++denom;
    }
  }
  ofstr.close();

  std::cerr << numer / (float)denom << " images in vis on the average" << std::endl;
}

void CMVS::CBundle::writeCameraCenters(void) {
  for (int i = 0; i < (int)_timages.size(); ++i) {
    char buffer[1024];
    sprintf(buffer, "%scenters-%04d.ply", _prefix.c_str(), i);

    std::ofstream ofstr;
    ofstr.open(buffer);
    ofstr << "ply" << std::endl
          << "format ascii 1.0" << std::endl
          << "element vertex " << (int)_timages[i].size() << std::endl
          << "property float x" << std::endl
          << "property float y" << std::endl
          << "property float z" << std::endl
          << "end_header" << std::endl;

    for (int j = 0; j < (int)_timages[i].size(); ++j) {
      const int c = _timages[i][j];
      ofstr << _pss._photos[c].OpticalCenter()[0] << ' '
            << _pss._photos[c].OpticalCenter()[1] << ' '
            << _pss._photos[c].OpticalCenter()[2] << std::endl;
    }
    ofstr.close();
  }

  {
    char buffer[1024];
    sprintf(buffer, "%scenters-all.ply", _prefix.c_str());
    std::ofstream ofstr;
    ofstr.open(buffer);

    ofstr << "ply" << std::endl
          << "format ascii 1.0" << std::endl
          << "element vertex " << _cnum << std::endl
          << "property float x" << std::endl
          << "property float y" << std::endl
          << "property float z" << std::endl
          << "property uchar diffuse_red" << std::endl
          << "property uchar diffuse_green" << std::endl
          << "property uchar diffuse_blue" << std::endl
          << "end_header" << std::endl;

    for (int c = 0; c < _cnum; ++c) {
      ofstr << _pss._photos[c].OpticalCenter()[0] << ' '
            << _pss._photos[c].OpticalCenter()[1] << ' '
            << _pss._photos[c].OpticalCenter()[2] << ' ';

      if (_removed[c]) {
        ofstr << "0 255 0" << std::endl;
      } else {
        ofstr << "255 0 255" << std::endl;
      }
    }
    ofstr.close();
  }
}

void CMVS::CBundle::writeGroups(void) {
  char buffer[1024];
  sprintf(buffer, "%sske.dat", _prefix.c_str());
  std::ofstream ofstr;
  ofstr.open(buffer);
  ofstr << "SKE" << std::endl
        << _cnum << ' ' << (int)_timages.size() << std::endl;

  for (int c = 0; c < (int)_timages.size(); ++c) {
    ofstr << (int)_timages[c].size() << ' ' << (int)_oimages[c].size() << std::endl;

    for (int i = 0; i < (int)_timages[c].size(); ++i) {
      ofstr << _timages[c][i] << ' ';
    }
    ofstr << std::endl;
    for (int i = 0; i < (int)_oimages[c].size(); ++i) {
      ofstr << _oimages[c][i] << ' ';
    }
    ofstr << std::endl;
  }
}

void CMVS::CBundle::startTimer(void) {
  time(&_tv);
}

time_t CMVS::CBundle::curTimer(void) {
  time(&_curtime);
  return _tv - _curtime;
}

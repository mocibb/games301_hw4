#ifndef TUTTE_H_
#define TUTTE_H_

#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <TinyAD/Scalar.hh>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"
#include "../PolyMesh/include/PolyMesh/PolyMesh_Base.h"

// #define NUMERIC_CHECK

enum TUTTE_WEIGHT {
  UNIFORM,
  FLOATER,
  MEANVALUE
};

class TutteEmbeding {
  using MVert = acamcad::polymesh::MVert;
public:
  TutteEmbeding(acamcad::polymesh::PolyMesh *xyzMesh,
                acamcad::polymesh::PolyMesh *uvMesh, 
                TUTTE_WEIGHT weightMethod_ = UNIFORM);

  void Solve();

  std::vector<size_t> GetBoundary();

protected:
  void ComputeBoundary();
  void ComputeBoundaryUV(std::unordered_map<int, std::vector<double>>& vert2uv);

  inline void UniformWeight(MVert* vi, std::vector<MVert *> vjs, std::unordered_map<int, double>& weights);
  inline void FloaterWeight(MVert* vi, std::vector<MVert *> vjs, std::unordered_map<int, double>& weights);
  inline void MeanValueWeight(MVert* vi, std::vector<MVert *> vjs, std::unordered_map<int, double>& weights);

  acamcad::polymesh::PolyMesh *xyzMesh_;
  acamcad::polymesh::PolyMesh *uvMesh_;

  std::vector<size_t> boundary_;

  TUTTE_WEIGHT weightMethod_;
};

#endif
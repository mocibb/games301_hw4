#ifndef BFF_H_
#define BFF_H_

#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"
#include "../PolyMesh/include/PolyMesh/PolyMesh_Base.h"

#define NUMERIC_CHECK

class BoundaryFlattenFirst {
  using MatrixMapType = std::unordered_map<size_t, Eigen::Matrix2d>;

public:
  BoundaryFlattenFirst(acamcad::polymesh::PolyMesh *xyzMesh,
                       acamcad::polymesh::PolyMesh *uvMesh,
                       int boundaryType);   

  void Solve();

protected:
  void PreCompute();

  bool extendHarmonic(const  Eigen::VectorXd& g,  Eigen::VectorXd& h);

  void closeLengths(const Eigen::VectorXd& lstar, const Eigen::MatrixXd& Ttilde, Eigen::VectorXd& ltilde) const;

  double computeTargetBoundaryLengths(const Eigen::VectorXd& u, Eigen::VectorXd& lstar) const;

  double computeTargetDualBoundaryLengths(const Eigen::VectorXd& lstar, Eigen::VectorXd& ldual) const;

  bool convertNeumannToDirichlet(const Eigen::VectorXd& phi, const Eigen::VectorXd& h, Eigen::VectorXd& g);

  acamcad::polymesh::PolyMesh *xyzMesh_;
  acamcad::polymesh::PolyMesh *uvMesh_;

  Eigen::SparseMatrix<double> A_;
  Eigen::SparseMatrix<double> Aii_;
  Eigen::SparseMatrix<double> Aib_;
  Eigen::SparseMatrix<double> Abb_;

  Eigen::VectorXd K_;
  Eigen::VectorXd k_;

  size_t N_;
  size_t iN_;
  size_t bN_;

  int boundaryType_;

  std::unordered_map<size_t, size_t> vertex2index_;

  std::vector<size_t> boundary_;

  Eigen::VectorXd boundary_length_;
};

#endif
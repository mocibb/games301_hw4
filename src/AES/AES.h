#ifndef AES_H_
#define AES_H_

#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <TinyAD/Scalar.hh>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"
#include "../PolyMesh/include/PolyMesh/PolyMesh_Base.h"

// #define NUMERIC_CHECK

#include "./AESUtil.h"

typedef Eigen::Vector<double, 6> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;

enum DeformationType { ARAP, SD };

class DeformationEnergy {
public:
  DeformationEnergy(DeformationType dtype) : dtype_(dtype), H_proj_(true) {}

  DeformationEnergy(const DeformationEnergy &energy) = default;

  double Cost(const Eigen::Matrix2d &source, const Eigen::Vector2d &ui,
              const Eigen::Vector2d &uj, const Eigen::Vector2d &uk);

  void ComputeDerivative(const Eigen::Matrix2d &source,
                         const Eigen::Vector2d &ui, const Eigen::Vector2d &uj,
                         const Eigen::Vector2d &uk, Vector6d &g, Matrix6d &H);

private:
  DeformationType dtype_;
  bool H_proj_;

#ifdef NUMERIC_CHECK
  static constexpr bool numeric_check_ = true;
#else
  static constexpr bool numeric_check_ = false;
#endif
  static constexpr bool output_detail_ = false;
};

class Energy2DSystem {
  using MatrixMapType = std::unordered_map<size_t, Eigen::Matrix2d>;

public:
  Energy2DSystem(acamcad::polymesh::PolyMesh *xyzMesh,
                 acamcad::polymesh::PolyMesh *uvMesh,
                 const DeformationEnergy &energy);

  void Solve();

protected:
  double Cost(const Eigen::VectorXd &x0);

  void FillGradientAndHessian(const Eigen::VectorXd &x0, Eigen::VectorXd &g,
                              Eigen::SparseMatrix<double> &H);

  void PreCompute();

  static constexpr int max_iters = 200;
  static constexpr double convergence_eps = 1e-2;

  acamcad::polymesh::PolyMesh *xyzMesh_;
  acamcad::polymesh::PolyMesh *uvMesh_;
  Eigen::VectorXd uv_;
  DeformationEnergy energy_;
  MatrixMapType sourceCoordinates_;
};

#endif
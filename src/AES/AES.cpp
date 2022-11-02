#include "./AES.h"
#include "../PolyMesh/include/Math/MPoint3.h"

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <TinyAD/Operations/SVD.hh>
#include <TinyAD/Utils/HessianProjection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>

using namespace acamcad::polymesh;

template <> Eigen::Vector<double, 6> warpData<double>(const Vector6d &d) {
  return d;
}

template <>
Eigen::Vector<TinyAD::Double<6>, 6>
warpData<TinyAD::Double<6>>(const Vector6d &d) {
  return TinyAD::Double<6>::make_active(d);
}

double DeformationEnergy::Cost(const Eigen::Matrix2d &source,
                               const Eigen::Vector2d &ui,
                               const Eigen::Vector2d &uj,
                               const Eigen::Vector2d &uk) {

  typedef std::conditional<numeric_check_, TinyAD::Double<6>, double>::type T;

  Eigen::Matrix2<T> target = createTarget<T>(ui, uj, uk);
  Eigen::Matrix2<T> F = target * source.inverse();

  // TODO: area
  T cost = 0.5 * source.determinant() * (F.squaredNorm() + F.inverse().squaredNorm());

#ifdef NUMERIC_CHECK
  if (output_detail_) {
    std::cout << "cost[n] = " << cost.val << std::endl;
    std::cout << "grad[n] = " << cost.grad.transpose() << std::endl;
    std::cout << "Hess[n] = " << cost.Hess << std::endl;
  }

  return cost.val;

#else
  if (output_detail_) {
    std::cout << "cost = " << cost << std::endl;
  }
  return cost;
#endif
}

void DeformationEnergy::ComputeDerivative(const Eigen::Matrix2d &source,
                                          const Eigen::Vector2d &ui,
                                          const Eigen::Vector2d &uj,
                                          const Eigen::Vector2d &uk,
                                          Vector6d &g, Matrix6d &H) {
  typedef std::conditional<numeric_check_, TinyAD::Double<4>, double>::type T;
  Eigen::Matrix2d target;
  target << uj - ui, uk - ui;

  Eigen::Matrix2d sourceInv = source.inverse();
  Eigen::Matrix2d F = target * sourceInv;

  Eigen::Matrix<double, 4, 6> pFpX;
  pFpX << -sourceInv(0, 0) - sourceInv(1, 0), 0, sourceInv(0, 0), 0,
      sourceInv(1, 0), 0, 0, -sourceInv(0, 0) - sourceInv(1, 0), 0,
      sourceInv(0, 0), 0, sourceInv(1, 0), -sourceInv(0, 1) - sourceInv(1, 1),
      0, sourceInv(0, 1), 0, sourceInv(1, 1), 0, 0,
      -sourceInv(0, 1) - sourceInv(1, 1), 0, sourceInv(0, 1), 0,
      sourceInv(1, 1);

  Eigen::Vector4d pCpF;
  Eigen::Matrix4d p2Cp2F;

#ifdef NUMERIC_CHECK
  Eigen::Vector4<T> vv = T::make_active({F(0, 0), F(1, 0), F(0, 1), F(1, 1)});
  Eigen::Matrix2<T> FF;
  FF << vv[0], vv[2], vv[1], vv[3];

  auto cost = 0.5 * (FF.squaredNorm() + FF.inverse().squaredNorm());

  if (H_proj_) {
    TinyAD::project_positive_definite(cost.Hess,
                                      TinyAD::default_hessian_projection_eps);
  }

  pCpF = cost.grad;
  p2Cp2F = cost.Hess;

  if (output_detail_) {
    std::cout << "pCpF2[n] = " << cost.grad << std::endl;
    std::cout << "p2Cp2F2[n] = " << cost.Hess << std::endl;
  }
#else
  // svd
  Eigen::Matrix2d U, V;
  Eigen::Vector2d sigma;
  TinyAD::svd(F, U, sigma, V);

  Eigen::Matrix2d twist;
  twist << 0, -1, 1, 0;
  Eigen::Matrix2d flip;
  flip << 0, 1, 1, 0;

  if (dtype_ == SD) {
    double I2 = F.squaredNorm();
    double I3 = F.determinant();
    double lam1 = 1 + 3 / (sigma[0] * sigma[0] * sigma[0] * sigma[0]);
    double lam2 = 1 + 3 / (sigma[1] * sigma[1] * sigma[1] * sigma[1]);
    double lam3 = 1 + 1 / (I3 * I3) - I2 / (I3 * I3 * I3);
    double lam4 = 1 + 1 / (I3 * I3) + I2 / (I3 * I3 * I3);

    if (H_proj_) {
      lam3 = std::max(lam3, 0.0);
      lam4 = std::max(lam4, 0.0);
    }

    Eigen::Matrix2d D1 = U * Eigen::Vector2d(1, 0).asDiagonal() * V.transpose();
    Eigen::Matrix2d D2 = U * Eigen::Vector2d(0, 1).asDiagonal() * V.transpose();
    Eigen::Matrix2d L = 1 / sqrt(2) * U * flip * V.transpose();
    Eigen::Matrix2d T = 1 / sqrt(2) * U * twist * V.transpose();

    Eigen::Matrix2d G = twist * F * twist.transpose();
    pCpF =
        (1 + 1 / (I3 * I3)) * Vectorize(F) - I2 / (I3 * I3 * I3) * Vectorize(G);

    //
    p2Cp2F = lam1 * Vectorize(D1) * Vectorize(D1).transpose() +
             lam2 * Vectorize(D2) * Vectorize(D2).transpose() +
             lam3 * Vectorize(L) * Vectorize(L).transpose() +
             lam4 * Vectorize(T) * Vectorize(T).transpose();

    if (output_detail_) {
      std::cout << "pCpF2 = " << pCpF.transpose() << std::endl;
      std::cout << "p2Cp2F2 = " << p2Cp2F << std::endl;
    }
  }
#endif

  g = source.determinant() * pFpX.transpose() * pCpF;
  H = source.determinant() * pFpX.transpose() * p2Cp2F * pFpX;

  if (output_detail_) {
    std::cout << "g = " << g.transpose() << std::endl;
    std::cout << "H = " << H << std::endl;
  }
}

Energy2DSystem::Energy2DSystem(acamcad::polymesh::PolyMesh *xyzMesh,
                               acamcad::polymesh::PolyMesh *uvMesh,
                               const DeformationEnergy &energy)
    : xyzMesh_(xyzMesh), uvMesh_(uvMesh), uv_(2 * uvMesh_->vertices().size()),
      energy_(energy) {
  PreCompute();
}

void Energy2DSystem::Solve() {

  Eigen::SparseMatrix<double> H;
  Eigen::VectorXd x = uv_;
  Eigen::VectorXd g;

  Eigen::SimplicialLDLT<decltype(H)> solver;

  bool analyzePattern = false;

  double oldCost = cost(x);

  for (int i = 0; i < max_iters; ++i) {

    fillGradientAndHessian(x, g, H);

    if (!analyzePattern) {
      solver.analyzePattern(H);
      analyzePattern = true;
    }

    solver.factorize(H);
    Eigen::VectorXd d = solver.solve(-g);

    if (TinyAD::newton_decrement(d, g) < convergence_eps) {
      std::cout << "Locally convergenced..." << std::endl;
      break;
    }

    // code copy from TinyAD::line_search
    double s = 1.0;
    constexpr double shrink = 0.8;
    constexpr double armijoConst = 1e-4;

    Eigen::VectorXd newX;
    for (int j = 0; j < 64; j++) {
      newX = x + s * d;
      const double newCost = cost(newX);
      if (newCost <= oldCost + armijoConst * s * d.dot(g)) {
        x = newX;
        std::cout << "iter " << i << ": " << "old cost = " << oldCost << ", new cost = " << newCost << std::endl;
        oldCost = newCost;
        break;
      }

      s *= shrink;
    }
  }

  for (auto vert : uvMesh_->vertices()) {
    vert->setPosition(x[2 * vert->index()], x[2 * vert->index() + 1], 0);
  }
}

void Energy2DSystem::PreCompute() {
  int n_vertices = uvMesh_->vertices().size();

  for (auto vert : uvMesh_->vertices()) {
    eigen_assert(vert->index() >= 0 && vert->index() < n_vertices);
    uv_[2 * vert->index()] = vert->x();
    uv_[2 * vert->index() + 1] = vert->y();
  }

  for (auto face : xyzMesh_->polyfaces()) {
    auto vertices = xyzMesh_->polygonVertices(face);

    auto xi = vertices[0];
    auto xj = vertices[1];
    auto xk = vertices[2];

    auto vij = xj->position() - xi->position();
    auto vik = xk->position() - xi->position();

    double cross_norm = norm(vij.cross(vik));
    Eigen::Matrix2d X;

    X << Eigen::Vector2d(norm(vij), 0),
        Eigen::Vector2d(vij.dot(vik), cross_norm) / norm(vij);

    sourceCoordinates_.emplace(face->index(), X);
  }
}

double Energy2DSystem::cost(const Eigen::VectorXd &x0) {
  double total = 0;
  for (auto face : xyzMesh_->polyfaces()) {
    auto vertices = xyzMesh_->polygonVertices(face);
    Vector6d gl;
    Matrix6d Hl;

    size_t i = vertices[0]->index();
    size_t j = vertices[1]->index();
    size_t k = vertices[2]->index();
    size_t toGlobals[6] = {2 * i,     2 * i + 1, 2 * j,
                           2 * j + 1, 2 * k,     2 * k + 1};

    Eigen::Vector2d ui = x0.segment(i * 2, 2);
    Eigen::Vector2d uj = x0.segment(j * 2, 2);
    Eigen::Vector2d uk = x0.segment(k * 2, 2);

    total += energy_.Cost(sourceCoordinates_.at(face->index()), ui, uj, uk);
  }

  return total;
}

void Energy2DSystem::fillGradientAndHessian(const Eigen::VectorXd &x0,
                                            Eigen::VectorXd &g,
                                            Eigen::SparseMatrix<double> &H) {

  std::vector<Eigen::Triplet<double>> H_triplets;
  g = Eigen::VectorXd::Zero(x0.size());
  H = Eigen::SparseMatrix<double>(x0.size(), x0.size());

  for (auto face : xyzMesh_->polyfaces()) {
    auto vertices = xyzMesh_->polygonVertices(face);
    Vector6d gl;
    Matrix6d Hl;

    size_t i = vertices[0]->index();
    size_t j = vertices[1]->index();
    size_t k = vertices[2]->index();
    size_t toGlobals[6] = {2 * i,     2 * i + 1, 2 * j,
                           2 * j + 1, 2 * k,     2 * k + 1};

    Eigen::Vector2d ui = x0.segment(i * 2, 2);
    Eigen::Vector2d uj = x0.segment(j * 2, 2);
    Eigen::Vector2d uk = x0.segment(k * 2, 2);

    auto cost = energy_.Cost(sourceCoordinates_.at(face->index()), ui, uj, uk);

    energy_.ComputeDerivative(sourceCoordinates_.at(face->index()), ui, uj, uk,
                              gl, Hl);

    // throw new std::runtime_error("check derivative");

    for (size_t n = 0; n < 6; n++) {
      g[toGlobals[n]] += gl[n];
      for (size_t m = 0; m < 6; m++) {
        H_triplets.push_back(
            Eigen::Triplet<double>(toGlobals[n], toGlobals[m], Hl(n, m)));
      }
    }
  }

  H.setFromTriplets(H_triplets.begin(), H_triplets.end());

  // std::cout << "finish H" << std::endl;
} 

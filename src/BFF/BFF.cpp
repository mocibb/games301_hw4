#include "./BFF.h"
#include "../PolyMesh/include/Math/MPoint3.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

using namespace acamcad::polymesh;
using acamcad::MVector3;

#define ERR throw new std::runtime_error("error")

BoundaryFlattenFirst::BoundaryFlattenFirst(acamcad::polymesh::PolyMesh *xyzMesh,
                                           acamcad::polymesh::PolyMesh *uvMesh,
                                           int boundaryType)
    : xyzMesh_(xyzMesh), uvMesh_(uvMesh), boundaryType_(boundaryType),
      A_(xyzMesh_->numVertices(), xyzMesh_->numVertices())
{
  PreCompute();
}

void BoundaryFlattenFirst::PreCompute()
{
  /// 计算边界
  boundary_.clear();

  for (const auto he : xyzMesh_->halfEdges())
  {
    if (xyzMesh_->isBoundary(he) &&
        std::find(boundary_.begin(), boundary_.end(), he->toVertex()->index()) ==
            boundary_.end())
    {
      {
        auto *cur = he;
        int cnt = 0;
        do
        {
          boundary_.push_back(cur->toVertex()->index());
          // std::cout << "v = " << cur->toVertex()->index() << " " << cur->toVertex()->position()[0] << " " << cur->toVertex()->position()[1] << " " << cur->toVertex()->position()[2] << std::endl;

          cur = cur->next();
          cnt++;
        } while (cur != he);
        // std::cout << "cnt = " << cnt << std::endl;
      }
    }
  }

  /// 边界长度
  boundary_length_ = Eigen::VectorXd::Zero(boundary_.size());
  for (size_t k = 0; k < boundary_.size(); k++)
  {
    size_t i = boundary_[k];
    size_t j = boundary_[(k + boundary_.size() - 1) % boundary_.size()];

    boundary_length_(k) = (xyzMesh_->vert(i)->position() - xyzMesh_->vert(j)->position()).norm();
    // std::cout << "l[" << k << "]=" << boundary_length_(k) << std::endl;
  }

  /// 调整index，把边界放到最后
  auto vert = xyzMesh_->vertices();
  N_ = vert.size();
  bN_ = boundary_.size();
  iN_ = N_ - bN_;

  auto bend = boundary_.end();
  int iidx = 0;
  int bidx = iN_;
  for (auto v : vert)
  {
    if (std::find(boundary_.begin(), bend, v->index()) == bend)
    {
      vertex2index_.emplace(v->index(), iidx);
      iidx++;
    }
  }

  for (size_t k = 0; k < boundary_.size(); k++)
  {
    vertex2index_.emplace(boundary_[k], bidx++);
  }

  /// 计算 cotan matrix
  std::vector<Eigen::Triplet<double>> ATriplets;

  for (auto f : xyzMesh_->polyfaces())
  {
    auto he = f->halfEdge();
    do
    {
      double w = 0;

      int i = vertex2index_[he->fromVertex()->index()];
      int j = vertex2index_[he->toVertex()->index()];

      if (!xyzMesh_->isBoundary(he))
      {
        const auto &a = he->fromVertex()->position();
        const auto &b = he->next()->fromVertex()->position();
        const auto &c = he->prev()->fromVertex()->position();

        MVector3 u = a - c;
        MVector3 v = b - c;

        double cotan = u.dot(v) / u.cross(v).norm();
        if (std::isinf(cotan) || std::isnan(cotan))
        {
          cotan = 0.0;
        }

        w = 0.5 * cotan;
      }

      ATriplets.emplace_back(i, i, w);
      ATriplets.emplace_back(j, j, w);
      ATriplets.emplace_back(i, j, -w);
      ATriplets.emplace_back(j, i, -w);

      he = he->next();
    } while (he != f->halfEdge());
  }

  Eigen::SparseMatrix<double> id(N_, N_);
  id.setIdentity();
  A_.setFromTriplets(ATriplets.begin(), ATriplets.end());
  A_ += id * 1e-8;

  Aii_ = A_.block(0, 0, iN_, iN_);
  Aib_ = A_.block(0, iN_, iN_, bN_);
  Abb_ = A_.block(iN_, iN_, bN_, bN_);

  A_solver_.analyzePattern(A_);
  A_solver_.factorize(A_);

  Aii_solver_.analyzePattern(Aii_);
  Aii_solver_.factorize(Aii_);

  /// 计算曲率
  K_ = Eigen::VectorXd(iN_);

  for (auto vi : vert)
  {
    if (std::find(boundary_.begin(), bend, vi->index()) == bend)
    {
      int i = vertex2index_[vi->index()];
      double angleDefect = 2 * M_PI;

      auto vjs = xyzMesh_->vertAdjacentVertices(vi);
      int n = vjs.size();

      for (int j = 0; j < n; j++)
      {
        // 向量节点i到节点j
        MVector3 ej = vjs[j]->position() - vi->position();
        // 向量节点i到节点j+1
        MVector3 ejp1 = vjs[(j + 1) % n]->position() - vi->position();

        angleDefect -= vectorAngle(ej, ejp1);
      }

      K_(i) = angleDefect;
      // std::cout << "K(" << i << ") = " << K_(i) << std::endl;
    }
  }

  k_ = Eigen::VectorXd(bN_);
  for (size_t k = 0; k < boundary_.size(); k++)
  {
    size_t i = vertex2index_.at(boundary_[k]) - iN_;
    double exteriorAngle = M_PI;

    auto vi = xyzMesh_->vert(boundary_[k]);

    auto faces = xyzMesh_->vertAdjacentPolygon(vi);

    for (auto f : faces)
    {
      auto vjs = xyzMesh_->polygonVertices(f);
      int n = vjs.size();

      MVector3 ej[2];
      int idx = 0;
      for (int j = 0; j < n; j++)
      {
        if (vjs[j] != vi)
        {
          ej[idx++] = vjs[j]->position() - vi->position();
        }
      }
      exteriorAngle -= vectorAngle(ej[0], ej[1]);
    }

    k_(i) = exteriorAngle;
    // std::cout << "k(" << i << ") = " << k_(i) << std::endl;
  }

  if (boundaryType_ == 2) {
    int step = bN_>>2;
    target_ = Eigen::VectorXd::Zero(bN_);
    target_(0) = M_PI/2;
    target_(step) = M_PI/2;
    target_(step*2) = M_PI/2;
    target_(step*3) = M_PI/2;
    // std::cout << "target = " << target_.transpose() << std::endl;
  }
}

void BoundaryFlattenFirst::Solve()
{
  Eigen::VectorXd u = Eigen::VectorXd::Zero(bN_);
  Eigen::VectorXd ktilde = Eigen::VectorXd::Zero(bN_);
  Eigen::VectorXd lstar;

  if (boundaryType_ == 1) {
    /// 圆盘的曲率是1，根据曲率计算编辑尺度因子u
    // 因为尺度因子u未知，所以ktilde
    for (int iter = 0; iter < 10; iter++)
    {
      computeTargetBoundaryLengths(u, lstar);

      Eigen::VectorXd ldual;
      double L = computeTargetDualBoundaryLengths(lstar, ldual);
      // std::cout << "L = " << L << std::endl;

      // set ktilde proportional to the most recent dual lengths
      for (size_t k = 0; k < boundary_.size(); k++)
      {
        size_t i = vertex2index_.at(boundary_[k]) - iN_;

        ktilde(i) = 2 * M_PI * ldual(i) / L;
      }

      // // compute target scale factors
      if (!convertNeumannToDirichlet(-K_, k_ - ktilde, u))
      {
        std::cout << "convertNeumannToDirichlet failed" << std::endl;
        throw new std::runtime_error("convertNeumannToDirichlet");
      };
    }
  } else if (boundaryType_ == 2) {
    ktilde = target_;
    if (!convertNeumannToDirichlet(-K_, k_ - target_, u)) {
      std::cout << "convertNeumannToDirichlet failed" << std::endl;
      throw new std::runtime_error("convertNeumannToDirichlet");
    }
  }

  // std::cout << "u = " << u.transpose() << std::endl;

  // 
  double finalL = computeTargetBoundaryLengths(u, lstar);
  // std::cout << "finalL = " << finalL << std::endl;
  // std::cout << "lstar = " << lstar.transpose() << std::endl;

  Eigen::VectorXd gammaRe, gammaIm;
  constructBestFitCurve(lstar, ktilde, gammaRe, gammaIm);
	// for (size_t i = 0; i < bN_; i++) {
	// 	std::cout << "gamma = " << gammaRe(i) << ", " << gammaIm(i) << std::endl;
	// }

  // extend
  Eigen::VectorXd flatteningRe, flatteningIm;
  if (!extendCurve(gammaRe, gammaIm, flatteningRe, flatteningIm)) {
    std::cout << "extendCurve failed" << std::endl;
    throw new std::runtime_error("extendCurve");
  }

	// for (size_t i = 0; i < N_; i++) {
	// 	std::cout << "flattening = " << flatteningRe(i) << ", " << flatteningIm(i) << std::endl;
	// }

  uvMesh_->clear();
  
  // 映射后的坐标用来构建uvMesh_
  std::unordered_map<size_t, MVert*> uv_vertices;

  for (auto vert : xyzMesh_->vertices()) {
    size_t i = vertex2index_.at(vert->index());
    uv_vertices.emplace( vert->index(), uvMesh_->addVertex(flatteningRe[i], flatteningIm[i], 0) );
  }

  for (const auto &fh : xyzMesh_->polyfaces()) {
    auto xyz_vertices = xyzMesh_->polygonVertices(fh);
    std::vector<MVert *> f_vertices;
    for (MVert *v : xyz_vertices) {
      f_vertices.push_back(uv_vertices.at(v->index()));
    }

    uvMesh_->addPolyFace(f_vertices);
  }
}

bool BoundaryFlattenFirst::extendCurve(const Eigen::VectorXd &gammaRe, const Eigen::VectorXd &gammaIm, Eigen::VectorXd &a, Eigen::VectorXd &b)
{
  std::cout << "extendCurve" << std::endl;
  // extend real component of gamma
  if (!extendHarmonic(gammaRe, a))
    return false;

  // extend imaginary component of gamma
  if (!extendHarmonic(gammaIm, b))
    return false;

  return true;
}

void BoundaryFlattenFirst::constructBestFitCurve(const Eigen::VectorXd &lstar, const Eigen::VectorXd &ktilde, Eigen::VectorXd &gammaRe, Eigen::VectorXd &gammaIm) const
{
  std::cout << "constructBestFitCurve" << std::endl;
  // compute tangents as cumulative sum of angles phi
  double phi = 0.0;
  Eigen::MatrixXd Ttilde(2, bN_);

  for (size_t i = 0; i < boundary_.size(); i++)
  {
    phi += ktilde(i);
    Ttilde(0, i) = cos(phi);
    Ttilde(1, i) = sin(phi);
  }

  // modify target lengths lstar to ensure gamma closes
  Eigen::VectorXd ltilde;
  closeLengths(lstar, Ttilde, ltilde);

  // compute gamma as cumulative sum of products ltilde*Ttilde
  double re = 0.0;
  double im = 0.0;
  gammaRe = Eigen::VectorXd(bN_);
  gammaIm = Eigen::VectorXd(bN_);

  // std::cout << "ltilde = " << ltilde.transpose() << std::endl;
  // for (size_t i = 0; i < boundary_.size(); i++) {
  //   std::cout << "Ttilde = " << Ttilde.col(i).transpose() << std::endl;
  // }
  
  for (size_t i = 0; i < boundary_.size(); i++)
  {
    gammaRe(i) = re;
    gammaIm(i) = im;
    re += ltilde(i) * Ttilde(0, i);
    im += ltilde(i) * Ttilde(1, i);
  }
}


double BoundaryFlattenFirst::computeTargetBoundaryLengths(const Eigen::VectorXd &u, Eigen::VectorXd &lstar) const
{
  double sum = 0.0;
  lstar = Eigen::VectorXd(bN_);
  for (size_t k = 0; k < boundary_.size(); k++)
  {
    size_t i = vertex2index_.at(boundary_[k]) - iN_;
    size_t j = vertex2index_.at(boundary_[(k + boundary_.size() + 1) % boundary_.size()]) - iN_;

    lstar(i) = exp(0.5 * (u(i) + u(j))) * boundary_length_(i);
    sum += lstar(i);
  }

  return sum;
}

double BoundaryFlattenFirst::computeTargetDualBoundaryLengths(const Eigen::VectorXd &lstar, Eigen::VectorXd &ldual) const
{
  double sum = 0.0;
  ldual = Eigen::VectorXd(bN_);

  for (size_t k = 0; k < boundary_.size(); k++)
  {
    size_t i = vertex2index_.at(boundary_[k]) - iN_;
    size_t j = vertex2index_.at(boundary_[(k + boundary_.size() + 1) % boundary_.size()]) - iN_;

    ldual(j) = 0.5 * (lstar(i) + lstar(j));
    sum += ldual(j);
  }

  return sum;
}

void BoundaryFlattenFirst::closeLengths(const Eigen::VectorXd &lstar, const Eigen::MatrixXd &Ttilde, Eigen::VectorXd &ltilde) const
{
  std::cout << "closeLengths" << std::endl;
  int N = lstar.rows();
  Eigen::SparseMatrix<double> A(N + 2, N + 2);
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N + 2);

  std::vector<Eigen::Triplet<double>> ATriplets;
  auto vert = xyzMesh_->vertices();

  for (size_t k = 0; k < boundary_.size(); k++)
  {
    size_t next_k = boundary_[(k + 1) % boundary_.size()];
    MVector3 ei = vert[k]->position() - vert[next_k]->position();
    size_t i = vertex2index_.at(boundary_[k]) - iN_;
    double l = ei.norm();

    ATriplets.emplace_back(i, i, 1 / l);
    ATriplets.emplace_back(i, N , Ttilde(0, i));
    ATriplets.emplace_back(i, N + 1, Ttilde(1, i));
    ATriplets.emplace_back(N, i, Ttilde(0, i));
    ATriplets.emplace_back(N + 1, i, Ttilde(1, i));

    rhs[i] = lstar[i] / l;
  }

  A.setFromTriplets(ATriplets.begin(), ATriplets.end());

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
  Eigen::VectorXd a = solver.solve(rhs);
  ltilde = a.head(N); 

}

bool BoundaryFlattenFirst::convertNeumannToDirichlet(const Eigen::VectorXd &phi, const Eigen::VectorXd &h, Eigen::VectorXd &g)
{
  Eigen::VectorXd rhs(phi.rows() + h.rows());
  rhs.head(phi.rows()) = phi;
  rhs.segment(phi.rows(), h.rows()) = -h;

  Eigen::VectorXd a = A_solver_.solve(rhs);

  g = a.segment(phi.rows(), h.rows());
  return true;
}

bool BoundaryFlattenFirst::extendHarmonic(const Eigen::VectorXd &g, Eigen::VectorXd &h)
{
  Eigen::VectorXd rhs = -(Aib_ * g);
  Eigen::VectorXd a = Aii_solver_.solve(rhs);
  
  h.resize(a.rows() + g.rows());
  
  h.head(a.rows()) = a;
  h.segment(a.rows(), g.rows()) = g;

  return true;
}

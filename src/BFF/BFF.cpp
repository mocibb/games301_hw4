#include "./BFF.h"
#include "../PolyMesh/include/Math/MPoint3.h"
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

using namespace acamcad::polymesh;
using acamcad::MVector3;

BoundaryFlattenFirst::BoundaryFlattenFirst(acamcad::polymesh::PolyMesh *xyzMesh,
                                           acamcad::polymesh::PolyMesh *uvMesh,
                                           int boundaryType)
    : xyzMesh_(xyzMesh), uvMesh_(uvMesh), boundaryType_(boundaryType),
    A_(xyzMesh_->numVertices(), xyzMesh_->numVertices()) {
  PreCompute();
}

void BoundaryFlattenFirst::PreCompute() {
  // 计算边界
  // 保存边界节点的index
  boundary_.clear();

  for (const auto he : xyzMesh_->halfEdges()) {
    if (xyzMesh_->isBoundary(he) &&
        std::find(boundary_.begin(), boundary_.end(), he->toVertex()->index()) ==
            boundary_.end()) {
      {
        auto* cur = he;
        do {
          boundary_.push_back(cur->toVertex()->index());
          cur = cur->next();
        } while (cur != he);
      }
    }
  }

  std::cout << "boundary = " << boundary_.size() << std::endl;

  std::cout << "++++++++++++++++" << std::endl;
  // 边界长度
  boundary_length_ = Eigen::VectorXd::Zero(boundary_.size());
  for (size_t k = 0; k < boundary_.size(); k++) {
    size_t i = boundary_[k];
    size_t j = boundary_[(k+boundary_.size()-1)%boundary_.size()];

    boundary_length_(k) = (xyzMesh_->vert(i)->position()-xyzMesh_->vert(j)->position()).norm();
    // std::cout << boundary_length_(k) << std::endl;
	}
  std::cout << "----------------" << std::endl;

  // 调整index，把边界放到最后
  auto vert = xyzMesh_->vertices();
  N_ = vert.size();
  iN_ = N_ - boundary_.size();
  bN_ = N_ - iN_;
  auto bend = boundary_.end();
  int iidx = 0;
  int bidx = iN_;
  for (auto v : vert) {
    if (std::find(boundary_.begin(), boundary_.end(), v->index()) == bend) {
      vertex2index_.emplace(v->index(), iidx);
      iidx++;
    } else {
      vertex2index_.emplace(v->index(), bidx);
      bidx++;
    }
  }

  // 计算 cotan matrix
  Eigen::SparseMatrix<double> A(N_, N_);
  std::vector<Eigen::Triplet<double>> ATriplets;

  for (auto f : xyzMesh_->polyfaces()) {
    auto he = f->halfEdge();
    do {
      int i = vertex2index_[he->fromVertex()->index()];
      int j = vertex2index_[he->toVertex()->index()];

      double w = 0;

      if (!xyzMesh_->isBoundary(he)) {
        const auto& a = he->fromVertex()->position();
        const auto& b = he->next()->fromVertex()->position();
        const auto& c = he->prev()->fromVertex()->position();

        MVector3 u = a - c;
        MVector3 v = b - c;

        double cotan = u.dot(v)/u.cross(v).norm();
        if (std::isinf(cotan) || std::isnan(cotan)) {
          cotan = 0.0;
        }

        w = 0.5*cotan;
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
	A_ += id*1e-8;
  
  Aii_ = A_.block(0, 0, iN_, iN_);
  Aib_ = A_.block(0, iN_, iN_, bN_);
  Abb_ = A_.block(iN_, iN_, bN_, bN_);


  // 计算曲率
  K_ = Eigen::VectorXd(iN_);

	for (auto vi : vert) {
		if (std::find(boundary_.begin(), boundary_.end(), vi->index()) == bend) {
			int i = vertex2index_[vi->index()];
			double angleDefect = 2*M_PI;

      auto vjs = xyzMesh_->vertAdjacentVertices(vi);
      int n = vjs.size();

      for (int j = 0; j < n; j++) {
        // 向量节点i到节点j
        MVector3 ej = vjs[j]->position() - vi->position();
        // 向量节点i到节点j+1
        MVector3 ejp1 = vjs[(j + 1) % n]->position() - vi->position();

        angleDefect -= vectorAngle(ej, ejp1);
      }

      K_(i) = angleDefect;
		}
	}

	k_ =  Eigen::VectorXd(bN_);
	for (size_t k = 0; k < boundary_.size(); k++) {
		size_t i = vertex2index_.at(boundary_[k])-iN_;
    double exteriorAngle = M_PI;

    auto vi = xyzMesh_->vert(boundary_[k]);

    auto faces = xyzMesh_->vertAdjacentPolygon(vi);

    for (auto f : faces) {
      auto vjs = xyzMesh_->polygonVertices(f);
      int n = vjs.size();

      MVector3 ej[2];
      int idx = 0;
      for (int j = 0; j < n; j++) {
        if (vjs[j] != vi) {
          ej[idx++] = vjs[j]->position() - vi->position();
        }
      }
      exteriorAngle -= vectorAngle(ej[0], ej[1]);
    }
    
    k_(i) = exteriorAngle;
	}

}

void BoundaryFlattenFirst::Solve() {
  Eigen::VectorXd u = Eigen::VectorXd::Zero(bN_);
  Eigen::VectorXd ktilde = Eigen::VectorXd::Zero(bN_);

  for (int iter = 0; iter < 10; iter++) {
		// compute target dual boundary edge lengths
		double L;
		Eigen::VectorXd lstar, ldual;
		computeTargetBoundaryLengths(u, lstar);

		L = computeTargetDualBoundaryLengths(lstar, ldual);
    std::cout << "L = " << L << std::endl;

		// set ktilde proportional to the most recent dual lengths
    for (size_t k = 0; k < boundary_.size(); k++) {

      size_t i = vertex2index_.at(boundary_[k])-iN_;

			ktilde(i) = 2*M_PI*ldual(i)/L;
		}

		// // compute target scale factors
		if (!convertNeumannToDirichlet(-K_, k_ - ktilde, u)) {
      throw new std::runtime_error("wrong");
    };
	}
}

double BoundaryFlattenFirst::computeTargetBoundaryLengths(const Eigen::VectorXd& u, Eigen::VectorXd& lstar) const
{
	double sum = 0.0;
	lstar = Eigen::VectorXd(bN_);
  for (size_t k = 0; k < boundary_.size(); k++) {
    size_t i = vertex2index_.at(boundary_[k])-iN_;
    size_t j = vertex2index_.at(boundary_[(k+boundary_.size()-1)%boundary_.size()])-iN_;

		lstar(i) = exp(0.5*(u(i) + u(j)))*boundary_length_(i);
		sum += lstar(i);
	}

	return sum;
}

double BoundaryFlattenFirst::computeTargetDualBoundaryLengths(const Eigen::VectorXd& lstar, Eigen::VectorXd& ldual) const
{
	double sum = 0.0;
	ldual = Eigen::VectorXd(bN_);

  for (size_t k = 0; k < boundary_.size(); k++) {
    size_t i = vertex2index_.at(boundary_[k])-iN_;
    size_t j = vertex2index_.at(boundary_[(k+boundary_.size()-1)%boundary_.size()])-iN_;

		ldual(j) = 0.5*(lstar(i) + lstar(j));
		sum += ldual(j);
	}

	return sum;
}


void BoundaryFlattenFirst::closeLengths(const Eigen::VectorXd& lstar, const Eigen::MatrixXd& Ttilde, Eigen::VectorXd& ltilde) const
{
  int N = lstar.rows();
  Eigen::SparseMatrix<double> A(N+2, N+2);
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N+2);

  std::vector<Eigen::Triplet<double>> ATriplets;
  auto vert = xyzMesh_->vertices();
  for (size_t i = 0; i < boundary_.size(); i++) {
    size_t next_i = boundary_[(i+1)%boundary_.size()];
    MVector3 ei = vert[i]->position() - vert[next_i]->position();
    double l = ei.norm();
      
    ATriplets.emplace_back(i, i, 1/l);
    ATriplets.emplace_back(i, N+1, Ttilde(i, 0));
    ATriplets.emplace_back(i, N+2, Ttilde(i, 1));
    ATriplets.emplace_back(N+1, i, Ttilde(i, 0));
    ATriplets.emplace_back(N+2, i, Ttilde(i, 1));
    rhs[i] = lstar[i]/l;
  }
  A.setFromTriplets(ATriplets.begin(), ATriplets.end());

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
  Eigen::VectorXd a = solver.solve(rhs);

  ltilde = a.head(N);
}


bool BoundaryFlattenFirst::convertNeumannToDirichlet(const Eigen::VectorXd& phi, const Eigen::VectorXd& h, Eigen::VectorXd& g)
{
  Eigen::VectorXd rhs(phi.rows()+h.rows());
  rhs.head(phi.rows()) = phi;
  rhs.segment(phi.rows(), h.rows()) = -h;

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A_);
  Eigen::VectorXd a = solver.solve(rhs);

	g = a.segment(phi.rows(), h.rows());
	return true;
}

bool BoundaryFlattenFirst::extendHarmonic(const  Eigen::VectorXd& g,  Eigen::VectorXd& h)
{
	Eigen::VectorXd rhs = -(Aib_*g);

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(Aii_);
  Eigen::VectorXd a = solver.solve(rhs);

  h.resize(a.rows()+g.rows());
  h.head(a.rows()) = a;
  h.segment(a.rows(), h.rows()) = g;

	return true;
}

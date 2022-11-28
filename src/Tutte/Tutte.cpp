#include "./Tutte.h"
#include "../PolyMesh/include/Math/MPoint3.h"

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <TinyAD/Operations/SVD.hh>
#include <TinyAD/Utils/HessianProjection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>

using namespace acamcad;

TutteEmbeding::TutteEmbeding(acamcad::polymesh::PolyMesh *xyzMesh,
                             acamcad::polymesh::PolyMesh *uvMesh,
                             TUTTE_WEIGHT weightMethod)
    : xyzMesh_(xyzMesh), uvMesh_(uvMesh), weightMethod_(weightMethod) {
  ComputeBoundary();
}

void TutteEmbeding::ComputeBoundary() {
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
}

void TutteEmbeding::ComputeBoundaryUV(std::unordered_map<int, std::vector<double>>& vert2uv) {
  if (weightMethod_ != MEANVALUE) {
    int i = 0;
    // 映射边界到2D， 100(半径)*exp(i theta[])
    double step = 2 * M_PI / boundary_.size();
    for (auto idx : boundary_) {
      vert2uv.emplace(idx, std::vector<double>());

      vert2uv[idx].push_back(1 * cos(i * step));
      vert2uv[idx].push_back(1 * sin(i * step));
      i++;
    }
  } else {
    //
    for (auto idx : boundary_) {
      auto vert = uvMesh_->vert(idx);
      vert2uv[idx].push_back(vert->position()[0]);
      vert2uv[idx].push_back(vert->position()[1]);
    }
    //
  }
}

inline void TutteEmbeding::UniformWeight(MVert * vi, std::vector<MVert *> vjs, std::unordered_map<int, double>& weights) {
  int n = vjs.size();
  for (int j = 0; j < n; j++) {
    weights.emplace(vjs[j]->index(), 1);
  }
}

inline void TutteEmbeding::FloaterWeight(MVert * vi, std::vector<MVert *> vjs, std::unordered_map<int, double>& weights) {
    double sum_angles = 0;
    std::vector<double> lengths;
    std::vector<double> angles;

    int n = vjs.size();

    for (int j = 0; j < n; j++) {
      // 向量节点i到节点j
      MVector3 ej = vjs[j]->position() - vi->position();
      // 向量节点i到节点j+1
      MVector3 ejp1 = vjs[(j + 1) % n]->position() - vi->position();

      lengths.push_back(norm(ej));
      angles.push_back(vectorAngle(ej, ejp1));
      sum_angles += angles.back();

      // 
      weights.emplace(vjs[j]->index(), 0);
    }

    // 从3D 映射到2D， vi设置到原点，其他点坐标为 length[j]* exp(i theta[j])
    std::vector<Eigen::Vector2d> vertices;
    double theta = 0;
    for (int j = 0; j < n; j++) {
      vertices.push_back(lengths[j] * Eigen::Vector2d(cos(theta), sin(theta)));
      theta += angles[j] / sum_angles * 2 * M_PI;
    }

    // 遍历每个相邻节点计算论文中p_{r(l)}和p_{r(l)+1}
    // p_{r(l)}和p_{r(l)+1}是过p_l和p直线和凸多边形相交的边。
    for (int l = 0; l < n; l++) {
      bool found = false;
      for (int r = 0; r < n && !found; r++) {
        // 选择不同于l的其他的顶点
        if (r == l || (r + 1) % n == l) {
          continue;
        }

        int k = (r + 1) % n;

        // 点(u1, v1)和(u2, v2)共线，那么u1*v2 = u2*v1;
        // 其中(u1, v1)是p_l - p, (u2, v2)是 t*p_{r(l)}+(1-t)*p_{r(l)+1}
        Eigen::Vector2d uu = vertices[l];
        double den = uu[0] * (vertices[r][1] - vertices[k][1]) -
                     uu[1] * (vertices[r][0] - vertices[k][0]);

        if (std::abs(den) < 1e-5) {
          continue;
        }
        double t = (vertices[k][0] * uu[1] - vertices[k][1] * uu[0]) / den;

        // 如果 t在(0, 1]间，计算barycenteric坐标
        if (t > 0 && t <= 1 + 1e-4) {
          Eigen::Vector2d v0 = vertices[r] - vertices[l];
          Eigen::Vector2d v1 = vertices[k] - vertices[l];
          Eigen::Vector2d v2 = -vertices[l];

          // Compute dot products
          double dot00 = v0.dot(v0);
          double dot01 = v0.dot(v1);
          double dot02 = v0.dot(v2);
          double dot11 = v1.dot(v1);
          double dot12 = v1.dot(v2);

          // Compute barycentric coordinates
          double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
          double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
          double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

          // 分别设置l,r,k的权重
          weights.at(vjs[l]->index()) += 1 - u - v;
          weights.at(vjs[r]->index()) += u;
          weights.at(vjs[k]->index()) += v;

          found = true;
        }
      }
    }

    // 归一化权重
    double sum_weights = 0;
    for (auto &w : weights) {
      w.second = w.second / n;
      sum_weights += w.second;
      // 权重为负，非同胚映射
      if (w.second < 0) {
        throw new std::runtime_error("weight is wrong, w = " + std::to_string(w.second));
      }
    }

    // 归一化的权重为1
    if (std::abs(sum_weights - 1) > 1e-3) {
      throw new std::runtime_error("sum_weights is sum_weights, w = " + std::to_string(sum_weights));
    }
}

inline void TutteEmbeding::MeanValueWeight(MVert * vi, std::vector<MVert *> vjs, std::unordered_map<int, double>& weights) {
    std::vector<double> angles;

    int n = vjs.size();

    for (int j = 0; j < n; j++) {
      // 向量节点i到节点j
      MVector3 ej = vjs[j]->position() - vi->position();
      // 向量节点i到节点j+1
      MVector3 ejp1 = vjs[(j + 1) % n]->position() - vi->position();
      
      // alpha[j]
      angles.push_back(vectorAngle(ej, ejp1));
    }

    double sum_weights = 0;
    for (int j = 0; j < n; j++) {
      // 向量节点i到节点j
      MVector3 ej = vjs[j]->position() - vi->position();

      double ajm1 = angles[(j-1+n)%n];
      double w = (tan(ajm1/2) + tan(angles[j]/2)) / norm(ej);
      weights.emplace(vjs[j]->index(), w);
      sum_weights += w;

      if (w < 0) {
        throw new std::runtime_error("w is negative.");
      }
    }

    for (auto &w : weights) {
      w.second = w.second / sum_weights;
    }

}



void TutteEmbeding::Solve() {
  if (boundary_.size() == 0) {
    std::cout << "ERROR: mesh is no Homeomorphic with disk." << std::endl;
    return;
  }

  using namespace acamcad;
  using acamcad::polymesh::MHalfedge;
  using acamcad::polymesh::MVert;

  // 二层map用来保存节点的权重W_ij
  std::unordered_map<int, std::unordered_map<int, double>> vertex_weights;

  std::cout << "Computing Weight..." << std::endl;

  // 计算权重
  for (const auto vi : xyzMesh_->vertices()) {
    // 边界节点不需要计算权重
    if (xyzMesh_->isBoundary(vi)) {
      continue;
    }
    // 邻居节点 
    std::vector<MVert *> vjs = xyzMesh_->vertAdjacentVertices(vi);
    
    vertex_weights.emplace(vi->index(), std::unordered_map<int, double>());
    std::unordered_map<int, double> &weights = vertex_weights.at(vi->index());

    switch (weightMethod_) {
      case UNIFORM:
        UniformWeight(vi, vjs, weights);
        break;
      case FLOATER:
        FloaterWeight(vi, vjs, weights);
        break;
      case MEANVALUE:
        MeanValueWeight(vi, vjs, weights);
        break;
    } 
  }

  std::unordered_map<int, std::vector<double>> vert2uv;
  ComputeBoundaryUV(vert2uv);

  std::cout << "Contructing matrix A..." << std::endl;
  // 构建矩阵临接权重矩阵A和b
  std::vector<Eigen::Triplet<float>> ATriplets;
  std::vector<double> u_vec, v_vec;
  Eigen::VectorXd ub, vb;

  for (const auto v : xyzMesh_->vertices()) {
    // 边界情况
    if (xyzMesh_->isBoundary(v)) {
      ATriplets.emplace_back(v->index(), v->index(), 1);
      u_vec.push_back(vert2uv[v->index()][0]);
      v_vec.push_back(vert2uv[v->index()][1]);
    } else {
      double total_weight = 0;
      const auto &weights = vertex_weights.at(v->index());

      for (const auto w : xyzMesh_->vertAdjacentVertices(v)) {
        ATriplets.emplace_back(v->index(), w->index(), weights.at(w->index()));
        total_weight += weights.at(w->index());
      }
      ATriplets.emplace_back(v->index(), v->index(), -total_weight);
      u_vec.push_back(0);
      v_vec.push_back(0);
    }
  }

  Eigen::SparseMatrix<double> A(xyzMesh_->numVertices(), xyzMesh_->numVertices());
  Eigen::Map<Eigen::VectorXd> uu(u_vec.data(), u_vec.size());
  Eigen::Map<Eigen::VectorXd> vv(v_vec.data(), v_vec.size());
  A.setFromTriplets(ATriplets.begin(), ATriplets.end());

  std::cout << "Solving Linear equation..." << std::endl;
  // LU算法求解
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
  ub = solver.solve(uu);
  vb = solver.solve(vv);

  // 
  uvMesh_->clear();
  
  // 映射后的坐标用来构建uvMesh_
  std::vector<MVert*> uv_vertices;

  for (auto vert : xyzMesh_->vertices()) {
    uv_vertices.push_back( uvMesh_->addVertex(ub[vert->index()], vb[vert->index()], 0) );
  }

  for (const auto &fh : xyzMesh_->polyfaces()) {
    auto xyz_vertices = xyzMesh_->polygonVertices(fh);
    std::vector<MVert *> f_vertices;
    for (MVert *v : xyz_vertices) {
      f_vertices.push_back(uv_vertices[v->index()]);
    }

    uvMesh_->addPolyFace(f_vertices);
  }

}


std::vector<size_t> TutteEmbeding::GetBoundary() {
  return boundary_;
}

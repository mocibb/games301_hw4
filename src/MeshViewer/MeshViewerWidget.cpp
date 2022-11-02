#include "MeshViewerWidget.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <QtCore>
#include <cfloat>
#include <unordered_map>
#include <unordered_set>
#include "../AES/AES.h"


MeshViewerWidget::MeshViewerWidget(QWidget *parent)
    : QGLViewerWidget(parent), ptMin(0, 0, 0), ptMax(0, 0, 0),
      isEnableLighting(true), isTwoSideLighting(false),
      isDrawBoundingBox(false), isDrawBoundary(false) {}

MeshViewerWidget::~MeshViewerWidget(void) {}

bool MeshViewerWidget::LoadMesh(const std::string &filename) {
  Clear();

  // bool read_OK = acamcad::polymesh::loadMesh(filename, polyMesh);
  bool read_OK = acamcad::polymesh::loadMesh(filename, xyzMesh);

  polyMesh = xyzMesh;
  std::cout << "Load mesh from file " << filename << std::endl;
  if (read_OK) {
    strMeshFileName = QString::fromStdString(filename);
    QFileInfo fi(strMeshFileName);
    strMeshPath = fi.path();
    strMeshBaseName = fi.baseName();
    UpdateMesh();
    update();
    return true;
  }
  return false;
}

void MeshViewerWidget::CalcUVMesh() {
  using namespace acamcad;
  using acamcad::polymesh::MHalfedge;
  using acamcad::polymesh::MVert;

  std::cout << "CalcUVMesh" << std::endl;
  // 清理uvMesh，每次load
  uvMesh->clear();
  // 计算边界
  // 保存边界节点的index
  std::vector<size_t> boundary;
  for (const auto he : xyzMesh->halfEdges()) {
    if (xyzMesh->isBoundary(he) &&
        std::find(boundary.begin(), boundary.end(), he->toVertex()->index()) ==
            boundary.end()) {
      {
        MHalfedge *cur = he;
        do {
          boundary.push_back(cur->toVertex()->index());
          cur = cur->next();
        } while (cur != he);
      }
    }
  }

  std::cout << "boundary = " << boundary.size() << std::endl;

  if (boundary.size() == 0) {
    std::cout << "ERROR: mesh is no Homeomorphic with disk." << std::endl;
    return;
  }

  // 二层map用来保存节点的权重W_ij
  std::unordered_map<int, std::unordered_map<int, double>> vertex_weights;

  std::cout << "Computing Weight..." << std::endl;
  // 计算权重
  for (const auto vi : xyzMesh->vertices()) {
    // 边界节点不需要计算权重
    if (xyzMesh->isBoundary(vi)) {
      continue;
    }
    // 
    std::vector<MVert *> vjs = xyzMesh->vertAdjacentVertices(vi);
    int n = vjs.size();

    double sum_angles = 0;
    std::vector<double> lengths;
    std::vector<double> angles;

    vertex_weights.emplace(vi->index(), std::unordered_map<int, double>());
    std::unordered_map<int, double> &weights = vertex_weights.at(vi->index());

    for (int j = 0; j < n; j++) {

      // 向量节点i到节点j
      MVector3 ej = vjs[j]->position() - vi->position();
      // 向量节点i到节点j+1
      MVector3 ejp1 = vjs[(j + 1) % n]->position() - vi->position();

      lengths.push_back(norm(ej));
      angles.push_back(vectorAngle(ej, ejp1));
      sum_angles += angles.back();

      // 如果是Floater权重设置初值为0，均匀权重时为1
      if (usingFloaterWeight) {
        weights.emplace(vjs[j]->index(), 0);
      } else {
        weights.emplace(vjs[j]->index(), 1);
      }
    }

    if (!usingFloaterWeight) {
      continue;
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


  // 映射边界到2D， 100(半径)*exp(i theta[])
  double step = 2 * M_PI / boundary.size();
  std::unordered_map<int, std::vector<double>> vert2uv;
  int i = 0;
  for (auto idx : boundary) {
    vert2uv.emplace(idx, std::vector<double>());

    vert2uv[idx].push_back(1 * cos(i * step));
    vert2uv[idx].push_back(1 * sin(i * step));
    i++;
  }

  std::cout << "Contructing matrix A..." << std::endl;
  // 构建矩阵临接权重矩阵A和b
  std::vector<Eigen::Triplet<float>> ATriplets;
  std::vector<double> u_vec, v_vec;
  Eigen::VectorXd ub, vb;

  for (const auto v : xyzMesh->vertices()) {
    // 边界情况
    if (xyzMesh->isBoundary(v)) {
      ATriplets.emplace_back(v->index(), v->index(), 1);
      u_vec.push_back(vert2uv[v->index()][0]);
      v_vec.push_back(vert2uv[v->index()][1]);
    } else {
      double total_weight = 0;
      const auto &weights = vertex_weights.at(v->index());

      for (const auto w : xyzMesh->vertAdjacentVertices(v)) {
        ATriplets.emplace_back(v->index(), w->index(), weights.at(w->index()));
        total_weight += weights.at(w->index());
      }
      ATriplets.emplace_back(v->index(), v->index(), -total_weight);
      u_vec.push_back(0);
      v_vec.push_back(0);
    }
  }

  Eigen::SparseMatrix<double> A(xyzMesh->numVertices(), xyzMesh->numVertices());
  Eigen::Map<Eigen::VectorXd> uu(u_vec.data(), u_vec.size());
  Eigen::Map<Eigen::VectorXd> vv(v_vec.data(), v_vec.size());
  A.setFromTriplets(ATriplets.begin(), ATriplets.end());

  std::cout << "Solving Linear equation..." << std::endl;
  // LU算法求解
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
  ub = solver.solve(uu);
  vb = solver.solve(vv);

  // 映射后的坐标用来构建uvMesh
  std::vector<MVert*> uv_vertices;

  for (auto vert : xyzMesh->vertices()) {
    uv_vertices.push_back( uvMesh->addVertex(ub[vert->index()], vb[vert->index()], 0) );
  }

  for (const auto &fh : xyzMesh->polyfaces()) {
    auto xyz_vertices = xyzMesh->polygonVertices(fh);
    std::vector<MVert *> f_vertices;
    for (MVert *v : xyz_vertices) {
      f_vertices.push_back(uv_vertices[v->index()]);
    }

    uvMesh->addPolyFace(f_vertices);
  }

  DeformationEnergy energy(DeformationType::SD); 
  Energy2DSystem system(xyzMesh, uvMesh, energy);
  system.Solve();
}

void MeshViewerWidget::Clear(void) { polyMesh->clear(); }

void MeshViewerWidget::UpdateMesh(void) {
  polyMesh->updateFacesNormal();
  polyMesh->updateMeshNormal();
  polyMesh->updateVerticesNormal();
  if (polyMesh->numVertices() == 0) {
    std::cerr << "ERROR: UpdateMesh() No vertices!" << std::endl;
    return;
  }
  ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
  ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;

  for (const auto &vh : polyMesh->vertices()) {
    auto p = vh->position();
    for (size_t i = 0; i < 3; i++) {
      ptMin[i] = ptMin[i] < p[i] ? ptMin[i] : p[i];
      ptMax[i] = ptMax[i] > p[i] ? ptMax[i] : p[i];
    }
  }

  double avelen = 0.0;
  double maxlen = 0.0;
  double minlen = DBL_MAX;
  for (const auto &eh : polyMesh->edges()) {
    double len = eh->length();
    maxlen = len > maxlen ? len : maxlen;
    minlen = len < minlen ? len : minlen;
    avelen += len;
  }

  SetScenePosition((ptMin + ptMax) * 0.5, (ptMin - ptMax).norm() * 0.5);
  std::cout << "Information of the input mesh:" << std::endl;
  std::cout << "  [V, E, F] = [" << polyMesh->numVertices() << ", "
            << polyMesh->numEdges() << ", " << polyMesh->numPolygons() << "]\n";
  std::cout << "  BoundingBox:\n";
  std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
  std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
  std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
  std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
  std::cout << "  Edge Length: [" << minlen << ", " << maxlen
            << "]; AVG: " << avelen / polyMesh->numEdges() << std::endl;
}

bool MeshViewerWidget::SaveMesh(const std::string &filename) {
  return acamcad::polymesh::writeMesh(filename, polyMesh);
}

bool MeshViewerWidget::ScreenShot() {
  update();
  QString filename =
      strMeshPath + "/" +
      QDateTime::currentDateTime().toString("yyyyMMddHHmmsszzz") +
      QString(".png");
  QImage image = grabFramebuffer();
  image.save(filename);
  std::cout << "Save screen shot to " << filename.toStdString() << std::endl;
  return true;
}

void MeshViewerWidget::SetDrawBoundingBox(bool b) {
  isDrawBoundingBox = b;
  update();
}
void MeshViewerWidget::SetDrawBoundary(bool b) {
  isDrawBoundary = b;

  UpdateMesh();
  update();
}

void MeshViewerWidget::SetDrawReparam(bool b) {
  usingFloaterWeight = false;
  if (b) {
    CalcUVMesh();
    polyMesh = uvMesh;
  } else {
    polyMesh = xyzMesh;
  }

  UpdateMesh();
  update();
}

void MeshViewerWidget::SetDrawReparam2(bool b) {
  usingFloaterWeight = true;
  if (b) {
    CalcUVMesh();
    polyMesh = uvMesh;
  } else {
    polyMesh = xyzMesh;
  }

  UpdateMesh();
  update();
}

void MeshViewerWidget::EnableLighting(bool b) {
  isEnableLighting = b;
  update();
}
void MeshViewerWidget::EnableDoubleSide(bool b) {
  isTwoSideLighting = b;
  update();
}

void MeshViewerWidget::ResetView(void) {
  ResetModelviewMatrix();
  ViewCenter();
  update();
}

void MeshViewerWidget::ViewCenter(void) {
  if (polyMesh->numVertices() != 0) {
    UpdateMesh();
  }
  update();
}

void MeshViewerWidget::CopyRotation(void) { CopyModelViewMatrix(); }

void MeshViewerWidget::LoadRotation(void) {
  LoadCopyModelViewMatrix();
  update();
}

void MeshViewerWidget::PrintMeshInfo(void) {
  std::cout << "Mesh Info:\n";
  std::cout << "  [V, E, F] = [" << polyMesh->numVertices() << ", "
            << polyMesh->numEdges() << ", " << polyMesh->numPolygons() << "]\n";
  std::cout << "  BoundingBox:\n";
  std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
  std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
  std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
  std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
}

void MeshViewerWidget::DrawScene(void) {
  glMatrixMode(GL_PROJECTION);
  glLoadMatrixd(&projectionmatrix[0]);
  glMatrixMode(GL_MODELVIEW);
  glLoadMatrixd(&modelviewmatrix[0]);
  // DrawAxis();
  if (isDrawBoundingBox)
    DrawBoundingBox();
  if (isDrawBoundary)
    DrawBoundary();
  if (isEnableLighting)
    glEnable(GL_LIGHTING);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);
  DrawSceneMesh();
  if (isEnableLighting)
    glDisable(GL_LIGHTING);
}

void MeshViewerWidget::DrawSceneMesh(void) {
  if (polyMesh->numVertices() == 0) {
    return;
  }
  SetMaterial();
  switch (drawmode) {
  case POINTS:
    DrawPoints();
    break;
  case WIREFRAME:
    DrawWireframe();
    break;
  case HIDDENLINES:
    DrawHiddenLines();
    break;
  case FLATLINES:
    DrawFlatLines();
    break;
  case FLAT:
    glColor3d(0.8, 0.8, 0.8);
    DrawFlat();
    break;
  default:
    break;
  }
}

void MeshViewerWidget::DrawPoints(void) const {
  glColor3d(1.0, 0.5, 0.5);
  glPointSize(5);
  glBegin(GL_POINTS);
  for (const auto &vh : polyMesh->vertices()) {
    // if (isDrawReparam) {
    //   glNormal3d(0, 0, 1);
    // } else {
    glNormal3dv(vh->normal().data());
    // }

    glVertex3dv(vh->position().data());
    // std::cout << "v[" << vh->index() << "] = (" << vh->nx() << ", " <<
    // vh->ny() << ", " << vh->nz() << ")" << std::endl; std::cout << "v[" <<
    // vh->index() << "] = (" << vh->x() << ", " << vh->y() << ", " << vh->z()
    // << ")" << std::endl;
  }
  glEnd();
}

void MeshViewerWidget::DrawWireframe(void) const {
  glColor3d(0.2, 0.2, 0.2);
  glBegin(GL_LINES);
  for (const auto &eh : polyMesh->edges()) {
    auto heh = eh->halfEdge();
    auto v0 = heh->fromVertex();
    auto v1 = heh->toVertex();
    glNormal3dv(v0->normal().data());
    glVertex3dv(v0->position().data());
    glNormal3dv(v1->normal().data());
    glVertex3dv(v1->position().data());
  }
  glEnd();
}

void MeshViewerWidget::DrawHiddenLines() const {
  glLineWidth(1.0);
  float backcolor[4];
  glGetFloatv(GL_COLOR_CLEAR_VALUE, backcolor);
  glColor4fv(backcolor);
  glDepthRange(0.01, 1.0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  if (glIsEnabled(GL_LIGHTING)) {
    glDisable(GL_LIGHTING);
    DrawFlat();
    glEnable(GL_LIGHTING);
  } else {
    DrawFlat();
  }
  glDepthRange(0.0, 1.0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glColor3d(.3, .3, .3);
  DrawFlat();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::DrawFlatLines(void) const {
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.5f, 2.0f);
  glShadeModel(GL_FLAT);
  // glColor3d(0.8, 0.8, 0.8);
  glColor3d(1.0, 1.0, 1.0);
  DrawFlat();
  glDisable(GL_POLYGON_OFFSET_FILL);
  if (glIsEnabled(GL_LIGHTING)) {
    glDisable(GL_LIGHTING);
    DrawWireframe();
    glEnable(GL_LIGHTING);
  } else {
    DrawWireframe();
  }
}

void MeshViewerWidget::DrawFlat(void) const {
  glBegin(GL_TRIANGLES);
  for (const auto &fh : polyMesh->polyfaces()) {
    glNormal3dv(fh->normal().data());
    for (const auto &fvh : polyMesh->polygonVertices(fh)) {
      glVertex3dv(fvh->position().data());
    }
  }
  glEnd();
}

void MeshViewerWidget::DrawBoundingBox(void) const {
  float linewidth;
  glGetFloatv(GL_LINE_WIDTH, &linewidth);
  glLineWidth(2.0f);
  glColor3d(.3, .7, .3);
  glBegin(GL_LINES);
  for (const auto &i : {0, 1}) {
    for (const auto &j : {0, 1}) {
      for (const auto &k : {0, 1}) {
        glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1],
                   k ? ptMin[2] : ptMax[2]);
        glVertex3d(~i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1],
                   k ? ptMin[2] : ptMax[2]);
        glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1],
                   k ? ptMin[2] : ptMax[2]);
        glVertex3d(i ? ptMin[0] : ptMax[0], ~j ? ptMin[1] : ptMax[1],
                   k ? ptMin[2] : ptMax[2]);
        glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1],
                   k ? ptMin[2] : ptMax[2]);
        glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1],
                   ~k ? ptMin[2] : ptMax[2]);
      }
    }
  }
  glEnd();
  glLineWidth(linewidth);
}

void MeshViewerWidget::DrawBoundary(void) const {
  float linewidth;
  glGetFloatv(GL_LINE_WIDTH, &linewidth);
  glLineWidth(10.0f);
  // glColor3d(0.1, 0.1, 0.1);
  glColor3d(0.0, 0.0, 1.0);
  glBegin(GL_LINES);

  for (const auto &eh : polyMesh->edges()) {
    if (polyMesh->isBoundary(eh)) {
      auto heh = eh->halfEdge();
      auto v0 = heh->fromVertex();
      auto v1 = heh->toVertex();
      glNormal3dv(v0->normal().data());
      glVertex3dv(v0->position().data());
      glNormal3dv(v1->normal().data());
      glVertex3dv(v1->position().data());
    }
  }

  glEnd();
  glLineWidth(linewidth);
}

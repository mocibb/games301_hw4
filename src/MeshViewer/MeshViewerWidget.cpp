#include "MeshViewerWidget.h"
#include <QtCore>
#include <cfloat>
#include <unordered_map>
#include <unordered_set>
#include "../AES/AES.h"
#include "../BFF/BFF.h"
#include "../Tutte/Tutte.h"


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
  std::cout << "CalcUVMesh" << std::endl;
  // 清理uvMesh，每次load
  uvMesh->clear();

  TutteEmbeding tutte(xyzMesh, uvMesh);
  tutte.Solve();

  BoundaryFlattenFirst bff(xyzMesh, uvMesh, 1);
  bff.Solve();

  // DeformationEnergy energy(DeformationType::LSCM); 
  // Energy2DSystem system(xyzMesh, uvMesh, energy);
  // auto boundary = tutte.GetBoundary();
  // size_t idx1 = boundary[0];
  // size_t idx2 = boundary[boundary.size()/2];
  // system.AddFix(idx1, Eigen::Vector2d(0, 0));
  // system.AddFix(idx2, Eigen::Vector2d(1, 0));
  // system.Solve();

  // TutteEmbeding tutte2(xyzMesh, uvMesh, MEANVALUE);
  // tutte2.Solve();
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

void MeshViewerWidget::DrawPoints(void) {
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

void MeshViewerWidget::DrawWireframe(void) {
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

void MeshViewerWidget::DrawHiddenLines() {
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

void MeshViewerWidget::DrawFlatLines(void) {
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

void MeshViewerWidget::DrawBoundingBox(void) {
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

void MeshViewerWidget::DrawBoundary(void) {
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

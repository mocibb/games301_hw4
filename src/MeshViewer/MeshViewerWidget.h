#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "../PolyMesh/include/PolyMesh/IOManager.h"

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
public:
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	bool CalcUVMesh(const std::string & filename);
	void Clear(void);
	void UpdateMesh(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	void SetDrawReparam(bool b);
	void SetDrawReparam2(bool b);
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
protected:
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);

private:
	void CalcUVMesh(void);
	void DrawPoints(void);
	void DrawWireframe(void);
	void DrawHiddenLines(void);
	void DrawFlatLines(void);
	void DrawFlat(void) const;
	void DrawBoundingBox(void);
	void DrawBoundary(void);
protected:
	acamcad::polymesh::PolyMesh* polyMesh = new acamcad::polymesh::PolyMesh();
	acamcad::polymesh::PolyMesh* xyzMesh = new acamcad::polymesh::PolyMesh();
	acamcad::polymesh::PolyMesh* uvMesh = new acamcad::polymesh::PolyMesh();
	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	acamcad::MPoint3 ptMin;
	acamcad::MPoint3 ptMax;
	
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;

	// usingFloaterWeight: false: 均匀权重， true: Floater权重
	bool usingFloaterWeight = true;
};

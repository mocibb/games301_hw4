#ifndef AES_UTIL_H_
#define AES_UTIL_H_

#include <Eigen/Dense>
#include <TinyAD/Scalar.hh>

template <typename Scalar, int k>
Eigen::Vector<Scalar, k> warpData(const Eigen::Vector<double, k> &d) {
  return Eigen::Vector<Scalar, k>::Ones();
}

template <typename Scalar, int rows, int cols>
Eigen::Matrix<double, rows, cols>
unwarpData(const Eigen::Matrix<Scalar, rows, cols> &d) {
  return Eigen::Matrix<double, rows, cols>::Ones();
}

template <int k, int rows, int cols>
Eigen::Matrix<double, rows, cols>
unwarpData(const Eigen::Matrix<TinyAD::Double<k>, rows, cols> &d) {
  return TinyAD::to_passive(d);
}

template <typename T>
Eigen::Matrix2<T> createTarget(const Eigen::Vector2d &ui,
                               const Eigen::Vector2d &uj,
                               const Eigen::Vector2d &uk) {

  Eigen::Matrix2<T> target;

#ifdef NUMERIC_CHECK
  Eigen::Vector<double, 6> uc;
  uc << ui, uj, uk;

  Eigen::Vector<T, 6> auc = warpData<T>(uc);
  Eigen::Vector2<T> a(auc[0], auc[1]);
  Eigen::Vector2<T> b(auc[2], auc[3]);
  Eigen::Vector2<T> c(auc[4], auc[5]);
  target << b - a, c - a;
#else
  target << uj - ui, uk - ui;
#endif

  return target;
}

inline Eigen::Vector4d Vectorize(const Eigen::Matrix2d &m) {
    Eigen::Vector4d v(m(0, 0), m(1, 0), m(0, 1), m(1, 1));
    return v;
}

#endif
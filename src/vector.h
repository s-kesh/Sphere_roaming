#ifndef VECTOR_H_
#define VECTOR_H_

typedef struct {
    double x, y, z;
} Vector3D;

typedef struct {
    double a, x, y, z;
} Quat;

double distance(const Vector3D *v1, const Vector3D *v2);
double distance_sq(const Vector3D *v1, const Vector3D *v2);
void cross_product(const Vector3D *v1,
                   const Vector3D *v2,
                   Vector3D *res);
double norm(const Vector3D *vec);
double norm_sq(const Vector3D *vec);
void normalize(const Vector3D *vec, Vector3D *res);
void scalar_multiply(const double scalar,
                     const Vector3D *vec,
                     Vector3D *res);
void add_vectors(const Vector3D *v1, const Vector3D *v2, Vector3D *res);
double dot_product(const Vector3D *v1, const Vector3D *v2);
void rotate_vector(const Vector3D *v, const Vector3D *axis,
                   const double angle, Vector3D* result);
void create_perpend_vector(const Vector3D *vec, Vector3D *res);

void quat_to_pos(const Quat *q, Vector3D *pos);
void quat_random(Quat *result);
void quat_add(const Quat *q1, const Quat *q2, Quat *result);
void quat_scalar_multiply(const Quat *q, const double scalar, Quat *result);
double quat_dot(const Quat *q1, const Quat *q2);
double quat_norm(const Quat *q);
double quat_norm_sq(const Quat *q);
void quat_product(const Quat *q1, const Quat *q2, Quat *result);
void quat_conjugate(const Quat *q, Quat *result);
void quat_inverse(const Quat *q, Quat *result);
void quat_rotate_vector(const Quat *q, const Vector3D *v, Quat *result);

#endif

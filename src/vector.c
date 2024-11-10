#include "vector.h"
#include <math.h>
#include <stdlib.h>


double distance_sq(const Vector3D *v1, const Vector3D *v2) {
    double xx = (v1->x - v2->x);
    double yy = (v1->y - v2->y);
    double zz = (v1->z - v2->z);

    return xx*xx + yy*yy + zz*zz;
}

double distance(const Vector3D *v1, const Vector3D *v2) {
    double xx = (v1->x - v2->x);
    double yy = (v1->y - v2->y);
    double zz = (v1->z - v2->z);

    return sqrt(xx*xx + yy*yy + zz*zz);
}

void cross_product(const Vector3D *v1,
                   const Vector3D *v2,
                   Vector3D *res) {
  res->x = v1->y * v2->z - v1->z * v2->y;
  res->y = v1->z * v2->x - v1->x * v2->z;
  res->z = v1->x * v2->y - v1->y * v2->x;
}

double norm(const Vector3D *v1) {
  return sqrt(v1->x*v1->x + v1->y*v1->y + v1->z*v1->z);
}

double norm_sq(const Vector3D *v1) {
  return v1->x*v1->x + v1->y*v1->y + v1->z*v1->z;
}

void normalize(const Vector3D *vec, Vector3D *res) {
    double nor = norm(vec);
    res->x = vec->x / nor;
    res->y = vec->y / nor;
    res->z = vec->z / nor;
}

void add_vectors(const Vector3D *v1,
                 const Vector3D *v2,
                 Vector3D *res) {
  res->x = v1->x + v2->x;
  res->y = v1->y + v2->y;
  res->z = v1->z + v2->z;
}

void scalar_multiply(const double scalar,
                     const Vector3D *vec,
                     Vector3D *res) {
  res->x = vec->x * scalar;
  res->y = vec->y * scalar;
  res->z = vec->z * scalar;
}

double dot_product(const Vector3D *v1,
                   const Vector3D *v2) {
  return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

void rotate_vector(const Vector3D *v,
                   const Vector3D *axis,
                   const double angle,
                   Vector3D *result) {
    double cosTheta = cos(angle);
    double sinTheta = sqrt(1 - cosTheta*cosTheta);
    double dot = dot_product(axis, v); // Compute dot product of axis and vector

    // Rodrigues' rotation formula
    result->x = cosTheta * v->x + sinTheta * (axis->y * v->z - axis->z * v->y) + (1 - cosTheta) * dot * axis->x;
    result->y = cosTheta * v->y + sinTheta * (axis->z * v->x - axis->x * v->z) + (1 - cosTheta) * dot * axis->y;
    result->z = cosTheta * v->z + sinTheta * (axis->x * v->y - axis->y * v->x) + (1 - cosTheta) * dot * axis->z;
}

void create_perpend_vector(const Vector3D *vec, Vector3D *res) {
  // Create a vector 'a' which is not parallel to vec
  Vector3D a = {0, 0, 0};
  if (vec->x != 0 || vec->y != 0) {
    a.x = -1 * vec->y;
    a.y = vec->x;
    a.z = 0;
  } else {
    a.x = 0;
    a.y = -1 * vec->z;
    a.z = vec->y;
  }
  // Take a cross product
  cross_product(&a, vec, res);
}



void quat_to_pos(const Quat *q, Vector3D *pos) {
    Vector3D vec = {0, 0, 1};
    Quat res = {0, 0, 0, 0};
    quat_rotate_vector(q, &vec, &res);
    pos->x = res.x;
    pos->y = res.y;
    pos->z = res.z;
}

void quat_random(Quat *result) {
    double x, y, z, u, v, w, s;
    do {
        x = 2*((double)rand() / RAND_MAX) - 1;
        y = 2*((double)rand() / RAND_MAX) - 1;
        z = x*x + y*y;
    } while (z > 1);
    do {
        u = 2*((double)rand() / RAND_MAX) - 1;
        v = 2*((double)rand() / RAND_MAX) - 1;
        w = u*u + v*v;
    } while (w > 1);
    s = sqrt((1-z)/w);
    result->a = x;
    result->x = y;
    result->y = s * u;
    result->z = s * v;
}

// Function to add two quaternions
void quat_add(const Quat *q1, const Quat *q2, Quat *result) {
    result->a = q1->a + q2->a;
    result->x = q1->x + q2->x;
    result->y = q1->y + q2->y;
    result->z = q1->z + q2->z;
}

// Function to multiply a quaternion by a scalar
void quat_scalar_multiply(const Quat *q, const double scalar, Quat *result) {
    result->a = q->a * scalar;
    result->x = q->x * scalar;
    result->y = q->y * scalar;
    result->z = q->z * scalar;
}

double quat_norm(const Quat *q) {
    return sqrt(q->a * q->a + q->x * q->x + q->y * q->y + q->z * q->z);
}

double quat_norm_sq(const Quat *q) {
    return q->a * q->a + q->x * q->x + q->y * q->y + q->z * q->z;
}

// Function to calculate the dot product of two quaternions (using vector part only)
double quat_dot(const Quat *q1, const Quat *q2) {
    return q1->x * q2->x + q1->y * q2->y + q1->z * q2->z;
}

// Function to calculate the quaternion product
void quat_product(const Quat *q1, const Quat *q2, Quat *result) {
    result->a = q1->a * q2->a - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z;
    result->x = q1->a * q2->x + q1->x * q2->a + q1->y * q2->z - q1->z * q2->y;
    result->y = q1->a * q2->y - q1->x * q2->z + q1->y * q2->a + q1->z * q2->x;
    result->z = q1->a * q2->z + q1->x * q2->y - q1->y * q2->x + q1->z * q2->a;
}

// Function to calculate the conjugate of a quaternion
void quat_conjugate(const Quat *q, Quat *result) {
    result->a = q->a;
    result->x = -q->x;
    result->y = -q->y;
    result->z = -q->z;
}

// Function to calculate the inverse of a quaternion
void quat_inverse(const Quat *q, Quat *result) {
    Quat conjugate;
    quat_conjugate(q, &conjugate);
    quat_scalar_multiply(&conjugate, 1.0 / quat_norm_sq(q), result);
}

// Function to rotate a vector using a quaternion
void quat_rotate_vector(const Quat *q, const Vector3D *v, Quat *result) {
    // The vector as a quaternion with 0 as the scalar part
    Quat vector_q = {0, v->x, v->y, v->z};

    // Apply the rotation: v' = q * v * q^-1
    Quat q_inv;
    quat_inverse(q, &q_inv);
    Quat temp_result;
    quat_product(q, &vector_q, &temp_result);
    quat_product(&temp_result, &q_inv, result);
}

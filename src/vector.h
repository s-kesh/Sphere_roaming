#ifndef VECTOR_H_
#define VECTOR_H_

/**
 * @brief Represents a 3D vector with x, y, and z components.
 */
typedef struct {
    double x; /**< X-component of the vector */
    double y; /**< Y-component of the vector */
    double z; /**< Z-component of the vector */
} Vector3D;

/**
 * @brief Represents a quaternion with scalar (a) and vector (x, y, z) components.
 */
typedef struct {
    double a; /**< Scalar component of the quaternion */
    double x; /**< X-component of the quaternion */
    double y; /**< Y-component of the quaternion */
    double z; /**< Z-component of the quaternion */
} Quat;

/**
 * @brief Computes the Euclidean distance between two 3D vectors.
 * @param[in] v1 Pointer to the first vector.
 * @param[in] v2 Pointer to the second vector.
 * @return The distance between the two vectors.
 */
double distance(const Vector3D *v1, const Vector3D *v2);

/**
 * @brief Computes the squared Euclidean distance between two 3D vectors.
 * @param[in] v1 Pointer to the first vector.
 * @param[in] v2 Pointer to the second vector.
 * @return The squared distance between the two vectors.
 */
double distance_sq(const Vector3D *v1, const Vector3D *v2);

/**
 * @brief Computes the cross product of two 3D vectors.
 * @param[in] v1 Pointer to the first vector.
 * @param[in] v2 Pointer to the second vector.
 * @param[out] res Pointer to the result vector.
 */
void cross_product(const Vector3D *v1, const Vector3D *v2, Vector3D *res);

/**
 * @brief Computes the norm (magnitude) of a 3D vector.
 * @param[in] vec Pointer to the vector.
 * @return The norm of the vector.
 */
double norm(const Vector3D *vec);

/**
 * @brief Computes the squared norm of a 3D vector.
 * @param[in] vec Pointer to the vector.
 * @return The squared norm of the vector.
 */
double norm_sq(const Vector3D *vec);

/**
 * @brief Normalizes a 3D vector.
 * @param[in] vec Pointer to the vector to normalize.
 * @param[out] res Pointer to the result vector (normalized).
 */
void normalize(const Vector3D *vec, Vector3D *res);

/**
 * @brief Multiplies a 3D vector by a scalar.
 * @param[in] scalar The scalar value.
 * @param[in] vec Pointer to the vector.
 * @param[out] res Pointer to the result vector.
 */
void scalar_multiply(const double scalar, const Vector3D *vec, Vector3D *res);

/**
 * @brief Adds two 3D vectors.
 * @param[in] v1 Pointer to the first vector.
 * @param[in] v2 Pointer to the second vector.
 * @param[out] res Pointer to the result vector.
 */
void add_vectors(const Vector3D *v1, const Vector3D *v2, Vector3D *res);

/**
 * @brief Computes the dot product of two 3D vectors.
 * @param[in] v1 Pointer to the first vector.
 * @param[in] v2 Pointer to the second vector.
 * @return The dot product of the two vectors.
 */
double dot_product(const Vector3D *v1, const Vector3D *v2);

/**
 * @brief Rotates a vector around a specified axis by a given angle.
 * @param[in] v Pointer to the vector to rotate.
 * @param[in] axis Pointer to the axis vector.
 * @param[in] angle Rotation angle in radians.
 * @param[out] result Pointer to the result vector.
 */
void rotate_vector(const Vector3D *v, const Vector3D *axis, const double angle, Vector3D *result);

/**
 * @brief Creates a vector perpendicular to the input vector.
 * @param[in] vec Pointer to the input vector.
 * @param[out] res Pointer to the perpendicular vector.
 */
void create_perpend_vector(const Vector3D *vec, Vector3D *res);

/**
 * @brief Converts a quaternion to a 3D position vector.
 * @param[in] q Pointer to the quaternion.
 * @param[out] pos Pointer to the resulting position vector.
 */
void quat_to_pos(const Quat *q, Vector3D *pos);

/**
 * @brief Generates a random quaternion.
 * @param[out] result Pointer to the resulting random quaternion.
 */
void quat_random(Quat *result);

/**
 * @brief Adds two quaternions.
 * @param[in] q1 Pointer to the first quaternion.
 * @param[in] q2 Pointer to the second quaternion.
 * @param[out] result Pointer to the resulting quaternion.
 */
void quat_add(const Quat *q1, const Quat *q2, Quat *result);

/**
 * @brief Multiplies a quaternion by a scalar.
 * @param[in] q Pointer to the quaternion.
 * @param[in] scalar The scalar value.
 * @param[out] result Pointer to the resulting quaternion.
 */
void quat_scalar_multiply(const Quat *q, const double scalar, Quat *result);

/**
 * @brief Computes the dot product of two quaternions.
 * @param[in] q1 Pointer to the first quaternion.
 * @param[in] q2 Pointer to the second quaternion.
 * @return The dot product of the two quaternions.
 */
double quat_dot(const Quat *q1, const Quat *q2);

/**
 * @brief Computes the norm of a quaternion.
 * @param[in] q Pointer to the quaternion.
 * @return The norm of the quaternion.
 */
double quat_norm(const Quat *q);

/**
 * @brief Computes the squared norm of a quaternion.
 * @param[in] q Pointer to the quaternion.
 * @return The squared norm of the quaternion.
 */
double quat_norm_sq(const Quat *q);

/**
 * @brief Computes the product of two quaternions.
 * @param[in] q1 Pointer to the first quaternion.
 * @param[in] q2 Pointer to the second quaternion.
 * @param[out] result Pointer to the resulting quaternion.
 */
void quat_product(const Quat *q1, const Quat *q2, Quat *result);

/**
 * @brief Computes the conjugate of a quaternion.
 * @param[in] q Pointer to the quaternion.
 * @param[out] result Pointer to the conjugated quaternion.
 */
void quat_conjugate(const Quat *q, Quat *result);

/**
 * @brief Computes the inverse of a quaternion.
 * @param[in] q Pointer to the quaternion.
 * @param[out] result Pointer to the inverted quaternion.
 */
void quat_inverse(const Quat *q, Quat *result);

/**
 * @brief Rotates a vector using a quaternion.
 * @param[in] q Pointer to the quaternion.
 * @param[in] v Pointer to the vector.
 * @param[out] result Pointer to the resulting quaternion after rotation.
 */
void quat_rotate_vector(const Quat *q, const Vector3D *v, Quat *result);

#endif

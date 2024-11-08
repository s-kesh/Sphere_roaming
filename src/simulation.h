#ifndef SIMULATION_H_
#define SIMULATION_H_

struct Vector3D {
    double x, y, z;
};

struct Particle_Pair {
    unsigned int success;
    unsigned int index;
    double time;
    double distance;

    double p1_velocity;
    double p2_velocity;

    struct Vector3D p1_pos;
    struct Vector3D p2_pos;

    struct Vector3D p1_angvel;
    struct Vector3D p2_angvel;

};

double distance_between_pair(struct Particle_Pair *particle);
double distance(struct Vector3D *v1, struct Vector3D *v2);
void cross_product(struct Vector3D *v1,
                   struct Vector3D *v2,
                   struct Vector3D *res);
double norm(struct Vector3D *vec);
void normalize(struct Vector3D *vec, struct Vector3D *res);
void scalar_multiply(double scalar,
                     struct Vector3D *vec,
                     struct Vector3D *res);
void add_vectors(struct Vector3D *v1, struct Vector3D *v2, struct Vector3D *res);
double dot_product(struct Vector3D *v1, struct Vector3D *v2);
void rotate_vector(struct Vector3D* v, struct Vector3D* axis, double angle, struct Vector3D* result);

int initialize_particle_pair(double radius,
                             double velocity,
                             struct Particle_Pair *particle,
                             int no_of_particles);
void simulate_particle(double radius,
                       int list_size,
                       double *rlist, double *force_list,
                       struct Particle_Pair *particle);
double find_force(double distance, int list_size,
                  double *rlist, double *force_list);
#endif

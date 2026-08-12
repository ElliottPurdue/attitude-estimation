/* Dense matrix arithmetic for small, fixed-size problems.
 *
 * Row-major, float, caller-allocated. No dynamic allocation and no size larger
 * than MAT_MAX, so every buffer an estimator needs can live on the stack or in
 * a struct and the whole library remains usable on a microcontroller.
 *
 * Deliberately not a general linear algebra package. It implements exactly the
 * operations a 6-state error-covariance filter performs and nothing else; a
 * fuller library would be more code to test and no more useful here.
 */

#ifndef MATRIX_H
#define MATRIX_H

/* Largest dimension any routine here will see: the 6-element error state. */
#define MAT_MAX 6

#define MAT_AT(a, cols, i, j) ((a)[(i) * (cols) + (j)])

void mat_zero(float *a, int rows, int cols);
void mat_identity(float *a, int n);
void mat_copy(float *dst, const float *src, int rows, int cols);

void mat_add(const float *a, const float *b, float *out, int rows, int cols);
void mat_sub(const float *a, const float *b, float *out, int rows, int cols);
void mat_scale(const float *a, float k, float *out, int rows, int cols);

/* out = a (n x m) * b (m x p), out is n x p. Aliasing out with a or b is not
 * permitted; the caller supplies distinct storage. */
void mat_mul(const float *a, const float *b, float *out, int n, int m, int p);

/* out = a (n x m) * b^T, where b is p x m. out is n x p.
 *
 * Provided because covariance updates are full of A*B^T and materialising the
 * transpose first would cost a copy of the largest matrix in the filter. */
void mat_mul_bt(const float *a, const float *b, float *out, int n, int m, int p);

void mat_transpose(const float *a, float *out, int rows, int cols);

/* Forces exact symmetry by averaging with the transpose.
 *
 * A covariance matrix is symmetric by definition, but the Joseph-form update
 * accumulates float asymmetry over thousands of steps, and an asymmetric P
 * eventually yields a negative innovation variance and a filter that diverges
 * without ever producing a NaN to notice. */
void mat_symmetrize(float *a, int n);

/* Inverts a 3x3 matrix. Returns 0 if the determinant is too small to trust,
 * leaving out untouched, so the caller can skip the update rather than inject
 * garbage into the state. */
int mat3_inverse(const float *a, float *out);

#endif /* MATRIX_H */

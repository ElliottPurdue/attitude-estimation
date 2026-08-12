#include "matrix.h"

#include <math.h>

void mat_zero(float *a, int rows, int cols)
{
    int count = rows * cols;
    for (int i = 0; i < count; ++i) {
        a[i] = 0.0f;
    }
}

void mat_identity(float *a, int n)
{
    mat_zero(a, n, n);
    for (int i = 0; i < n; ++i) {
        MAT_AT(a, n, i, i) = 1.0f;
    }
}

void mat_copy(float *dst, const float *src, int rows, int cols)
{
    int count = rows * cols;
    for (int i = 0; i < count; ++i) {
        dst[i] = src[i];
    }
}

void mat_add(const float *a, const float *b, float *out, int rows, int cols)
{
    int count = rows * cols;
    for (int i = 0; i < count; ++i) {
        out[i] = a[i] + b[i];
    }
}

void mat_sub(const float *a, const float *b, float *out, int rows, int cols)
{
    int count = rows * cols;
    for (int i = 0; i < count; ++i) {
        out[i] = a[i] - b[i];
    }
}

void mat_scale(const float *a, float k, float *out, int rows, int cols)
{
    int count = rows * cols;
    for (int i = 0; i < count; ++i) {
        out[i] = a[i] * k;
    }
}

void mat_mul(const float *a, const float *b, float *out, int n, int m, int p)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < m; ++k) {
                sum += MAT_AT(a, m, i, k) * MAT_AT(b, p, k, j);
            }
            MAT_AT(out, p, i, j) = sum;
        }
    }
}

void mat_mul_bt(const float *a, const float *b, float *out, int n, int m, int p)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < m; ++k) {
                sum += MAT_AT(a, m, i, k) * MAT_AT(b, m, j, k);
            }
            MAT_AT(out, p, i, j) = sum;
        }
    }
}

void mat_transpose(const float *a, float *out, int rows, int cols)
{
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            MAT_AT(out, rows, j, i) = MAT_AT(a, cols, i, j);
        }
    }
}

void mat_symmetrize(float *a, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            float mean = 0.5f * (MAT_AT(a, n, i, j) + MAT_AT(a, n, j, i));
            MAT_AT(a, n, i, j) = mean;
            MAT_AT(a, n, j, i) = mean;
        }
    }
}

int mat3_inverse(const float *a, float *out)
{
    float a00 = MAT_AT(a, 3, 0, 0), a01 = MAT_AT(a, 3, 0, 1), a02 = MAT_AT(a, 3, 0, 2);
    float a10 = MAT_AT(a, 3, 1, 0), a11 = MAT_AT(a, 3, 1, 1), a12 = MAT_AT(a, 3, 1, 2);
    float a20 = MAT_AT(a, 3, 2, 0), a21 = MAT_AT(a, 3, 2, 1), a22 = MAT_AT(a, 3, 2, 2);

    float c00 =  (a11 * a22 - a12 * a21);
    float c01 = -(a10 * a22 - a12 * a20);
    float c02 =  (a10 * a21 - a11 * a20);

    float det = a00 * c00 + a01 * c01 + a02 * c02;

    /* The innovation covariance being near-singular means the measurement
     * carries no usable information along some direction. Reporting failure
     * lets the caller skip the update; inverting anyway would produce an
     * enormous gain and destroy the state. */
    if (fabsf(det) < 1e-20f) {
        return 0;
    }

    float inv_det = 1.0f / det;

    MAT_AT(out, 3, 0, 0) = c00 * inv_det;
    MAT_AT(out, 3, 0, 1) = -(a01 * a22 - a02 * a21) * inv_det;
    MAT_AT(out, 3, 0, 2) =  (a01 * a12 - a02 * a11) * inv_det;

    MAT_AT(out, 3, 1, 0) = c01 * inv_det;
    MAT_AT(out, 3, 1, 1) =  (a00 * a22 - a02 * a20) * inv_det;
    MAT_AT(out, 3, 1, 2) = -(a00 * a12 - a02 * a10) * inv_det;

    MAT_AT(out, 3, 2, 0) = c02 * inv_det;
    MAT_AT(out, 3, 2, 1) = -(a00 * a21 - a01 * a20) * inv_det;
    MAT_AT(out, 3, 2, 2) =  (a00 * a11 - a01 * a10) * inv_det;

    return 1;
}

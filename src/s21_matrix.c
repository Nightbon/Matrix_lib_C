#include "s21_matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPS 1e-7

void s21_get_minor(matrix_t *A, matrix_t *minor, int skip_row, int skip_col);
int s21_incorrect_matrix(const matrix_t *result);
int s21_size_matrix_not_eq(const matrix_t *A, const matrix_t *B);

int s21_incorrect_matrix(const matrix_t *result) {
  return (!result || !result->matrix || result->rows <= 0 ||
          result->columns <= 0);
}

int s21_size_matrix_not_eq(const matrix_t *A, const matrix_t *B) {
  return (A->rows != B->rows || A->columns != B->columns);
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (!result || rows <= 0 || columns <= 0) {
    return INCORRECT_MATRIX;
  }

  int status = OK;

  result->matrix = (double **)malloc(rows * sizeof(double *) +
                                     rows * columns * sizeof(double));

  if (!result->matrix) {
    status = INCORRECT_MATRIX;
  } else {
    memset(result->matrix, 0,
           (rows * sizeof(double *) + rows * columns * sizeof(double)));
    double *data = (double *)(result->matrix + rows);

    for (int i = 0; i < rows; i++) {
      result->matrix[i] = data + i * columns;
    }

    result->rows = rows;
    result->columns = columns;
  }

  return status;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->matrix != NULL) {
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(const matrix_t *A, const matrix_t *B) {
  if (s21_incorrect_matrix(A) || s21_incorrect_matrix(B) ||
      s21_size_matrix_not_eq(A, B)) {
    return FAILURE;
  }

  int status = SUCCESS;

  for (int i = 0; i < A->rows && status == SUCCESS; i++) {
    for (int j = 0; j < A->columns && status == SUCCESS; j++) {
      status = !(fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS);
    }
  }

  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (s21_incorrect_matrix(A) || s21_incorrect_matrix(B) || !result) {
    return INCORRECT_MATRIX;
  } else if (s21_size_matrix_not_eq(A, B)) {
    return CALCULATION_ERROR;
  }

  int status = s21_create_matrix(A->rows, A->columns, result);

  for (int i = 0; i < A->rows && status == OK; i++) {
    for (int j = 0; j < A->columns && status == OK; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      if (isinf(result->matrix[i][j]) || isnan(result->matrix[i][j])) {
        status = CALCULATION_ERROR;
        s21_remove_matrix(result);
      }
    }
  }

  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (s21_incorrect_matrix(A) || s21_incorrect_matrix(B)) {
    return INCORRECT_MATRIX;
  } else if (s21_size_matrix_not_eq(A, B)) {
    return CALCULATION_ERROR;
  }

  int status = s21_create_matrix(A->rows, A->columns, result);

  for (int i = 0; i < A->rows && status == OK; i++) {
    for (int j = 0; j < A->columns && status == OK; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      if (isinf(result->matrix[i][j]) || isnan(result->matrix[i][j])) {
        status = CALCULATION_ERROR;
        s21_remove_matrix(result);
      };
    }
  }

  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (s21_incorrect_matrix(A) || !result) {
    return INCORRECT_MATRIX;
  }

  int status = s21_create_matrix(A->rows, A->columns, result);

  for (int i = 0; i < A->rows && status == OK; i++) {
    for (int j = 0; j < A->columns && status == OK; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
      if (isinf(result->matrix[i][j]) || isnan(result->matrix[i][j])) {
        status = CALCULATION_ERROR;
        s21_remove_matrix(result);
      }
    }
  }

  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (s21_incorrect_matrix(A) || s21_incorrect_matrix(B)) {
    return INCORRECT_MATRIX;
  } else if (A->columns != B->rows) {
    return CALCULATION_ERROR;
  }

  int status = s21_create_matrix(A->rows, B->columns, result);
  int inf_nan = 0;
  for (int i = 0; i < A->rows && status == OK; i++) {
    for (int j = 0; j < B->columns && status == OK; j++) {
      result->matrix[i][j] = 0;
      for (int k = 0; k < A->columns && status == OK; k++) {
        result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        inf_nan = (isinf(result->matrix[i][j]) || isnan(result->matrix[i][j]));
        status = inf_nan * 2;
      }
    }
  }

  if (status != OK) {
    s21_remove_matrix(result);
  }

  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (s21_incorrect_matrix(A)) {
    return INCORRECT_MATRIX;
  }

  int status = s21_create_matrix(A->columns, A->rows, result);

  for (int i = 0; i < A->rows && status == OK; i++) {
    for (int j = 0; j < A->columns && status == OK; j++) {
      result->matrix[j][i] = A->matrix[i][j];
      if (isinf(A->matrix[i][j]) || isnan(A->matrix[i][j])) {
        status = CALCULATION_ERROR;
        s21_remove_matrix(result);
      }
    }
  }

  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (s21_incorrect_matrix(A) || !result) {
    return INCORRECT_MATRIX;
  } else if (A->rows != A->columns) {
    return CALCULATION_ERROR;
  }

  int status = s21_create_matrix(A->rows, A->columns, result);

  if (A->rows == 1 && status == OK) {
    result->matrix[0][0] = 1;
    if (isinf(A->matrix[0][0]) || isnan(A->matrix[0][0])) {
      status = CALCULATION_ERROR;
      s21_remove_matrix(result);
    }
  } else {
    matrix_t minor;
    if (status == OK) {
      status = s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
    }
    int inf_nan = 0;

    for (int i = 0; i < A->rows && A->rows > 1 && status == OK; i++) {
      int sign = (i % 2 == 0) ? 1 : -1;

      for (int j = 0; j < A->columns && status == OK; j++) {
        s21_get_minor(A, &minor, i, j);

        double det_minor = 0;
        status = s21_determinant(&minor, &det_minor);

        result->matrix[i][j] = sign * det_minor;
        sign = -sign;

        inf_nan = (isinf(result->matrix[i][j]) || isnan(result->matrix[i][j]));
        status += (inf_nan * 2);
      }
    }
    s21_remove_matrix(&minor);
  }
  if (status != OK) {
    status = CALCULATION_ERROR;
    s21_remove_matrix(result);
  }

  return status;
}

void s21_get_minor(matrix_t *A, matrix_t *minor, int skip_row, int skip_col) {
  int minor_i = 0;

  for (int i = 0; i < A->rows; i++) {
    int minor_j = 0;

    for (int j = 0; j < A->columns && i != skip_row; j++) {
      if (j != skip_col) {
        minor->matrix[minor_i][minor_j] = A->matrix[i][j];
        minor_j++;
      }
    }
    minor_i += (i != skip_row);
  }
}

int s21_determinant(matrix_t *A, double *result) {
  if (s21_incorrect_matrix(A) || !result) {
    return INCORRECT_MATRIX;
  } else if (A->rows != A->columns) {
    return CALCULATION_ERROR;
  }

  int size = A->rows;
  int status = OK;
  if (size == 1) {
    *result = A->matrix[0][0];
    if (isinf(*result) || isnan(*result)) {
      status = CALCULATION_ERROR;
    }
  } else if (size == 2) {
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];

    if (isinf(*result) || isnan(*result)) {
      status = CALCULATION_ERROR;
    }

  } else {
    double det = 0, sign = 1;
    matrix_t minor;
    status = s21_create_matrix(size - 1, size - 1, &minor);

    for (int col = 0; col < size && status == OK; col++) {
      s21_get_minor(A, &minor, 0, col);
      double det_minor = 0;
      status = s21_determinant(&minor, &det_minor);

      det += sign * A->matrix[0][col] * det_minor;
      sign = -sign;
    }
    s21_remove_matrix(&minor);
    *result = det;
    if (isinf(*result) || isnan(*result)) {
      status = CALCULATION_ERROR;
    }
  }
  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (s21_incorrect_matrix(A) || !result) {
    return INCORRECT_MATRIX;
  } else if (A->rows != A->columns) {
    return CALCULATION_ERROR;
  }

  double det = 0;
  matrix_t complements = {0};
  matrix_t transposed = {0};

  int status = s21_determinant(A, &det);

  if (status == OK && fabs(det) < EPS) {
    status = CALCULATION_ERROR;
  }

  if (status == OK) {
    status = s21_calc_complements(A, &complements);
  }

  if (status == OK) {
    status = s21_transpose(&complements, &transposed);
  }

  if (status == OK) {
    status = s21_mult_number(&transposed, 1.0 / det, result);
  }

  s21_remove_matrix(&complements);
  s21_remove_matrix(&transposed);

  return status;
}
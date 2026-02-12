#include "s21_matrix.h"

#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define INCORRECT_MATRIX 1
#define S21_NAN 0.0 / 0.0L
#define S21_INFINITY 1.0 / 0.0L
#define S21_INFINITY_NEGATIVE -1.0 / 0.0L
#define MAX_DOUBLE 1.7976931348623158e+308
#define MIN_DOUBLE -1.7976931348623158e+308

Suite *suite_create_matrix();
Suite *suite_remove_matrix();
Suite *suite_eq_matrix();
Suite *suite_sum_matrix();
Suite *suite_sub_matrix();
Suite *suite_mult_number();
Suite *suite_mult_matrix();
Suite *suite_transpose();
Suite *suite_calc_complements();
Suite *suite_determinant();
Suite *suite_inverse_matrix();

START_TEST(calc_complements_1) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &result);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;

  result.matrix[0][0] = 0;
  result.matrix[0][1] = 10;
  result.matrix[0][2] = -20;
  result.matrix[1][0] = 4;
  result.matrix[1][1] = -14;
  result.matrix[1][2] = 8;
  result.matrix[2][0] = -8;
  result.matrix[2][1] = -2;
  result.matrix[2][2] = 4;

  res = OK;

  s21_res = s21_calc_complements(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(calc_complements_2) {
  int res, s21_res;
  matrix_t A, s21_result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = S21_NAN;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;

  res = CALCULATION_ERROR;

  s21_res = s21_calc_complements(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(calc_complements_3) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &result);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;

  result.matrix[0][0] = 4;
  result.matrix[0][1] = -3;
  result.matrix[1][0] = -2;
  result.matrix[1][1] = 1;

  res = OK;

  s21_res = s21_calc_complements(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(calc_complements_4) {
  int res, s21_res;
  matrix_t A = {.rows = 0}, s21_result;
  s21_create_matrix(0, 3, &A);

  res = INCORRECT_MATRIX;

  s21_res = s21_calc_complements(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
}
END_TEST

START_TEST(calc_complements_5) {
  int res, s21_res;
  matrix_t A, s21_result;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 1;

  res = OK;

  s21_res = s21_calc_complements(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);

  s21_remove_matrix(&A);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(calc_complements_5_5) {
  int res, s21_res;
  matrix_t A, s21_result;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = S21_NAN;
  ;

  res = CALCULATION_ERROR;

  s21_res = s21_calc_complements(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);

  s21_remove_matrix(&A);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(calc_complements_6) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &result);

  A.matrix[0][0] = 2;
  result.matrix[0][0] = 1;

  res = OK;
  s21_res = s21_calc_complements(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(calc_complements_7) {
  int res, s21_res;
  matrix_t A, s21_result;
  s21_create_matrix(1, 2, &A);
  A.matrix[0][0] = 1;

  res = CALCULATION_ERROR;

  s21_res = s21_calc_complements(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(create_matrix_1) {
  int res, s21_res, rows, columns;
  matrix_t s21_mtrx = {.rows = 0, .columns = 0, .matrix = 0};

  rows = 2;
  columns = 3;

  s21_res = s21_create_matrix(rows, columns, &s21_mtrx);

  res = OK;
  ck_assert_int_eq(rows, s21_mtrx.rows);
  ck_assert_int_eq(columns, s21_mtrx.columns);
  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&s21_mtrx);
}
END_TEST

START_TEST(create_matrix_2) {
  int res, s21_res, rows, columns;
  matrix_t s21_mtrx = {.rows = 0, .columns = 0, .matrix = 0};

  rows = -2;
  columns = 3;

  s21_res = s21_create_matrix(rows, columns, &s21_mtrx);

  res = INCORRECT_MATRIX;
  ck_assert_int_eq(FAILURE, s21_mtrx.rows);
  ck_assert_int_eq(FAILURE, s21_mtrx.columns);
  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&s21_mtrx);
}
END_TEST

START_TEST(create_matrix_3) {
  int res, s21_res, rows, columns;
  matrix_t s21_mtrx = {.rows = 0, .columns = 0, .matrix = 0};

  rows = 2;
  columns = -3;

  s21_res = s21_create_matrix(rows, columns, &s21_mtrx);

  res = INCORRECT_MATRIX;
  ck_assert_int_eq(FAILURE, s21_mtrx.rows);
  ck_assert_int_eq(FAILURE, s21_mtrx.columns);
  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&s21_mtrx);
}
END_TEST

START_TEST(create_matrix_4) {
  int res, s21_res, rows, columns;
  matrix_t s21_mtrx = {.rows = 0, .columns = 0, .matrix = NULL};

  rows = 0;
  columns = 0;

  s21_res = s21_create_matrix(rows, columns, &s21_mtrx);

  res = INCORRECT_MATRIX;
  ck_assert_int_eq(FAILURE, s21_mtrx.rows);
  ck_assert_int_eq(FAILURE, s21_mtrx.columns);
  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&s21_mtrx);
}
END_TEST

START_TEST(create_matrix_5) {
  int res, s21_res, rows, columns;
  matrix_t s21_mtrx;

  rows = 100000;
  columns = 100000;

  s21_res = s21_create_matrix(rows, columns, &s21_mtrx);

  res = INCORRECT_MATRIX;
  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&s21_mtrx);
}
END_TEST

START_TEST(determinant_1) {
  int res, s21_res;
  double det, s21_det;
  matrix_t A;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 4;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 1;

  det = 0;

  res = OK;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  ck_assert_double_eq(det, s21_det);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_2) {
  int res, s21_res;
  double det, s21_det;
  matrix_t A;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = 4;

  det = 4;

  res = OK;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  ck_assert_double_eq(det, s21_det);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_3) {
  int res, s21_res;
  double det, s21_det;
  matrix_t A;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  A.matrix[2][0] = 7;
  A.matrix[2][1] = 8;
  A.matrix[2][2] = 9;

  det = 0;

  res = OK;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  ck_assert_double_eq(det, s21_det);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_4) {
  int res, s21_res;
  double s21_det;
  matrix_t A;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = S21_INFINITY;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  A.matrix[2][0] = 7;
  A.matrix[2][1] = 8;
  A.matrix[2][2] = 9;

  res = CALCULATION_ERROR;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_5) {
  int res, s21_res;
  double s21_det;
  matrix_t A;
  s21_create_matrix(3, 2, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 7;
  A.matrix[2][1] = 8;

  res = CALCULATION_ERROR;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_6) {
  int res, s21_res;
  double s21_det;
  matrix_t A = {.rows = 0, .columns = 2};
  s21_create_matrix(0, 2, &A);

  res = INCORRECT_MATRIX;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
}
END_TEST

START_TEST(determinant_7) {
  int res, s21_res;
  double s21_det;
  matrix_t A;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = MAX_DOUBLE;
  A.matrix[0][1] = S21_INFINITY;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = MAX_DOUBLE;

  res = CALCULATION_ERROR;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_8) {
  int res, s21_res;
  double s21_det;
  matrix_t A;
  s21_create_matrix(1, 1, &A);

  A.matrix[0][0] = S21_NAN;

  res = CALCULATION_ERROR;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_9) {
  int res, s21_res;
  double s21_det;
  matrix_t A;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  A.matrix[2][0] = 7;
  A.matrix[2][1] = 8;
  A.matrix[2][2] = S21_NAN;

  res = CALCULATION_ERROR;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_10) {
  int res, s21_res;
  double det, s21_det;
  matrix_t A;
  s21_create_matrix(5, 5, &A);

  A.matrix[0][0] = 0;
  A.matrix[0][1] = 12;
  A.matrix[0][2] = -12;
  A.matrix[0][3] = 12;
  A.matrix[0][4] = 6;
  A.matrix[1][0] = -3;
  A.matrix[1][1] = -9;
  A.matrix[1][2] = 9;
  A.matrix[1][3] = 9;
  A.matrix[1][4] = -6;
  A.matrix[2][0] = 0;
  A.matrix[2][1] = 0;
  A.matrix[2][2] = -2;
  A.matrix[2][3] = 4;
  A.matrix[2][4] = -2;
  A.matrix[3][0] = -3;
  A.matrix[3][1] = -17;
  A.matrix[3][2] = 13;
  A.matrix[3][3] = 3;
  A.matrix[3][4] = -8;
  A.matrix[4][0] = 0;
  A.matrix[4][1] = 0;
  A.matrix[4][2] = 4;
  A.matrix[4][3] = -8;
  A.matrix[4][4] = 0;

  det = -1728;

  res = OK;
  s21_res = s21_determinant(&A, &s21_det);

  ck_assert_int_eq(res, s21_res);
  ck_assert_double_eq(det, s21_det);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_11_by_py) {
  int res, s21_res;
  double det = 0, s21_det;
  matrix_t A;
  s21_create_matrix(6, 6, &A);
  A.matrix[0][0] = -7;
  A.matrix[0][1] = -8;
  A.matrix[0][2] = -6;
  A.matrix[0][3] = 7;
  A.matrix[0][4] = -5;
  A.matrix[0][5] = -8;
  A.matrix[1][0] = -2;
  A.matrix[1][1] = 8;
  A.matrix[1][2] = -2;
  A.matrix[1][3] = -4;
  A.matrix[1][4] = -1;
  A.matrix[1][5] = -6;
  A.matrix[2][0] = 8;
  A.matrix[2][1] = 0;
  A.matrix[2][2] = 4;
  A.matrix[2][3] = 8;
  A.matrix[2][4] = 9;
  A.matrix[2][5] = 0;
  A.matrix[3][0] = 2;
  A.matrix[3][1] = -7;
  A.matrix[3][2] = 0;
  A.matrix[3][3] = -2;
  A.matrix[3][4] = -1;
  A.matrix[3][5] = 1;
  A.matrix[4][0] = -8;
  A.matrix[4][1] = -1;
  A.matrix[4][2] = -5;
  A.matrix[4][3] = 9;
  A.matrix[4][4] = -1;
  A.matrix[4][5] = -10;
  A.matrix[5][0] = 6;
  A.matrix[5][1] = -4;
  A.matrix[5][2] = 8;
  A.matrix[5][3] = -7;
  A.matrix[5][4] = 8;
  A.matrix[5][5] = 4;

  res = OK;
  s21_res = s21_determinant(&A, &s21_det);
  det = 49020.0;

  ck_assert_int_eq(res, s21_res);
  ck_assert_double_eq_tol(det, s21_det, 1e-7);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_12_by_py) {
  int res, s21_res;
  double det = 0, s21_det;
  matrix_t A;
  s21_create_matrix(7, 7, &A);
  A.matrix[0][0] = -4.9985940449323385;
  A.matrix[0][1] = 5.700687399942708;
  A.matrix[0][2] = 7.443694005987482;
  A.matrix[0][3] = 4.321760989799609;
  A.matrix[0][4] = -0.7414567976537618;
  A.matrix[0][5] = 8.480508699017202;
  A.matrix[0][6] = 2.424682601240471;
  A.matrix[1][0] = -9.163212891555027;
  A.matrix[1][1] = -1.526454501014352;
  A.matrix[1][2] = 0.6344953597552637;
  A.matrix[1][3] = 1.413120548934803;
  A.matrix[1][4] = -8.532181574430304;
  A.matrix[1][5] = 5.45830281562182;
  A.matrix[1][6] = -7.617065682373765;
  A.matrix[2][0] = -2.9569584388342713;
  A.matrix[2][1] = 2.2691348015726165;
  A.matrix[2][2] = 7.670409007436496;
  A.matrix[2][3] = -5.998100089328121;
  A.matrix[2][4] = 4.04094899369222;
  A.matrix[2][5] = -5.582134113145195;
  A.matrix[2][6] = 8.46038488714884;
  A.matrix[3][0] = -1.2505946301994078;
  A.matrix[3][1] = 9.34109771789961;
  A.matrix[3][2] = -7.264788325030737;
  A.matrix[3][3] = -2.47209357135325;
  A.matrix[3][4] = 0.13549205071232018;
  A.matrix[3][5] = -8.356432153751658;
  A.matrix[3][6] = -3.6115751447148665;
  A.matrix[4][0] = -3.0892656613932;
  A.matrix[4][1] = 2.50525092322147;
  A.matrix[4][2] = 1.9808594396727206;
  A.matrix[4][3] = 8.951969411854533;
  A.matrix[4][4] = 5.045085593439952;
  A.matrix[4][5] = -3.4360518262291464;
  A.matrix[4][6] = 6.5067803357850424;
  A.matrix[5][0] = -4.3685562631196575;
  A.matrix[5][1] = 5.036435920934573;
  A.matrix[5][2] = -1.0648117742084677;
  A.matrix[5][3] = -1.5163681119921175;
  A.matrix[5][4] = -9.323970791908682;
  A.matrix[5][5] = -6.090641059768329;
  A.matrix[5][6] = 4.4604803665288895;
  A.matrix[6][0] = -2.985078413456541;
  A.matrix[6][1] = -4.74491926609907;
  A.matrix[6][2] = -2.1923356957386493;
  A.matrix[6][3] = 1.0833504024609666;
  A.matrix[6][4] = -7.769154207079822;
  A.matrix[6][5] = -2.84241391649513;
  A.matrix[6][6] = -1.5800423090636087;

  res = OK;
  s21_res = s21_determinant(&A, &s21_det);
  det = -2877209.7459725398;

  ck_assert_int_eq(res, s21_res);
  ck_assert_double_eq_tol(det, s21_det, 1e-7);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant_13_by_py) {
  int res, s21_res;
  double det = 0, s21_det;
  matrix_t A;
  s21_create_matrix(9, 9, &A);
  A.matrix[0][0] = -5;
  A.matrix[0][1] = 0;
  A.matrix[0][2] = 2;
  A.matrix[0][3] = -3;
  A.matrix[0][4] = 2;
  A.matrix[0][5] = 1;
  A.matrix[0][6] = 2;
  A.matrix[0][7] = -3;
  A.matrix[0][8] = -5;
  A.matrix[1][0] = -3;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = -5;
  A.matrix[1][3] = 2;
  A.matrix[1][4] = -5;
  A.matrix[1][5] = 4;
  A.matrix[1][6] = -2;
  A.matrix[1][7] = -2;
  A.matrix[1][8] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = -4;
  A.matrix[2][3] = 3;
  A.matrix[2][4] = 3;
  A.matrix[2][5] = -5;
  A.matrix[2][6] = 0;
  A.matrix[2][7] = 0;
  A.matrix[2][8] = -4;
  A.matrix[3][0] = 1;
  A.matrix[3][1] = -4;
  A.matrix[3][2] = 1;
  A.matrix[3][3] = 0;
  A.matrix[3][4] = 3;
  A.matrix[3][5] = 4;
  A.matrix[3][6] = 1;
  A.matrix[3][7] = -5;
  A.matrix[3][8] = 0;
  A.matrix[4][0] = 1;
  A.matrix[4][1] = -3;
  A.matrix[4][2] = 0;
  A.matrix[4][3] = -2;
  A.matrix[4][4] = -4;
  A.matrix[4][5] = 2;
  A.matrix[4][6] = 3;
  A.matrix[4][7] = 2;
  A.matrix[4][8] = -3;
  A.matrix[5][0] = -2;
  A.matrix[5][1] = 1;
  A.matrix[5][2] = 1;
  A.matrix[5][3] = 4;
  A.matrix[5][4] = 4;
  A.matrix[5][5] = 2;
  A.matrix[5][6] = -3;
  A.matrix[5][7] = 3;
  A.matrix[5][8] = 4;
  A.matrix[6][0] = 3;
  A.matrix[6][1] = -4;
  A.matrix[6][2] = -1;
  A.matrix[6][3] = 0;
  A.matrix[6][4] = 1;
  A.matrix[6][5] = -4;
  A.matrix[6][6] = -3;
  A.matrix[6][7] = 0;
  A.matrix[6][8] = -2;
  A.matrix[7][0] = -3;
  A.matrix[7][1] = 2;
  A.matrix[7][2] = -2;
  A.matrix[7][3] = 1;
  A.matrix[7][4] = 0;
  A.matrix[7][5] = 0;
  A.matrix[7][6] = -1;
  A.matrix[7][7] = 2;
  A.matrix[7][8] = 4;
  A.matrix[8][0] = -1;
  A.matrix[8][1] = 2;
  A.matrix[8][2] = -3;
  A.matrix[8][3] = -3;
  A.matrix[8][4] = 2;
  A.matrix[8][5] = 0;
  A.matrix[8][6] = -1;
  A.matrix[8][7] = -3;
  A.matrix[8][8] = -3;

  res = OK;
  s21_res = s21_determinant(&A, &s21_det);
  det = 12013892.999999976;

  ck_assert_int_eq(res, s21_res);
  ck_assert_double_eq_tol(det, s21_det, 1e-7);

  s21_remove_matrix(&A);
}
END_TEST

START_TEST(eq_matrix_1) {
  int res, s21_res;
  matrix_t A, B;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[0][2] = 3;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;
  B.matrix[1][2] = 6;

  res = SUCCESS;
  s21_res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(eq_matrix_2) {
  int res, s21_res;
  matrix_t A, B;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[0][2] = 3;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 8;
  B.matrix[1][2] = 6;

  res = FAILURE;
  s21_res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(eq_matrix_3) {
  int res, s21_res;
  matrix_t A, B;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = FAILURE;
  s21_res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(eq_matrix_4) {
  int res, s21_res;
  matrix_t A, B;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[0][2] = 3;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;
  B.matrix[1][2] = 6;

  res = FAILURE;
  s21_res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(eq_matrix_5) {
  int res, s21_res;
  matrix_t A, B;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[0][2] = 3;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;
  B.matrix[1][2] = 6;

  res = FAILURE;
  s21_res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(eq_matrix_6) {
  int res, s21_res;
  matrix_t A, B;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(3, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 4;
  B.matrix[2][0] = 5;
  B.matrix[2][1] = 6;

  res = FAILURE;
  s21_res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(eq_matrix_7) {
  int res, s21_res;
  matrix_t A, B;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = S21_INFINITY;
  A.matrix[0][1] = S21_INFINITY_NEGATIVE;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  B.matrix[0][0] = S21_INFINITY;
  B.matrix[0][1] = S21_INFINITY_NEGATIVE;
  B.matrix[0][2] = 3;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;
  B.matrix[1][2] = 6;

  res = SUCCESS;
  s21_res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(eq_matrix_8) {
  int res, s21_res;
  matrix_t A, B;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 0.0 / 0.0L;
  A.matrix[1][2] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[0][2] = 3;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 0.0 / 0.0L;
  B.matrix[1][2] = 6;

  res = SUCCESS;
  s21_res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(zero_matrix) {
  matrix_t A = {0};
  matrix_t B = {0};
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
}

START_TEST(zero_matrix_1) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(0, 0, &A);
  s21_create_matrix(0, 0, &B);
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_1) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 1;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_2) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 2;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_3) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_4) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_5) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.01;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1.01;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_6) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.01;
  A.matrix[0][1] = -2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = -4;
  B.matrix[0][0] = 1.01;
  B.matrix[0][1] = -2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = -4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_7) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.00000000234;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_8) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.0001;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(inverse_matrix_1) {
  int res, s21_res;
  matrix_t A, result, s21_result = {.rows = 0, .columns = 0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &result);

  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  result.matrix[0][0] = 1;
  result.matrix[0][1] = -1;
  result.matrix[0][2] = 1;
  result.matrix[1][0] = -38;
  result.matrix[1][1] = 41;
  result.matrix[1][2] = -34;
  result.matrix[2][0] = 27;
  result.matrix[2][1] = -29;
  result.matrix[2][2] = 24;

  res = OK;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(inverse_matrix_2) {
  int res, s21_res;
  matrix_t A, s21_result = {.rows = 0, .columns = 0};
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = S21_NAN;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;

  res = CALCULATION_ERROR;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(inverse_matrix_3) {
  int res, s21_res;
  matrix_t A, result;
  matrix_t s21_result;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &result);

  A.matrix[0][0] = 5;
  result.matrix[0][0] = 0.2;

  res = OK;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(inverse_matrix_4) {
  int res, s21_res;
  matrix_t A;
  matrix_t s21_result;
  s21_create_matrix(1, 2, &A);

  A.matrix[0][0] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(inverse_matrix_5) {
  int res, s21_res;
  matrix_t A;
  matrix_t s21_result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = S21_INFINITY;

  res = CALCULATION_ERROR;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST
START_TEST(s21_inverse_1) {
  matrix_t A, C;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &C);
  C.matrix[0][0] = 1.0;
  C.matrix[0][1] = -1.0;
  C.matrix[0][2] = 1.0;
  C.matrix[1][0] = -38.0;
  C.matrix[1][1] = 41.0;
  C.matrix[1][2] = -34.0;
  C.matrix[2][0] = 27.0;
  C.matrix[2][1] = -29.0;
  C.matrix[2][2] = 24.0;
  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 5.0;
  A.matrix[0][2] = 7.0;
  A.matrix[1][0] = 6.0;
  A.matrix[1][1] = 3.0;
  A.matrix[1][2] = 4.0;
  A.matrix[2][0] = 5.0;
  A.matrix[2][1] = -2.0;
  A.matrix[2][2] = -3.0;
  matrix_t B;
  s21_inverse_matrix(&A, &B);
  int res = s21_eq_matrix(&B, &C);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(test_incorrect) {
  matrix_t m = {0};
  matrix_t result = {0};
  int code = s21_inverse_matrix(&m, &result);
  ck_assert_int_eq(code, INCORRECT_MATRIX);
}
END_TEST

START_TEST(determinant) {
  const int size = 2;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 1;
  m.matrix[0][1] = 1;
  m.matrix[1][0] = 1;
  m.matrix[1][1] = 1;

  matrix_t result = {0};
  int code = s21_inverse_matrix(&m, &result);
  ck_assert_int_eq(code, CALCULATION_ERROR);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(inverse) {
  const int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  m.matrix[0][0] = 2;
  m.matrix[0][1] = 5;
  m.matrix[0][2] = 7;
  m.matrix[1][0] = 6;
  m.matrix[1][1] = 3;
  m.matrix[1][2] = 4;
  m.matrix[2][0] = 5;
  m.matrix[2][1] = -2;
  m.matrix[2][2] = -3;

  matrix_t res = {0};
  s21_inverse_matrix(&m, &res);

  matrix_t expected = {0};
  s21_create_matrix(size, size, &expected);
  expected.matrix[0][0] = 1;
  expected.matrix[0][1] = -1;
  expected.matrix[0][2] = 1;
  expected.matrix[1][0] = -38;
  expected.matrix[1][1] = 41;
  expected.matrix[1][2] = -34;
  expected.matrix[2][0] = 27;
  expected.matrix[2][1] = -29;
  expected.matrix[2][2] = 24;

  ck_assert_int_eq(s21_eq_matrix(&expected, &res), SUCCESS);

  s21_remove_matrix(&expected);
  s21_remove_matrix(&res);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(inverse_matrix_6_by_py) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(7, 7, &A);
  s21_create_matrix(7, 7, &result);
  A.matrix[0][0] = -355;
  A.matrix[0][1] = 471;
  A.matrix[0][2] = 633;
  A.matrix[0][3] = 618;
  A.matrix[0][4] = 341;
  A.matrix[0][5] = -696;
  A.matrix[0][6] = 968;
  A.matrix[1][0] = 644;
  A.matrix[1][1] = -160;
  A.matrix[1][2] = 335;
  A.matrix[1][3] = -476;
  A.matrix[1][4] = -256;
  A.matrix[1][5] = 446;
  A.matrix[1][6] = -917;
  A.matrix[2][0] = 460;
  A.matrix[2][1] = -287;
  A.matrix[2][2] = 868;
  A.matrix[2][3] = 434;
  A.matrix[2][4] = 188;
  A.matrix[2][5] = 103;
  A.matrix[2][6] = -752;
  A.matrix[3][0] = -461;
  A.matrix[3][1] = 196;
  A.matrix[3][2] = -448;
  A.matrix[3][3] = 862;
  A.matrix[3][4] = 971;
  A.matrix[3][5] = 193;
  A.matrix[3][6] = -858;
  A.matrix[4][0] = 395;
  A.matrix[4][1] = -620;
  A.matrix[4][2] = 852;
  A.matrix[4][3] = 931;
  A.matrix[4][4] = -860;
  A.matrix[4][5] = -177;
  A.matrix[4][6] = 867;
  A.matrix[5][0] = -411;
  A.matrix[5][1] = 440;
  A.matrix[5][2] = 661;
  A.matrix[5][3] = 216;
  A.matrix[5][4] = 407;
  A.matrix[5][5] = -274;
  A.matrix[5][6] = 701;
  A.matrix[6][0] = 432;
  A.matrix[6][1] = 129;
  A.matrix[6][2] = 742;
  A.matrix[6][3] = 399;
  A.matrix[6][4] = 96;
  A.matrix[6][5] = 755;
  A.matrix[6][6] = -788;
  result.matrix[0][0] = 0.0015534522165479704;
  result.matrix[0][1] = -0.01337361185688095;
  result.matrix[0][2] = 0.0064333704407399635;
  result.matrix[0][3] = -0.007005360873922197;
  result.matrix[0][4] = -0.005094002217970681;
  result.matrix[0][5] = -0.00817126403853511;
  result.matrix[0][6] = 0.006085650841082952;
  result.matrix[1][0] = 0.0020911500612235453;
  result.matrix[1][1] = 0.004310652787912465;
  result.matrix[1][2] = -0.0029402868249920363;
  result.matrix[1][3] = 0.0016340453178715247;
  result.matrix[1][4] = 0.0007938148142016136;
  result.matrix[1][5] = 0.000186407739329903;
  result.matrix[1][6] = -0.0003815244522944591;
  result.matrix[2][0] = -0.0005062111273943297;
  result.matrix[2][1] = 0.003903608010625649;
  result.matrix[2][2] = -0.0012833024451957397;
  result.matrix[2][3] = 0.0016769173449664655;
  result.matrix[2][4] = 0.0014213454460315713;
  result.matrix[2][5] = 0.0028648426972597307;
  result.matrix[2][6] = -0.001653313867724643;
  result.matrix[3][0] = 0.0007168461800648925;
  result.matrix[3][1] = 0.0005251697854402448;
  result.matrix[3][2] = -0.0006111174992759948;
  result.matrix[3][3] = 0.00071306187557715;
  result.matrix[3][4] = 0.0005491851342214948;
  result.matrix[3][5] = -0.0005697013870154505;
  result.matrix[3][6] = 0.00017368339890440346;
  result.matrix[4][0] = -0.0005899378987715673;
  result.matrix[4][1] = -0.012761141801191147;
  result.matrix[4][2] = 0.0066188212018627145;
  result.matrix[4][3] = -0.005890695749071904;
  result.matrix[4][4] = -0.004628896463455082;
  result.matrix[4][5] = -0.004970535384558988;
  result.matrix[4][6] = 0.004708340122072248;
  result.matrix[5][0] = -0.0014300814823182252;
  result.matrix[5][1] = -0.004832471891314581;
  result.matrix[5][2] = 0.0014787963761522544;
  result.matrix[5][3] = -0.0020505473668024753;
  result.matrix[5][4] = -0.0012281189215908297;
  result.matrix[5][5] = -0.0008385971471130454;
  result.matrix[5][6] = 0.0025910366880151827;
  result.matrix[6][0] = -0.0003617802952620184;
  result.matrix[6][1] = -0.008869152430347502;
  result.matrix[6][2] = 0.0037509782142171365;
  result.matrix[6][3] = -0.004315240795755499;
  result.matrix[6][4] = -0.00278686231461204;
  result.matrix[6][5] = -0.0034490473113122674;
  result.matrix[6][6] = 0.0035920795006134784;

  res = OK;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq_tol(result.matrix[0][0], s21_result.matrix[0][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][1], s21_result.matrix[0][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][2], s21_result.matrix[0][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][3], s21_result.matrix[0][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][4], s21_result.matrix[0][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][5], s21_result.matrix[0][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][6], s21_result.matrix[0][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][0], s21_result.matrix[1][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][1], s21_result.matrix[1][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][2], s21_result.matrix[1][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][3], s21_result.matrix[1][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][4], s21_result.matrix[1][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][5], s21_result.matrix[1][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][6], s21_result.matrix[1][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][0], s21_result.matrix[2][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][1], s21_result.matrix[2][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][2], s21_result.matrix[2][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][3], s21_result.matrix[2][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][4], s21_result.matrix[2][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][5], s21_result.matrix[2][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][6], s21_result.matrix[2][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][0], s21_result.matrix[3][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][1], s21_result.matrix[3][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][2], s21_result.matrix[3][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][3], s21_result.matrix[3][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][4], s21_result.matrix[3][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][5], s21_result.matrix[3][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][6], s21_result.matrix[3][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][0], s21_result.matrix[4][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][1], s21_result.matrix[4][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][2], s21_result.matrix[4][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][3], s21_result.matrix[4][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][4], s21_result.matrix[4][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][5], s21_result.matrix[4][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][6], s21_result.matrix[4][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][0], s21_result.matrix[5][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][1], s21_result.matrix[5][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][2], s21_result.matrix[5][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][3], s21_result.matrix[5][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][4], s21_result.matrix[5][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][5], s21_result.matrix[5][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][6], s21_result.matrix[5][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][0], s21_result.matrix[6][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][1], s21_result.matrix[6][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][2], s21_result.matrix[6][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][3], s21_result.matrix[6][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][4], s21_result.matrix[6][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][5], s21_result.matrix[6][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][6], s21_result.matrix[6][6], 1e-7);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(inverse_matrix_7_by_py) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(9, 9, &A);
  s21_create_matrix(9, 9, &result);
  A.matrix[0][0] = -371;
  A.matrix[0][1] = -673;
  A.matrix[0][2] = -498;
  A.matrix[0][3] = 630;
  A.matrix[0][4] = -12;
  A.matrix[0][5] = 333;
  A.matrix[0][6] = -664;
  A.matrix[0][7] = -888;
  A.matrix[0][8] = 16;
  A.matrix[1][0] = -146;
  A.matrix[1][1] = -731;
  A.matrix[1][2] = -260;
  A.matrix[1][3] = 173;
  A.matrix[1][4] = -19;
  A.matrix[1][5] = -408;
  A.matrix[1][6] = 648;
  A.matrix[1][7] = 138;
  A.matrix[1][8] = -586;
  A.matrix[2][0] = 74;
  A.matrix[2][1] = 707;
  A.matrix[2][2] = 564;
  A.matrix[2][3] = 979;
  A.matrix[2][4] = -149;
  A.matrix[2][5] = -37;
  A.matrix[2][6] = -455;
  A.matrix[2][7] = -955;
  A.matrix[2][8] = 966;
  A.matrix[3][0] = -717;
  A.matrix[3][1] = -509;
  A.matrix[3][2] = -95;
  A.matrix[3][3] = 350;
  A.matrix[3][4] = -647;
  A.matrix[3][5] = 802;
  A.matrix[3][6] = -898;
  A.matrix[3][7] = 84;
  A.matrix[3][8] = -191;
  A.matrix[4][0] = -57;
  A.matrix[4][1] = 798;
  A.matrix[4][2] = 23;
  A.matrix[4][3] = 240;
  A.matrix[4][4] = 949;
  A.matrix[4][5] = 103;
  A.matrix[4][6] = -558;
  A.matrix[4][7] = 366;
  A.matrix[4][8] = -121;
  A.matrix[5][0] = -963;
  A.matrix[5][1] = -158;
  A.matrix[5][2] = 89;
  A.matrix[5][3] = -353;
  A.matrix[5][4] = 487;
  A.matrix[5][5] = 298;
  A.matrix[5][6] = 770;
  A.matrix[5][7] = 149;
  A.matrix[5][8] = -947;
  A.matrix[6][0] = -57;
  A.matrix[6][1] = 694;
  A.matrix[6][2] = -483;
  A.matrix[6][3] = 430;
  A.matrix[6][4] = -139;
  A.matrix[6][5] = 795;
  A.matrix[6][6] = 323;
  A.matrix[6][7] = 626;
  A.matrix[6][8] = -756;
  A.matrix[7][0] = 54;
  A.matrix[7][1] = -106;
  A.matrix[7][2] = 681;
  A.matrix[7][3] = -424;
  A.matrix[7][4] = 697;
  A.matrix[7][5] = 635;
  A.matrix[7][6] = 152;
  A.matrix[7][7] = -299;
  A.matrix[7][8] = 826;
  A.matrix[8][0] = -378;
  A.matrix[8][1] = -324;
  A.matrix[8][2] = 327;
  A.matrix[8][3] = 945;
  A.matrix[8][4] = 762;
  A.matrix[8][5] = 651;
  A.matrix[8][6] = 808;
  A.matrix[8][7] = -581;
  A.matrix[8][8] = 590;
  result.matrix[0][0] = 0.0004871879707275386;
  result.matrix[0][1] = 0.0027664936987099485;
  result.matrix[0][2] = 0.001044193333692249;
  result.matrix[0][3] = -0.0002649294058620894;
  result.matrix[0][4] = 0.00022858719814165988;
  result.matrix[0][5] = -0.0006713609057233434;
  result.matrix[0][6] = 0.001283326983203877;
  result.matrix[0][7] = 0.0023046091234232;
  result.matrix[0][8] = -0.0016736502378716614;
  result.matrix[1][0] = -0.00010454149204021965;
  result.matrix[1][1] = -0.0012521396036668487;
  result.matrix[1][2] = -2.619244763389076e-05;
  result.matrix[1][3] = -0.00034016642623276434;
  result.matrix[1][4] = -0.00010345064374981778;
  result.matrix[1][5] = 0.000340089048967443;
  result.matrix[1][6] = 5.438530817879496e-06;
  result.matrix[1][7] = -0.0007461063019002741;
  result.matrix[1][8] = 0.0002681205241934487;
  result.matrix[2][0] = -0.0003376425480454883;
  result.matrix[2][1] = 0.003086722220201918;
  result.matrix[2][2] = 0.0017330859465075133;
  result.matrix[2][3] = 0.0005571231241316904;
  result.matrix[2][4] = 0.000411160138382318;
  result.matrix[2][5] = 0.0001055369900197317;
  result.matrix[2][6] = 0.0006852160724666511;
  result.matrix[2][7] = 0.0023086948107837895;
  result.matrix[2][8] = -0.0016827011227792015;
  result.matrix[3][0] = -0.00013369716650763687;
  result.matrix[3][1] = 0.0011048385969511596;
  result.matrix[3][2] = 0.000505521217338403;
  result.matrix[3][3] = 0.00028796709222087094;
  result.matrix[3][4] = 0.0003633529011658943;
  result.matrix[3][5] = -0.00033646849193091894;
  result.matrix[3][6] = 0.00027825894100760125;
  result.matrix[3][7] = 0.00015545912890434558;
  result.matrix[3][8] = 3.9876863139052035e-05;
  result.matrix[4][0] = 0.0002601886222516985;
  result.matrix[4][1] = 0.00029557501884028183;
  result.matrix[4][2] = -0.0001634174239609346;
  result.matrix[4][3] = -0.00023222197906015325;
  result.matrix[4][4] = 0.0006150831373482673;
  result.matrix[4][5] = -1.5390409449917066e-05;
  result.matrix[4][6] = -8.439645253517553e-05;
  result.matrix[4][7] = 0.0003072583240159787;
  result.matrix[4][8] = 4.20373300786936e-05;
  result.matrix[5][0] = 0.0002966598188950942;
  result.matrix[5][1] = 0.00045769703415441326;
  result.matrix[5][2] = 0.0001785699062545127;
  result.matrix[5][3] = 0.0001312481235925557;
  result.matrix[5][4] = -0.00014439513775315254;
  result.matrix[5][5] = -8.796988474187099e-05;
  result.matrix[5][6] = 0.0007987097688771858;
  result.matrix[5][7] = 0.0010407862469558708;
  result.matrix[5][8] = -0.0004078142036552484;
  result.matrix[6][0] = -0.0002321475395181368;
  result.matrix[6][1] = -0.0006116815634127405;
  result.matrix[6][2] = -0.0002882015469154899;
  result.matrix[6][3] = -0.00032324809175015066;
  result.matrix[6][4] = -0.000486668189324292;
  result.matrix[6][5] = 0.00015185729323988597;
  result.matrix[6][6] = -2.201140261974975e-06;
  result.matrix[6][7] = -0.000518891035817604;
  result.matrix[6][8] = 0.0006335478179758123;
  result.matrix[7][0] = -0.0008847217715884261;
  result.matrix[7][1] = -0.0005834491722682291;
  result.matrix[7][2] = -0.0007714429821417575;
  result.matrix[7][3] = 0.0005547318626044673;
  result.matrix[7][4] = 0.00030454005130477386;
  result.matrix[7][5] = -0.00042642140901182925;
  result.matrix[7][6] = -0.00044663758150360205;
  result.matrix[7][7] = -0.0008460647923868803;
  result.matrix[7][8] = 0.0008773590146968456;
  result.matrix[8][0] = -0.0005606760416002568;
  result.matrix[8][1] = -0.003019192293291078;
  result.matrix[8][2] = -0.0014665836990178966;
  result.matrix[8][3] = 1.750659091072698e-05;
  result.matrix[8][4] = -0.00038890907411891077;
  result.matrix[8][5] = -0.00027387980672137626;
  result.matrix[8][6] = -0.0012093704195625097;
  result.matrix[8][7] = -0.0021295434749275016;
  result.matrix[8][8] = 0.0020306515984399955;

  res = OK;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq_tol(result.matrix[0][0], s21_result.matrix[0][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][1], s21_result.matrix[0][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][2], s21_result.matrix[0][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][3], s21_result.matrix[0][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][4], s21_result.matrix[0][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][5], s21_result.matrix[0][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][6], s21_result.matrix[0][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][7], s21_result.matrix[0][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][8], s21_result.matrix[0][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][0], s21_result.matrix[1][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][1], s21_result.matrix[1][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][2], s21_result.matrix[1][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][3], s21_result.matrix[1][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][4], s21_result.matrix[1][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][5], s21_result.matrix[1][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][6], s21_result.matrix[1][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][7], s21_result.matrix[1][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][8], s21_result.matrix[1][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][0], s21_result.matrix[2][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][1], s21_result.matrix[2][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][2], s21_result.matrix[2][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][3], s21_result.matrix[2][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][4], s21_result.matrix[2][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][5], s21_result.matrix[2][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][6], s21_result.matrix[2][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][7], s21_result.matrix[2][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][8], s21_result.matrix[2][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][0], s21_result.matrix[3][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][1], s21_result.matrix[3][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][2], s21_result.matrix[3][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][3], s21_result.matrix[3][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][4], s21_result.matrix[3][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][5], s21_result.matrix[3][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][6], s21_result.matrix[3][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][7], s21_result.matrix[3][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][8], s21_result.matrix[3][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][0], s21_result.matrix[4][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][1], s21_result.matrix[4][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][2], s21_result.matrix[4][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][3], s21_result.matrix[4][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][4], s21_result.matrix[4][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][5], s21_result.matrix[4][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][6], s21_result.matrix[4][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][7], s21_result.matrix[4][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][8], s21_result.matrix[4][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][0], s21_result.matrix[5][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][1], s21_result.matrix[5][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][2], s21_result.matrix[5][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][3], s21_result.matrix[5][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][4], s21_result.matrix[5][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][5], s21_result.matrix[5][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][6], s21_result.matrix[5][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][7], s21_result.matrix[5][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][8], s21_result.matrix[5][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][0], s21_result.matrix[6][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][1], s21_result.matrix[6][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][2], s21_result.matrix[6][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][3], s21_result.matrix[6][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][4], s21_result.matrix[6][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][5], s21_result.matrix[6][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][6], s21_result.matrix[6][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][7], s21_result.matrix[6][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][8], s21_result.matrix[6][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][0], s21_result.matrix[7][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][1], s21_result.matrix[7][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][2], s21_result.matrix[7][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][3], s21_result.matrix[7][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][4], s21_result.matrix[7][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][5], s21_result.matrix[7][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][6], s21_result.matrix[7][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][7], s21_result.matrix[7][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][8], s21_result.matrix[7][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][0], s21_result.matrix[8][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][1], s21_result.matrix[8][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][2], s21_result.matrix[8][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][3], s21_result.matrix[8][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][4], s21_result.matrix[8][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][5], s21_result.matrix[8][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][6], s21_result.matrix[8][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][7], s21_result.matrix[8][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][8], s21_result.matrix[8][8], 1e-7);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}

START_TEST(inverse_matrix_8_by_py) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &result);
  A.matrix[0][0] = -771;
  A.matrix[0][1] = -941;
  A.matrix[0][2] = 224;
  A.matrix[0][3] = -373;
  A.matrix[0][4] = -804;
  A.matrix[1][0] = -55;
  A.matrix[1][1] = 89;
  A.matrix[1][2] = -821;
  A.matrix[1][3] = 697;
  A.matrix[1][4] = -946;
  A.matrix[2][0] = -636;
  A.matrix[2][1] = 219;
  A.matrix[2][2] = 994;
  A.matrix[2][3] = -397;
  A.matrix[2][4] = 186;
  A.matrix[3][0] = -213;
  A.matrix[3][1] = 866;
  A.matrix[3][2] = -555;
  A.matrix[3][3] = 885;
  A.matrix[3][4] = 362;
  A.matrix[4][0] = 863;
  A.matrix[4][1] = -803;
  A.matrix[4][2] = -129;
  A.matrix[4][3] = 134;
  A.matrix[4][4] = 292;
  result.matrix[0][0] = -0.0009708806880338015;
  result.matrix[0][1] = 0.0006947852833967972;
  result.matrix[0][2] = 0.0003690826367722159;
  result.matrix[0][3] = -0.0008510248564653199;
  result.matrix[0][4] = 0.0003976042209418066;
  result.matrix[1][0] = -0.0007900198216592012;
  result.matrix[1][1] = 0.000311222544301572;
  result.matrix[1][2] = 0.00010691384716138141;
  result.matrix[1][3] = -0.0004223993557150103;
  result.matrix[1][4] = -0.0007114274606434718;
  result.matrix[2][0] = -0.0005215837317217009;
  result.matrix[2][1] = 0.000993502986411834;
  result.matrix[2][2] = 0.0019032535679091725;
  result.matrix[2][3] = -0.00028911222153355827;
  result.matrix[2][4] = 0.0009286094705664024;
  result.matrix[3][0] = 2.647206574047746e-05;
  result.matrix[3][1] = 0.0009799961151052592;
  result.matrix[3][2] = 0.0014277809513793278;
  result.matrix[3][3] = 0.0008072448054973678;
  result.matrix[3][4] = 0.0013375684561586201;
  result.matrix[4][0] = 0.00045428958486133127;
  result.matrix[4][1] = -0.0012083753105282598;
  result.matrix[4][2] = -0.0006111966900286977;
  result.matrix[4][3] = 0.0008554160547121283;
  result.matrix[4][4] = 8.954710619288947e-05;

  res = OK;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq_tol(result.matrix[0][0], s21_result.matrix[0][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][1], s21_result.matrix[0][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][2], s21_result.matrix[0][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][3], s21_result.matrix[0][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][4], s21_result.matrix[0][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][0], s21_result.matrix[1][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][1], s21_result.matrix[1][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][2], s21_result.matrix[1][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][3], s21_result.matrix[1][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][4], s21_result.matrix[1][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][0], s21_result.matrix[2][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][1], s21_result.matrix[2][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][2], s21_result.matrix[2][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][3], s21_result.matrix[2][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][4], s21_result.matrix[2][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][0], s21_result.matrix[3][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][1], s21_result.matrix[3][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][2], s21_result.matrix[3][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][3], s21_result.matrix[3][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][4], s21_result.matrix[3][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][0], s21_result.matrix[4][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][1], s21_result.matrix[4][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][2], s21_result.matrix[4][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][3], s21_result.matrix[4][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][4], s21_result.matrix[4][4], 1e-7);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(inverse_matrix_9_by_py) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(6, 6, &A);
  s21_create_matrix(6, 6, &result);
  A.matrix[0][0] = -472;
  A.matrix[0][1] = 7;
  A.matrix[0][2] = -908;
  A.matrix[0][3] = 641;
  A.matrix[0][4] = -23;
  A.matrix[0][5] = 21;
  A.matrix[1][0] = -261;
  A.matrix[1][1] = -631;
  A.matrix[1][2] = 33;
  A.matrix[1][3] = 933;
  A.matrix[1][4] = 386;
  A.matrix[1][5] = -125;
  A.matrix[2][0] = 408;
  A.matrix[2][1] = -604;
  A.matrix[2][2] = 462;
  A.matrix[2][3] = -756;
  A.matrix[2][4] = 776;
  A.matrix[2][5] = 925;
  A.matrix[3][0] = 815;
  A.matrix[3][1] = -720;
  A.matrix[3][2] = -766;
  A.matrix[3][3] = -6;
  A.matrix[3][4] = -662;
  A.matrix[3][5] = 760;
  A.matrix[4][0] = 128;
  A.matrix[4][1] = 318;
  A.matrix[4][2] = 973;
  A.matrix[4][3] = -504;
  A.matrix[4][4] = 765;
  A.matrix[4][5] = -932;
  A.matrix[5][0] = -310;
  A.matrix[5][1] = 60;
  A.matrix[5][2] = 323;
  A.matrix[5][3] = 562;
  A.matrix[5][4] = 41;
  A.matrix[5][5] = 672;
  result.matrix[0][0] = 0.0011401694641858808;
  result.matrix[0][1] = -0.0006096384114763006;
  result.matrix[0][2] = -0.0003679044510756896;
  result.matrix[0][3] = 0.0020274235558828293;
  result.matrix[0][4] = 0.0023952279468206356;
  result.matrix[0][5] = 0.0013864193467685315;
  result.matrix[1][0] = 0.001427714846037204;
  result.matrix[1][1] = -0.0012608642145355553;
  result.matrix[1][2] = -0.0003457172356645252;
  result.matrix[1][3] = 0.0009204151853084899;
  result.matrix[1][4] = 0.0017420608897182168;
  result.matrix[1][5] = 0.0015718506142691753;
  result.matrix[2][0] = -0.0012138733562126784;
  result.matrix[2][1] = 0.0003757750682592097;
  result.matrix[2][2] = -0.0001563727949053855;
  result.matrix[2][3] = -0.0003697448404942149;
  result.matrix[2][4] = -0.00039761482581197535;
  result.matrix[2][5] = 0.0001897873815761778;
  result.matrix[3][0] = 0.0007118005956384697;
  result.matrix[3][1] = 0.00010030007770297647;
  result.matrix[3][2] = -0.0004821264943809499;
  result.matrix[3][3] = 0.0009750600467711139;
  result.matrix[3][4] = 0.0012355170423396912;
  result.matrix[3][5] = 0.001270852607769923;
  result.matrix[4][0] = 0.0015820918877013318;
  result.matrix[4][1] = -0.0002917635161086027;
  result.matrix[4][2] = 0.00046555161028331854;
  result.matrix[4][3] = 0.0006080123637365229;
  result.matrix[4][4] = 0.0015196740062498983;
  result.matrix[4][5] = 0.0006754719772585931;
  result.matrix[5][0] = 0.00029013844501117744;
  result.matrix[5][1] = -0.00041535378160991755;
  result.matrix[5][2] = 0.00031111390012101935;
  result.matrix[5][3] = 0.00017826149030689844;
  result.matrix[5][4] = 1.4522897049731703e-05;
  result.matrix[5][5] = 0.0007920597020218597;

  res = OK;
  s21_res = s21_inverse_matrix(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq_tol(result.matrix[0][0], s21_result.matrix[0][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][1], s21_result.matrix[0][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][2], s21_result.matrix[0][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][3], s21_result.matrix[0][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][4], s21_result.matrix[0][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][5], s21_result.matrix[0][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][0], s21_result.matrix[1][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][1], s21_result.matrix[1][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][2], s21_result.matrix[1][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][3], s21_result.matrix[1][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][4], s21_result.matrix[1][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][5], s21_result.matrix[1][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][0], s21_result.matrix[2][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][1], s21_result.matrix[2][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][2], s21_result.matrix[2][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][3], s21_result.matrix[2][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][4], s21_result.matrix[2][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][5], s21_result.matrix[2][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][0], s21_result.matrix[3][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][1], s21_result.matrix[3][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][2], s21_result.matrix[3][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][3], s21_result.matrix[3][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][4], s21_result.matrix[3][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][5], s21_result.matrix[3][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][0], s21_result.matrix[4][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][1], s21_result.matrix[4][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][2], s21_result.matrix[4][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][3], s21_result.matrix[4][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][4], s21_result.matrix[4][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][5], s21_result.matrix[4][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][0], s21_result.matrix[5][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][1], s21_result.matrix[5][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][2], s21_result.matrix[5][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][3], s21_result.matrix[5][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][4], s21_result.matrix[5][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][5], s21_result.matrix[5][5], 1e-7);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_1) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &B);
  s21_create_matrix(3, 3, &result);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = -1;
  B.matrix[0][2] = 1;
  B.matrix[1][0] = 2;
  B.matrix[1][1] = 3;
  B.matrix[1][2] = 4;
  result.matrix[0][0] = 9;
  result.matrix[0][1] = 11;
  result.matrix[0][2] = 17;
  result.matrix[1][0] = 12;
  result.matrix[1][1] = 13;
  result.matrix[1][2] = 22;
  result.matrix[2][0] = 15;
  result.matrix[2][1] = 15;
  result.matrix[2][2] = 27;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_2) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = S21_NAN;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = -1;
  B.matrix[0][2] = 1;
  B.matrix[1][0] = 2;
  B.matrix[1][1] = 3;
  B.matrix[1][2] = 4;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(mult_matrix_3) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = S21_INFINITY;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = -1;
  B.matrix[0][2] = 1;
  B.matrix[1][0] = 2;
  B.matrix[1][1] = 3;
  B.matrix[1][2] = 4;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(mult_matrix_4) {
  int res, s21_res;
  matrix_t A, B = {.rows = 0, .columns = 0}, s21_result;
  s21_create_matrix(3, 2, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = S21_INFINITY;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;

  res = INCORRECT_MATRIX;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(mult_matrix_5) {
  int res, s21_res;
  matrix_t A, B = {.rows = 3, .columns = 0}, s21_result;
  s21_create_matrix(3, 2, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;

  res = INCORRECT_MATRIX;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(mult_matrix_5_5) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(3, 1, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(mult_matrix_6) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = MAX_DOUBLE;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = -1;
  B.matrix[0][2] = 1;
  B.matrix[1][0] = 2;
  B.matrix[1][1] = 3;
  B.matrix[1][2] = 4;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_7) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 1.7976931348623158e+250;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = -1;
  B.matrix[0][2] = 1;
  B.matrix[1][0] = 2;
  B.matrix[1][1] = 3;
  B.matrix[1][2] = 4;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_8) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = -1;
  B.matrix[1][0] = 2;
  B.matrix[1][1] = 3;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_9) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &B);

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_10) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_create_matrix(3, 3, &result);

  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = -1;
  B.matrix[0][2] = 1;
  B.matrix[1][0] = -38;
  B.matrix[1][1] = 41;
  B.matrix[1][2] = -34;
  B.matrix[2][0] = 27;
  B.matrix[2][1] = -29;
  B.matrix[2][2] = 24;
  result.matrix[0][0] = 1;
  result.matrix[0][1] = 0;
  result.matrix[0][2] = 0;
  result.matrix[1][0] = 0;
  result.matrix[1][1] = 1;
  result.matrix[1][2] = 0;
  result.matrix[2][0] = 0;
  result.matrix[2][1] = 0;
  result.matrix[2][2] = 1;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_11_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_create_matrix(5, 5, &result);
  A.matrix[0][0] = 231;
  A.matrix[0][1] = 424;
  A.matrix[0][2] = 885;
  A.matrix[0][3] = -582;
  A.matrix[0][4] = 50;
  A.matrix[1][0] = -184;
  A.matrix[1][1] = 219;
  A.matrix[1][2] = -169;
  A.matrix[1][3] = 653;
  A.matrix[1][4] = -968;
  A.matrix[2][0] = -853;
  A.matrix[2][1] = -34;
  A.matrix[2][2] = 27;
  A.matrix[2][3] = -790;
  A.matrix[2][4] = -568;
  A.matrix[3][0] = -655;
  A.matrix[3][1] = 937;
  A.matrix[3][2] = 327;
  A.matrix[3][3] = 345;
  A.matrix[3][4] = 20;
  A.matrix[4][0] = 16;
  A.matrix[4][1] = 833;
  A.matrix[4][2] = -784;
  A.matrix[4][3] = 693;
  A.matrix[4][4] = 814;
  B.matrix[0][0] = -488;
  B.matrix[0][1] = -967;
  B.matrix[0][2] = 558;
  B.matrix[0][3] = -152;
  B.matrix[0][4] = 764;
  B.matrix[1][0] = 815;
  B.matrix[1][1] = 872;
  B.matrix[1][2] = 146;
  B.matrix[1][3] = 239;
  B.matrix[1][4] = -695;
  B.matrix[2][0] = 390;
  B.matrix[2][1] = 121;
  B.matrix[2][2] = 760;
  B.matrix[2][3] = 833;
  B.matrix[2][4] = -211;
  B.matrix[3][0] = -739;
  B.matrix[3][1] = -991;
  B.matrix[3][2] = 212;
  B.matrix[3][3] = 127;
  B.matrix[3][4] = -235;
  B.matrix[4][0] = -567;
  B.matrix[4][1] = 449;
  B.matrix[4][2] = 980;
  B.matrix[4][3] = -546;
  B.matrix[4][4] = 858;
  result.matrix[0][0] = 979730;
  result.matrix[0][1] = 852648;
  result.matrix[0][2] = 789018;
  result.matrix[0][3] = 702215;
  result.matrix[0][4] = -125261;
  result.matrix[1][0] = 268656;
  result.matrix[1][1] = -733308;
  result.matrix[1][2] = -1009342;
  result.matrix[1][3] = 550991;
  result.matrix[1][4] = -1241121;
  result.matrix[2][0] = 1304950;
  result.matrix[2][1] = 1326328;
  result.matrix[2][2] = -1184538;
  result.matrix[2][3] = 353819;
  result.matrix[2][4] = -935453;
  result.matrix[3][0] = 944530;
  result.matrix[3][1] = 1157101;
  result.matrix[3][2] = 112572;
  result.matrix[3][3] = 628789;
  result.matrix[3][4] = -1284547;
  result.matrix[4][0] = -608338;
  result.matrix[4][1] = 294763;
  result.matrix[4][2] = 479342;
  result.matrix[4][3] = -812850;
  result.matrix[4][4] = 134270;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[1][4], s21_result.matrix[1][4]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[2][4], s21_result.matrix[2][4]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[3][4], s21_result.matrix[3][4]);
  ck_assert_double_eq(result.matrix[4][0], s21_result.matrix[4][0]);
  ck_assert_double_eq(result.matrix[4][1], s21_result.matrix[4][1]);
  ck_assert_double_eq(result.matrix[4][2], s21_result.matrix[4][2]);
  ck_assert_double_eq(result.matrix[4][3], s21_result.matrix[4][3]);
  ck_assert_double_eq(result.matrix[4][4], s21_result.matrix[4][4]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_12_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(6, 6, &A);
  s21_create_matrix(6, 6, &B);
  s21_create_matrix(6, 6, &result);
  A.matrix[0][0] = 5.6095471332361155;
  A.matrix[0][1] = 37.43952600833608;
  A.matrix[0][2] = -18.127315142480626;
  A.matrix[0][3] = 35.335318165995986;
  A.matrix[0][4] = 11.101190391221675;
  A.matrix[0][5] = -3.4336630433232695;
  A.matrix[1][0] = 13.620604241139493;
  A.matrix[1][1] = 7.206036983722099;
  A.matrix[1][2] = 32.20890163947585;
  A.matrix[1][3] = 45.334490176220804;
  A.matrix[1][4] = -10.635663129603422;
  A.matrix[1][5] = 20.14967883234523;
  A.matrix[2][0] = -40.19820943453073;
  A.matrix[2][1] = 31.756986368568253;
  A.matrix[2][2] = -13.610979043515925;
  A.matrix[2][3] = -12.949001397465068;
  A.matrix[2][4] = -29.00305971473926;
  A.matrix[2][5] = -32.11502053565622;
  A.matrix[3][0] = 5.548763684815239;
  A.matrix[3][1] = -23.91391143816388;
  A.matrix[3][2] = -36.866000582825045;
  A.matrix[3][3] = 33.36508045012911;
  A.matrix[3][4] = -24.58359182369753;
  A.matrix[3][5] = -22.668180366194925;
  A.matrix[4][0] = -1.1838729040141722;
  A.matrix[4][1] = -47.663411205910855;
  A.matrix[4][2] = 43.350216691268294;
  A.matrix[4][3] = 33.6977866723797;
  A.matrix[4][4] = 48.624268605433386;
  A.matrix[4][5] = -10.514735240560404;
  A.matrix[5][0] = -25.24121630199927;
  A.matrix[5][1] = 6.075824986142404;
  A.matrix[5][2] = -8.672306619346884;
  A.matrix[5][3] = -19.360795477815625;
  A.matrix[5][4] = 43.54131259525847;
  A.matrix[5][5] = -33.40211002484375;
  B.matrix[0][0] = -49.42488480378568;
  B.matrix[0][1] = 30.469199286849815;
  B.matrix[0][2] = 16.495777683165727;
  B.matrix[0][3] = -41.291519599139896;
  B.matrix[0][4] = 30.617748005537624;
  B.matrix[0][5] = 40.320033916304695;
  B.matrix[1][0] = 49.89407023120788;
  B.matrix[1][1] = 10.814704765562162;
  B.matrix[1][2] = -18.573033101451003;
  B.matrix[1][3] = -5.792173280751055;
  B.matrix[1][4] = -1.7991482562525305;
  B.matrix[1][5] = -7.151113545802361;
  B.matrix[2][0] = 26.480500354702233;
  B.matrix[2][1] = 12.41497222551098;
  B.matrix[2][2] = -42.964311633381506;
  B.matrix[2][3] = 41.87500161688772;
  B.matrix[2][4] = -18.565650534022215;
  B.matrix[2][5] = -8.29394008788055;
  B.matrix[3][0] = -49.90484362433498;
  B.matrix[3][1] = -38.64116213279413;
  B.matrix[3][2] = -5.9587143801188205;
  B.matrix[3][3] = -46.19883926658897;
  B.matrix[3][4] = -37.92095444282739;
  B.matrix[3][5] = 48.62247989187037;
  B.matrix[4][0] = 40.04254157114378;
  B.matrix[4][1] = -43.22212183750769;
  B.matrix[4][2] = 32.954233845861076;
  B.matrix[4][3] = -25.949455797437324;
  B.matrix[4][4] = 33.23371240892765;
  B.matrix[4][5] = 46.46137678670543;
  B.matrix[5][0] = 7.3799573003452785;
  B.matrix[5][1] = -46.81164373400299;
  B.matrix[5][2] = -43.26849878492067;
  B.matrix[5][3] = -26.129541288557;
  B.matrix[5][4] = -41.9443680220865;
  B.matrix[5][5] = 31.658666700271866;
  result.matrix[0][0] = -233.48519224278047;
  result.matrix[0][1] = -1333.7136349729008;
  result.matrix[0][2] = 479.84350427150486;
  result.matrix[0][3] = -3038.364791964061;
  result.matrix[0][4] = -386.0545545107691;
  result.matrix[0][5] = 2233.9518875232957;
  result.matrix[1][0] = -2000.3363055717496;
  result.matrix[1][1] = -1342.5083575432543;
  result.matrix[1][2] = -2785.4605538451297;
  result.matrix[1][3] = -1600.319273247961;
  result.matrix[1][4] = -3111.6669633180773;
  result.matrix[1][5] = 2578.5530870144657;
  result.matrix[2][0] = 2458.705828696442;
  result.matrix[2][1] = 2206.950393731068;
  result.matrix[2][2] = -157.18342408305023;
  result.matrix[2][3] = 3095.9366206920763;
  result.matrix[2][4] = -161.01410612210393;
  result.matrix[2][5] = -4728.8557150974;
  result.matrix[3][0] = -5260.3983473585395;
  result.matrix[3][1] = 287.17842237397593;
  result.matrix[3][2] = 2091.479094833976;
  result.matrix[3][3] = -1945.5552042322436;
  result.matrix[3][4] = -234.08061299123978;
  result.matrix[3][5] = 462.96289996824856;
  result.matrix[4][0] = -983.9151502654364;
  result.matrix[4][1] = -2924.8993037712394;
  result.matrix[4][2] = 859.7498316919222;
  result.matrix[4][3] = -403.5776882554062;
  result.matrix[4][4] = 23.827697533112023;
  result.matrix[4][5] = 3498.3064909481;
  result.matrix[5][0] = 3784.240972791308;
  result.matrix[5][1] = -381.25444855822195;
  result.matrix[5][2] = 2838.8748976849542;
  result.matrix[5][3] = 1281.257817608506;
  result.matrix[5][4] = 2959.4962035871413;
  result.matrix[5][5] = -965.0948471181445;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq_tol(result.matrix[0][0], s21_result.matrix[0][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][1], s21_result.matrix[0][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][2], s21_result.matrix[0][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][3], s21_result.matrix[0][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][4], s21_result.matrix[0][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][5], s21_result.matrix[0][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][0], s21_result.matrix[1][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][1], s21_result.matrix[1][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][2], s21_result.matrix[1][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][3], s21_result.matrix[1][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][4], s21_result.matrix[1][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][5], s21_result.matrix[1][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][0], s21_result.matrix[2][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][1], s21_result.matrix[2][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][2], s21_result.matrix[2][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][3], s21_result.matrix[2][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][4], s21_result.matrix[2][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][5], s21_result.matrix[2][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][0], s21_result.matrix[3][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][1], s21_result.matrix[3][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][2], s21_result.matrix[3][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][3], s21_result.matrix[3][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][4], s21_result.matrix[3][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][5], s21_result.matrix[3][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][0], s21_result.matrix[4][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][1], s21_result.matrix[4][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][2], s21_result.matrix[4][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][3], s21_result.matrix[4][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][4], s21_result.matrix[4][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][5], s21_result.matrix[4][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][0], s21_result.matrix[5][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][1], s21_result.matrix[5][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][2], s21_result.matrix[5][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][3], s21_result.matrix[5][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][4], s21_result.matrix[5][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][5], s21_result.matrix[5][5], 1e-7);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_13_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(10, 10, &A);
  s21_create_matrix(10, 10, &B);
  s21_create_matrix(10, 10, &result);
  A.matrix[0][0] = 17.25016941485092;
  A.matrix[0][1] = -41.661805396118126;
  A.matrix[0][2] = -2.7994533308625997;
  A.matrix[0][3] = -32.87709695259722;
  A.matrix[0][4] = 28.244911449672163;
  A.matrix[0][5] = -44.17217562351856;
  A.matrix[0][6] = 22.42478289854417;
  A.matrix[0][7] = 0.7344719545392326;
  A.matrix[0][8] = -40.18088702639143;
  A.matrix[0][9] = 44.21985455899666;
  A.matrix[1][0] = 11.162153442919445;
  A.matrix[1][1] = -2.47606189944843;
  A.matrix[1][2] = -44.4346604758096;
  A.matrix[1][3] = -25.17318703395725;
  A.matrix[1][4] = 24.242519619208565;
  A.matrix[1][5] = -36.896675107030184;
  A.matrix[1][6] = -14.04959023909691;
  A.matrix[1][7] = -25.522084705932347;
  A.matrix[1][8] = 30.08727254935549;
  A.matrix[1][9] = 34.209515297199175;
  A.matrix[2][0] = -11.65722353170118;
  A.matrix[2][1] = 23.391828930537667;
  A.matrix[2][2] = 13.382070670033134;
  A.matrix[2][3] = -46.450072697520774;
  A.matrix[2][4] = 39.54793156687638;
  A.matrix[2][5] = -3.227634385894966;
  A.matrix[2][6] = -31.978107154355694;
  A.matrix[2][7] = 31.196422796163198;
  A.matrix[2][8] = 25.684653721668;
  A.matrix[2][9] = 13.269610025315393;
  A.matrix[3][0] = -1.2287645001028833;
  A.matrix[3][1] = -39.20866627808692;
  A.matrix[3][2] = -44.648037547884044;
  A.matrix[3][3] = -21.743674697102403;
  A.matrix[3][4] = -46.070671186061034;
  A.matrix[3][5] = 15.338060186118428;
  A.matrix[3][6] = -38.25406689467753;
  A.matrix[3][7] = 17.5234970487155;
  A.matrix[3][8] = -48.88411103630185;
  A.matrix[3][9] = -27.78783311546633;
  A.matrix[4][0] = 13.305669509980165;
  A.matrix[4][1] = 42.1081195177305;
  A.matrix[4][2] = -15.472879187686841;
  A.matrix[4][3] = 27.520233977773536;
  A.matrix[4][4] = 38.178929570696276;
  A.matrix[4][5] = 3.101139382002807;
  A.matrix[4][6] = 44.31393716189071;
  A.matrix[4][7] = -24.15267074421286;
  A.matrix[4][8] = -11.517826208484077;
  A.matrix[4][9] = -47.68188962730429;
  A.matrix[5][0] = -0.728729069374363;
  A.matrix[5][1] = -1.8527831119884484;
  A.matrix[5][2] = -17.532045495947283;
  A.matrix[5][3] = 8.810135451633647;
  A.matrix[5][4] = -26.04553786458243;
  A.matrix[5][5] = -0.5167408480910319;
  A.matrix[5][6] = 8.121600791713034;
  A.matrix[5][7] = 41.10345417872734;
  A.matrix[5][8] = -25.573354759766033;
  A.matrix[5][9] = 29.97663671147039;
  A.matrix[6][0] = -28.68646290135374;
  A.matrix[6][1] = -45.52976527954842;
  A.matrix[6][2] = 10.145214232880646;
  A.matrix[6][3] = 14.048528433361913;
  A.matrix[6][4] = -6.368683752636591;
  A.matrix[6][5] = 39.30176691253265;
  A.matrix[6][6] = 8.010489518644711;
  A.matrix[6][7] = 10.670884641167213;
  A.matrix[6][8] = 29.923642090488055;
  A.matrix[6][9] = -7.697963609341388;
  A.matrix[7][0] = -37.1966593074431;
  A.matrix[7][1] = -19.338950681607468;
  A.matrix[7][2] = 13.452815244311042;
  A.matrix[7][3] = 22.459334774054888;
  A.matrix[7][4] = -21.875254001951358;
  A.matrix[7][5] = 38.65982059256684;
  A.matrix[7][6] = -11.503229649694456;
  A.matrix[7][7] = -38.121316937288626;
  A.matrix[7][8] = -6.721061354391332;
  A.matrix[7][9] = -37.8695179577359;
  A.matrix[8][0] = 26.878946896587028;
  A.matrix[8][1] = -33.95931901488449;
  A.matrix[8][2] = -21.176995681413203;
  A.matrix[8][3] = -24.48546453364701;
  A.matrix[8][4] = -2.3838105663006908;
  A.matrix[8][5] = 19.304951637560222;
  A.matrix[8][6] = -12.110458721179377;
  A.matrix[8][7] = 43.613153171391986;
  A.matrix[8][8] = 22.004749204268364;
  A.matrix[8][9] = 5.641667100031599;
  A.matrix[9][0] = 20.416830110738314;
  A.matrix[9][1] = 26.30729509483448;
  A.matrix[9][2] = 36.67642195100953;
  A.matrix[9][3] = 37.92677030835209;
  A.matrix[9][4] = -38.23705307501064;
  A.matrix[9][5] = 21.679695911062236;
  A.matrix[9][6] = -0.11142794490972752;
  A.matrix[9][7] = 1.5188297196060851;
  A.matrix[9][8] = -45.34887016887828;
  A.matrix[9][9] = -23.749907413741447;
  B.matrix[0][0] = 33.81573136633548;
  B.matrix[0][1] = 41.84829408642901;
  B.matrix[0][2] = -8.079110430829425;
  B.matrix[0][3] = -26.29947484618072;
  B.matrix[0][4] = -5.922204757900451;
  B.matrix[0][5] = 37.449785547534454;
  B.matrix[0][6] = 20.848778536266977;
  B.matrix[0][7] = -11.60153744705875;
  B.matrix[0][8] = -23.36773176417686;
  B.matrix[0][9] = -14.878934237426865;
  B.matrix[1][0] = 44.14095354683313;
  B.matrix[1][1] = 9.732505119952732;
  B.matrix[1][2] = 17.622766334517173;
  B.matrix[1][3] = -35.20869942345994;
  B.matrix[1][4] = -26.234065712173553;
  B.matrix[1][5] = 41.322036248909185;
  B.matrix[1][6] = 23.209979342773423;
  B.matrix[1][7] = -19.88982321839407;
  B.matrix[1][8] = 6.490351142327893;
  B.matrix[1][9] = -19.625237904182484;
  B.matrix[2][0] = -41.1579549077362;
  B.matrix[2][1] = -16.267381456034474;
  B.matrix[2][2] = -6.6297331165394935;
  B.matrix[2][3] = -39.8619649876964;
  B.matrix[2][4] = -43.80110436534219;
  B.matrix[2][5] = -35.12121349897421;
  B.matrix[2][6] = -45.41478599882156;
  B.matrix[2][7] = -44.440586304466656;
  B.matrix[2][8] = -10.97724334968674;
  B.matrix[2][9] = -38.2844521651874;
  B.matrix[3][0] = -12.265402953504989;
  B.matrix[3][1] = -23.140314900916298;
  B.matrix[3][2] = 43.442725252354556;
  B.matrix[3][3] = -7.885143345050034;
  B.matrix[3][4] = 26.46616624722924;
  B.matrix[3][5] = -15.362406017285867;
  B.matrix[3][6] = -32.94114589965665;
  B.matrix[3][7] = -12.821940795681812;
  B.matrix[3][8] = -29.015102721175154;
  B.matrix[3][9] = -38.92816775362622;
  B.matrix[4][0] = -15.954991997307339;
  B.matrix[4][1] = 39.82411727725801;
  B.matrix[4][2] = -11.971434651600653;
  B.matrix[4][3] = 25.647538289591225;
  B.matrix[4][4] = -39.74055550764352;
  B.matrix[4][5] = -48.33105796062736;
  B.matrix[4][6] = -43.055028558396614;
  B.matrix[4][7] = -13.188456852367867;
  B.matrix[4][8] = 47.1040169444923;
  B.matrix[4][9] = 17.726049451888922;
  B.matrix[5][0] = -11.825405843267468;
  B.matrix[5][1] = -2.5596141709228952;
  B.matrix[5][2] = -42.95357566106056;
  B.matrix[5][3] = -23.20858497257575;
  B.matrix[5][4] = -0.43972985455947433;
  B.matrix[5][5] = -39.550516696009296;
  B.matrix[5][6] = -30.989281344642876;
  B.matrix[5][7] = 8.276401197966964;
  B.matrix[5][8] = -17.452459573718237;
  B.matrix[5][9] = 8.480709870294527;
  B.matrix[6][0] = 29.129394637446374;
  B.matrix[6][1] = 7.080641540550491;
  B.matrix[6][2] = -43.13316947070807;
  B.matrix[6][3] = -16.810734233818557;
  B.matrix[6][4] = -35.035371807729796;
  B.matrix[6][5] = 30.449620433974964;
  B.matrix[6][6] = -0.43854689783140477;
  B.matrix[6][7] = 28.823506281344073;
  B.matrix[6][8] = -48.28799213228038;
  B.matrix[6][9] = 48.728039122598766;
  B.matrix[7][0] = 8.646497086936678;
  B.matrix[7][1] = 35.47078388725123;
  B.matrix[7][2] = 43.35125537517978;
  B.matrix[7][3] = -34.82842583378818;
  B.matrix[7][4] = -15.554119286005976;
  B.matrix[7][5] = -31.636398942204277;
  B.matrix[7][6] = -3.071550419825733;
  B.matrix[7][7] = 16.39328398744643;
  B.matrix[7][8] = 40.804203140279505;
  B.matrix[7][9] = 2.753916942414987;
  B.matrix[8][0] = 48.79817794401711;
  B.matrix[8][1] = -25.012082123667806;
  B.matrix[8][2] = 9.111333567548831;
  B.matrix[8][3] = -34.60519763900569;
  B.matrix[8][4] = -25.365763915813524;
  B.matrix[8][5] = 43.13060906249754;
  B.matrix[8][6] = -27.59291285561824;
  B.matrix[8][7] = -27.674019198752013;
  B.matrix[8][8] = 46.76596472568405;
  B.matrix[8][9] = -20.125363147448876;
  B.matrix[9][0] = 23.71256213291262;
  B.matrix[9][1] = 1.4075956626242676;
  B.matrix[9][2] = -14.666915921948954;
  B.matrix[9][3] = -18.325683128114438;
  B.matrix[9][4] = 3.1313851289074615;
  B.matrix[9][5] = 7.236485287816146;
  B.matrix[9][6] = 49.53512678382398;
  B.matrix[9][7] = 31.338250230434134;
  B.matrix[9][8] = -16.81273857632973;
  B.matrix[9][9] = -34.165689741072406;
  result.matrix[0][0] = -918.1046041014649;
  result.matrix[0][1] = 3612.720439730745;
  result.matrix[0][2] = -2674.1346041500956;
  result.matrix[0][3] = 3311.156985870982;
  result.matrix[0][4] = -499.153089751363;
  result.matrix[0][5] = -843.6577835941653;
  result.matrix[0][6] = 4042.6531790392164;
  result.matrix[0][7] = 3592.5223370488197;
  result.matrix[0][8] = -1262.9055526764766;
  result.matrix[0][9] = 2466.630164874007;
  result.matrix[1][0] = 4104.755281804881;
  result.matrix[1][1] = 1099.087408773979;
  result.matrix[1][2] = -366.2148153953606;
  result.matrix[1][3] = 2698.4369950243486;
  result.matrix[1][4] = 564.858116109614;
  result.matrix[1][5] = 4475.50140222753;
  result.matrix[1][6] = 4071.042176057693;
  result.matrix[1][7] = 1008.210017259075;
  result.matrix[1][8] = 3196.048618944885;
  result.matrix[1][9] = 151.22548656466083;
  result.matrix[2][0] = 970.7290780029748;
  result.matrix[2][1] = 2436.61483534538;
  result.matrix[2][2] = 836.082787548592;
  result.matrix[2][3] = -1275.9152754433094;
  result.matrix[2][4] = -3905.19522599763;
  result.matrix[2][5] = -1766.9545741064983;
  result.matrix[2][6] = -513.6555951283046;
  result.matrix[2][7] = -1582.6945923933038;
  result.matrix[2][8] = 7339.448923238498;
  result.matrix[2][9] = -758.6726696551815;
  result.matrix[3][0] = -3121.441448714525;
  result.matrix[3][1] = 456.74677365484547;
  result.matrix[3][2] = 934.9168641107657;
  result.matrix[3][3] = 4060.0770146931004;
  result.matrix[3][4] = 6460.822133955466;
  result.matrix[3][5] = -2172.7442927839347;
  result.matrix[3][6] = 3251.883736096103;
  result.matrix[3][7] = 3458.2841586218055;
  result.matrix[3][8] = -799.2364201912213;
  result.matrix[3][9] = 2774.3693851657285;
  result.matrix[4][0] = 1351.3955136720929;
  result.matrix[4][1] = 1772.042789493972;
  result.matrix[4][2] = -1021.6086339869975;
  result.matrix[4][3] = 843.1261696152362;
  result.matrix[4][4] = -2330.0309659382087;
  result.matrix[4][5] = 1662.6832147439588;
  result.matrix[4][6] = -2678.3776551260967;
  result.matrix[4][7] = -1429.1641451959335;
  result.matrix[4][8] = -1784.3618958315114;
  result.matrix[4][9] = 3153.465727235642;
  result.matrix[5][0] = 983.6321310182961;
  result.matrix[5][1] = 1194.2006936964026;
  result.matrix[5][2] = 1565.1074965424832;
  result.matrix[5][3] = -1174.6892551657586;
  result.matrix[5][4] = 1907.9912913235894;
  result.matrix[5][5] = -283.3384044310458;
  result.matrix[5][6] = 3645.933586973754;
  result.matrix[5][7] = 3605.747411270377;
  result.matrix[5][8] = -1690.935045091598;
  result.matrix[5][9] = -91.1747117423669;
  result.matrix[6][0] = -2329.507900908044;
  result.matrix[6][1] = -2812.010385164651;
  result.matrix[6][2] = -1136.8340394319844;
  result.matrix[6][3] = -633.9349938115989;
  result.matrix[6][4] = 297.80191206092263;
  result.matrix[6][5] = -3633.169182653901;
  result.matrix[6][6] = -4765.35882050055;
  result.matrix[6][7] = 353.1407709755484;
  result.matrix[6][8] = 447.3826695248056;
  result.matrix[6][9] = 685.9885584380266;
  result.matrix[7][0] = -4939.44050000074;
  result.matrix[7][1] = -4772.348707061755;
  result.matrix[7][2] = -1214.728098587039;
  result.matrix[7][3] = 1935.166676221408;
  result.matrix[7][4] = 2632.9889767011928;
  result.matrix[7][5] = -3189.5755270907894;
  result.matrix[7][6] = -4399.636338633409;
  result.matrix[7][7] = -1418.4341180343129;
  result.matrix[7][8] = -2438.4397570610226;
  result.matrix[7][9] = 247.3338667700683;
  result.matrix[8][0] = 1623.5059569472587;
  result.matrix[8][1] = 2479.8775812057834;
  result.matrix[8][2] = -8.817509041930748;
  result.matrix[8][3] = -1163.4480235052768;
  result.matrix[8][4] = 302.9219912397158;
  result.matrix[8][5] = -683.6666132121283;
  result.matrix[8][6] = 588.5511323058486;
  result.matrix[8][7] = 1743.6283772385746;
  result.matrix[8][8] = 2943.8105841881184;
  result.matrix[8][9] = 1046.3008090663407;
  result.matrix[9][0] = -2535.6004982825993;
  result.matrix[9][1] = -788.1454046849368;
  result.matrix[9][2] = 1235.4754134290363;
  result.matrix[9][3] = -2754.576541499052;
  result.matrix[9][4] = 1152.4954750837721;
  result.matrix[9][5] = -1207.7294139063242;
  result.matrix[9][6] = -834.0497150082367;
  result.matrix[9][7] = -1660.2211554711537;
  result.matrix[9][8] = -5643.016794441232;
  result.matrix[9][9] = -2471.7091645098694;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq_tol(result.matrix[0][0], s21_result.matrix[0][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][1], s21_result.matrix[0][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][2], s21_result.matrix[0][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][3], s21_result.matrix[0][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][4], s21_result.matrix[0][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][5], s21_result.matrix[0][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][6], s21_result.matrix[0][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][7], s21_result.matrix[0][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][8], s21_result.matrix[0][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][9], s21_result.matrix[0][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][0], s21_result.matrix[1][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][1], s21_result.matrix[1][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][2], s21_result.matrix[1][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][3], s21_result.matrix[1][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][4], s21_result.matrix[1][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][5], s21_result.matrix[1][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][6], s21_result.matrix[1][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][7], s21_result.matrix[1][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][8], s21_result.matrix[1][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][9], s21_result.matrix[1][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][0], s21_result.matrix[2][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][1], s21_result.matrix[2][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][2], s21_result.matrix[2][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][3], s21_result.matrix[2][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][4], s21_result.matrix[2][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][5], s21_result.matrix[2][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][6], s21_result.matrix[2][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][7], s21_result.matrix[2][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][8], s21_result.matrix[2][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][9], s21_result.matrix[2][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][0], s21_result.matrix[3][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][1], s21_result.matrix[3][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][2], s21_result.matrix[3][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][3], s21_result.matrix[3][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][4], s21_result.matrix[3][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][5], s21_result.matrix[3][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][6], s21_result.matrix[3][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][7], s21_result.matrix[3][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][8], s21_result.matrix[3][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][9], s21_result.matrix[3][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][0], s21_result.matrix[4][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][1], s21_result.matrix[4][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][2], s21_result.matrix[4][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][3], s21_result.matrix[4][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][4], s21_result.matrix[4][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][5], s21_result.matrix[4][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][6], s21_result.matrix[4][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][7], s21_result.matrix[4][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][8], s21_result.matrix[4][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[4][9], s21_result.matrix[4][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][0], s21_result.matrix[5][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][1], s21_result.matrix[5][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][2], s21_result.matrix[5][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][3], s21_result.matrix[5][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][4], s21_result.matrix[5][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][5], s21_result.matrix[5][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][6], s21_result.matrix[5][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][7], s21_result.matrix[5][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][8], s21_result.matrix[5][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[5][9], s21_result.matrix[5][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][0], s21_result.matrix[6][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][1], s21_result.matrix[6][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][2], s21_result.matrix[6][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][3], s21_result.matrix[6][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][4], s21_result.matrix[6][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][5], s21_result.matrix[6][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][6], s21_result.matrix[6][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][7], s21_result.matrix[6][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][8], s21_result.matrix[6][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[6][9], s21_result.matrix[6][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][0], s21_result.matrix[7][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][1], s21_result.matrix[7][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][2], s21_result.matrix[7][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][3], s21_result.matrix[7][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][4], s21_result.matrix[7][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][5], s21_result.matrix[7][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][6], s21_result.matrix[7][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][7], s21_result.matrix[7][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][8], s21_result.matrix[7][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[7][9], s21_result.matrix[7][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][0], s21_result.matrix[8][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][1], s21_result.matrix[8][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][2], s21_result.matrix[8][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][3], s21_result.matrix[8][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][4], s21_result.matrix[8][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][5], s21_result.matrix[8][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][6], s21_result.matrix[8][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][7], s21_result.matrix[8][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][8], s21_result.matrix[8][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[8][9], s21_result.matrix[8][9], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][0], s21_result.matrix[9][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][1], s21_result.matrix[9][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][2], s21_result.matrix[9][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][3], s21_result.matrix[9][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][4], s21_result.matrix[9][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][5], s21_result.matrix[9][5], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][6], s21_result.matrix[9][6], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][7], s21_result.matrix[9][7], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][8], s21_result.matrix[9][8], 1e-7);
  ck_assert_double_eq_tol(result.matrix[9][9], s21_result.matrix[9][9], 1e-7);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_14_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(2, 4, &A);
  s21_create_matrix(4, 3, &B);
  s21_create_matrix(2, 3, &result);
  A.matrix[0][0] = 10.742662612833183;
  A.matrix[0][1] = 30.155782744556753;
  A.matrix[0][2] = 41.84970815320477;
  A.matrix[0][3] = 36.94363535621355;
  A.matrix[1][0] = 38.20755232913024;
  A.matrix[1][1] = 22.95392816655287;
  A.matrix[1][2] = 33.66501540339439;
  A.matrix[1][3] = -47.08131916142866;
  B.matrix[0][0] = -22.66416635699748;
  B.matrix[0][1] = -15.326436092743212;
  B.matrix[0][2] = -10.770089780426698;
  B.matrix[1][0] = -47.02969634805594;
  B.matrix[1][1] = -39.68393659300598;
  B.matrix[1][2] = -37.993152074955994;
  B.matrix[2][0] = 4.037231933870871;
  B.matrix[2][1] = -29.2329830198149;
  B.matrix[2][2] = 28.97588773135331;
  B.matrix[3][0] = 35.09477266259352;
  B.matrix[3][1] = -14.299572500125487;
  B.matrix[3][2] = 34.1055350961101;
  result.matrix[0][0] = -196.20533585343634;
  result.matrix[0][1] = -3113.0169023701624;
  result.matrix[0][2] = 1211.202216682073;
  result.matrix[1][0] = -3461.853311133673;
  result.matrix[1][1] = -1807.3739258532933;
  result.matrix[1][2] = -1913.8507286527733;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq_tol(result.matrix[0][0], s21_result.matrix[0][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][1], s21_result.matrix[0][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][2], s21_result.matrix[0][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][0], s21_result.matrix[1][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][1], s21_result.matrix[1][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][2], s21_result.matrix[1][2], 1e-7);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_matrix_15_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(4, 6, &A);
  s21_create_matrix(6, 5, &B);
  s21_create_matrix(4, 5, &result);
  A.matrix[0][0] = -10.193039775492018;
  A.matrix[0][1] = -20.344964868133637;
  A.matrix[0][2] = -34.895408254277456;
  A.matrix[0][3] = 45.31661238201161;
  A.matrix[0][4] = -41.32328676916061;
  A.matrix[0][5] = 32.75963517166471;
  A.matrix[1][0] = -15.153274630715686;
  A.matrix[1][1] = 15.844063367374773;
  A.matrix[1][2] = -41.3312159768898;
  A.matrix[1][3] = 33.490315612360504;
  A.matrix[1][4] = -47.068464020964996;
  A.matrix[1][5] = -47.28368459886788;
  A.matrix[2][0] = 10.943856724951841;
  A.matrix[2][1] = -31.049912731392798;
  A.matrix[2][2] = 39.39491790777085;
  A.matrix[2][3] = 10.580651257133606;
  A.matrix[2][4] = -19.432715003476915;
  A.matrix[2][5] = -14.149769026701657;
  A.matrix[3][0] = -4.591768558728564;
  A.matrix[3][1] = -20.647992645709007;
  A.matrix[3][2] = 49.23555916539776;
  A.matrix[3][3] = -32.91438481886316;
  A.matrix[3][4] = -3.6706459221992107;
  A.matrix[3][5] = -48.71121269339168;
  B.matrix[0][0] = -8.224648728349749;
  B.matrix[0][1] = -5.663883027810261;
  B.matrix[0][2] = 7.161742735823504;
  B.matrix[0][3] = 41.37901221767817;
  B.matrix[0][4] = 28.731810935001274;
  B.matrix[1][0] = 4.218769825997606;
  B.matrix[1][1] = 26.085371394819447;
  B.matrix[1][2] = 36.87781359774895;
  B.matrix[1][3] = -11.80731791282085;
  B.matrix[1][4] = 2.2048899663491848;
  B.matrix[2][0] = -17.452754556817275;
  B.matrix[2][1] = -17.23839493010381;
  B.matrix[2][2] = -29.50567486308587;
  B.matrix[2][3] = 4.084524511673684;
  B.matrix[2][4] = -13.846553875484734;
  B.matrix[3][0] = 0.07756835526489836;
  B.matrix[3][1] = -41.10778463331286;
  B.matrix[3][2] = 35.11323452136778;
  B.matrix[3][3] = -19.94743483338521;
  B.matrix[3][4] = 49.252881339122105;
  B.matrix[4][0] = 4.442233761158124;
  B.matrix[4][1] = 3.6847502412205446;
  B.matrix[4][2] = 3.708115120871914;
  B.matrix[4][3] = -33.964457185542386;
  B.matrix[4][4] = -47.67142289755205;
  B.matrix[5][0] = 46.301420319017076;
  B.matrix[5][1] = -43.607948762655184;
  B.matrix[5][2] = -15.83019553115918;
  B.matrix[5][3] = -8.24776052441281;
  B.matrix[5][4] = 20.192623114642224;
  result.matrix[0][0] = 1943.7895162142638;
  result.matrix[0][1] = -3315.1449759301504;
  result.matrix[0][2] = 1125.7247230112305;
  result.matrix[0][3] = -94.71039352752491;
  result.matrix[0][4] = 5008.874827640314;
  result.matrix[1][0] = -1482.9767007476626;
  result.matrix[1][1] = 1723.404755935982;
  result.matrix[1][2] = 3445.2039770654146;
  result.matrix[1][3] = 337.67165603899684;
  result.matrix[1][4] = 3110.381875142775;
  result.matrix[2][0] = -1649.2099873235504;
  result.matrix[2][1] = -1440.5478146067883;
  result.matrix[2][2] = -1705.593691527527;
  result.matrix[2][3] = 1546.0403514408426;
  result.matrix[2][4] = 862.2831107745083;
  result.matrix[3][0] = -3182.8968895866306;
  result.matrix[3][1] = 2102.3627661137502;
  result.matrix[3][2] = -2645.299958935148;
  result.matrix[3][3] = 1437.8858740660096;
  result.matrix[3][4] = -3288.9497356139145;

  res = OK;
  s21_res = s21_mult_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq_tol(result.matrix[0][0], s21_result.matrix[0][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][1], s21_result.matrix[0][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][2], s21_result.matrix[0][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][3], s21_result.matrix[0][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[0][4], s21_result.matrix[0][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][0], s21_result.matrix[1][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][1], s21_result.matrix[1][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][2], s21_result.matrix[1][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][3], s21_result.matrix[1][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[1][4], s21_result.matrix[1][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][0], s21_result.matrix[2][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][1], s21_result.matrix[2][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][2], s21_result.matrix[2][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][3], s21_result.matrix[2][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[2][4], s21_result.matrix[2][4], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][0], s21_result.matrix[3][0], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][1], s21_result.matrix[3][1], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][2], s21_result.matrix[3][2], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][3], s21_result.matrix[3][3], 1e-7);
  ck_assert_double_eq_tol(result.matrix[3][4], s21_result.matrix[3][4], 1e-7);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_number_1) {
  int res, s21_res;
  double number = 2;
  matrix_t A, result, s21_result;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &result);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 4;
  result.matrix[0][0] = 2;
  result.matrix[0][1] = 4;
  result.matrix[0][2] = 6;
  result.matrix[1][0] = 0;
  result.matrix[1][1] = 8;
  result.matrix[1][2] = 4;
  result.matrix[2][0] = 4;
  result.matrix[2][1] = 6;
  result.matrix[2][2] = 8;
  result.rows = 3;
  result.columns = 3;

  res = OK;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_number_2) {
  int res, s21_res;
  double number = -2;
  matrix_t A, result, s21_result;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &result);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 4;
  result.matrix[0][0] = -2;
  result.matrix[0][1] = -4;
  result.matrix[0][2] = -6;
  result.matrix[1][0] = -0;
  result.matrix[1][1] = -8;
  result.matrix[1][2] = -4;
  result.matrix[2][0] = -4;
  result.matrix[2][1] = -6;
  result.matrix[2][2] = -8;
  result.rows = 3;
  result.columns = 3;

  res = OK;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_number_3) {
  int res, s21_res;
  double number = 2;
  matrix_t A, s21_result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = S21_INFINITY;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 4;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(mult_number_4) {
  int res, s21_res;
  double number = S21_INFINITY;
  matrix_t A, s21_result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 4;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(mult_number_5) {
  int res, s21_res;
  double number = 2;
  matrix_t A = {.columns = 0, .rows = 10}, s21_result;

  res = INCORRECT_MATRIX;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
}
END_TEST

START_TEST(mult_number_6) {
  int res, s21_res;
  double number = S21_NAN;
  matrix_t A, s21_result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 4;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(mult_number_7) {
  int res, s21_res;
  double number = MAX_DOUBLE;
  matrix_t A, s21_result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 4;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(mult_number_8) {
  int res, s21_res;
  double number = 5;
  matrix_t A, s21_result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = MIN_DOUBLE;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 4;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(mult_number_9) {
  int res, s21_res;
  double number = 5;
  matrix_t A, s21_result;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = S21_NAN;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 4;

  res = CALCULATION_ERROR;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(mult_number_10) {
  int res, s21_res;
  double number = 2;
  matrix_t A, result, s21_result;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &result);

  A.matrix[0][0] = 1;
  result.matrix[0][0] = 2;
  result.rows = 1;
  result.columns = 1;

  res = OK;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_number_11_by_py) {
  int res, s21_res;
  double number = 758.9435117173717;
  matrix_t A, result, s21_result;
  s21_create_matrix(7, 11, &A);
  s21_create_matrix(7, 11, &result);
  A.matrix[0][0] = 70;
  A.matrix[0][1] = -231;
  A.matrix[0][2] = -111;
  A.matrix[0][3] = 939;
  A.matrix[0][4] = 941;
  A.matrix[0][5] = -145;
  A.matrix[0][6] = 195;
  A.matrix[0][7] = 805;
  A.matrix[0][8] = 613;
  A.matrix[0][9] = 412;
  A.matrix[0][10] = 639;
  A.matrix[1][0] = 364;
  A.matrix[1][1] = -626;
  A.matrix[1][2] = 807;
  A.matrix[1][3] = -271;
  A.matrix[1][4] = -137;
  A.matrix[1][5] = -138;
  A.matrix[1][6] = -587;
  A.matrix[1][7] = -773;
  A.matrix[1][8] = -9;
  A.matrix[1][9] = -764;
  A.matrix[1][10] = 67;
  A.matrix[2][0] = -968;
  A.matrix[2][1] = -999;
  A.matrix[2][2] = -73;
  A.matrix[2][3] = 464;
  A.matrix[2][4] = 81;
  A.matrix[2][5] = -537;
  A.matrix[2][6] = -5;
  A.matrix[2][7] = 159;
  A.matrix[2][8] = -605;
  A.matrix[2][9] = -272;
  A.matrix[2][10] = -406;
  A.matrix[3][0] = 684;
  A.matrix[3][1] = -34;
  A.matrix[3][2] = -833;
  A.matrix[3][3] = 224;
  A.matrix[3][4] = -363;
  A.matrix[3][5] = 365;
  A.matrix[3][6] = -485;
  A.matrix[3][7] = 271;
  A.matrix[3][8] = 328;
  A.matrix[3][9] = -560;
  A.matrix[3][10] = -351;
  A.matrix[4][0] = -506;
  A.matrix[4][1] = 706;
  A.matrix[4][2] = -796;
  A.matrix[4][3] = 100;
  A.matrix[4][4] = 434;
  A.matrix[4][5] = 487;
  A.matrix[4][6] = 870;
  A.matrix[4][7] = -895;
  A.matrix[4][8] = 441;
  A.matrix[4][9] = 331;
  A.matrix[4][10] = 429;
  A.matrix[5][0] = -266;
  A.matrix[5][1] = 585;
  A.matrix[5][2] = 827;
  A.matrix[5][3] = 402;
  A.matrix[5][4] = 600;
  A.matrix[5][5] = -398;
  A.matrix[5][6] = -654;
  A.matrix[5][7] = 174;
  A.matrix[5][8] = -868;
  A.matrix[5][9] = 489;
  A.matrix[5][10] = 88;
  A.matrix[6][0] = 979;
  A.matrix[6][1] = 529;
  A.matrix[6][2] = -173;
  A.matrix[6][3] = 23;
  A.matrix[6][4] = 576;
  A.matrix[6][5] = -712;
  A.matrix[6][6] = -151;
  A.matrix[6][7] = 648;
  A.matrix[6][8] = -465;
  A.matrix[6][9] = -966;
  A.matrix[6][10] = 916;
  result.matrix[0][0] = 53126.04582021602;
  result.matrix[0][1] = -175315.95120671287;
  result.matrix[0][2] = -84242.72980062827;
  result.matrix[0][3] = 712647.957502612;
  result.matrix[0][4] = 714165.8445260468;
  result.matrix[0][5] = -110046.8091990189;
  result.matrix[0][6] = 147993.9847848875;
  result.matrix[0][7] = 610949.5269324842;
  result.matrix[0][8] = 465232.3726827489;
  result.matrix[0][9] = 312684.72682755714;
  result.matrix[0][10] = 484964.9039874005;
  result.matrix[1][0] = 276255.4382651233;
  result.matrix[1][1] = -475098.6383350747;
  result.matrix[1][2] = 612467.413955919;
  result.matrix[1][3] = -205673.69167540775;
  result.matrix[1][4] = -103975.26110527993;
  result.matrix[1][5] = -104734.2046169973;
  result.matrix[1][6] = -445499.8413780972;
  result.matrix[1][7] = -586663.3345575284;
  result.matrix[1][8] = -6830.4916054563455;
  result.matrix[1][9] = -579832.842952072;
  result.matrix[1][10] = 50849.2152850639;
  result.matrix[2][0] = -734657.3193424159;
  result.matrix[2][1] = -758184.5682056544;
  result.matrix[2][2] = -55402.876355368135;
  result.matrix[2][3] = 352149.7894368605;
  result.matrix[2][4] = 61474.42444910711;
  result.matrix[2][5] = -407552.6657922286;
  result.matrix[2][6] = -3794.7175585868586;
  result.matrix[2][7] = 120672.0183630621;
  result.matrix[2][8] = -459160.8245890099;
  result.matrix[2][9] = -206432.63518712512;
  result.matrix[2][10] = -308131.0657572529;
  result.matrix[3][0] = 519117.36201468226;
  result.matrix[3][1] = -25804.07939839064;
  result.matrix[3][2] = -632199.9452605706;
  result.matrix[3][3] = 170003.34662469127;
  result.matrix[3][4] = -275496.4947534059;
  result.matrix[3][5] = 277014.38177684066;
  result.matrix[3][6] = -368087.6031829253;
  result.matrix[3][7] = 205673.69167540775;
  result.matrix[3][8] = 248933.47184329794;
  result.matrix[3][9] = -425008.3665617282;
  result.matrix[3][10] = -266389.17261279747;
  result.matrix[4][0] = -384025.4169289901;
  result.matrix[4][1] = 535814.1192724644;
  result.matrix[4][2] = -604119.0353270279;
  result.matrix[4][3] = 75894.35117173717;
  result.matrix[4][4] = 329381.48408533935;
  result.matrix[4][5] = 369605.49020636006;
  result.matrix[4][6] = 660280.8551941135;
  result.matrix[4][7] = -679254.4429870477;
  result.matrix[4][8] = 334694.08866736095;
  result.matrix[4][9] = 251210.30237845005;
  result.matrix[4][10] = 325586.7665267525;
  result.matrix[5][0] = -201878.97411682087;
  result.matrix[5][1] = 443981.95435466245;
  result.matrix[5][2] = 627646.2841902664;
  result.matrix[5][3] = 305095.29171038343;
  result.matrix[5][4] = 455366.10703042307;
  result.matrix[5][5] = -302059.51766351395;
  result.matrix[5][6] = -496349.0566631611;
  result.matrix[5][7] = 132056.17103882268;
  result.matrix[5][8] = -658762.9681706787;
  result.matrix[5][9] = 371123.3772297948;
  result.matrix[5][10] = 66787.02903112871;
  result.matrix[6][0] = 743005.697971307;
  result.matrix[6][1] = 401481.11769848963;
  result.matrix[6][2] = -131297.2275271053;
  result.matrix[6][3] = 17455.70076949955;
  result.matrix[6][4] = 437151.4627492061;
  result.matrix[6][5] = -540367.7803427686;
  result.matrix[6][6] = -114600.47026932314;
  result.matrix[6][7] = 491795.39559285686;
  result.matrix[6][8] = -352908.73294857785;
  result.matrix[6][9] = -733139.432318981;
  result.matrix[6][10] = 695192.2567331125;

  res = OK;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);
  ck_assert_double_eq(result.matrix[0][5], s21_result.matrix[0][5]);
  ck_assert_double_eq(result.matrix[0][6], s21_result.matrix[0][6]);
  ck_assert_double_eq(result.matrix[0][7], s21_result.matrix[0][7]);
  ck_assert_double_eq(result.matrix[0][8], s21_result.matrix[0][8]);
  ck_assert_double_eq(result.matrix[0][9], s21_result.matrix[0][9]);
  ck_assert_double_eq(result.matrix[0][10], s21_result.matrix[0][10]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[1][4], s21_result.matrix[1][4]);
  ck_assert_double_eq(result.matrix[1][5], s21_result.matrix[1][5]);
  ck_assert_double_eq(result.matrix[1][6], s21_result.matrix[1][6]);
  ck_assert_double_eq(result.matrix[1][7], s21_result.matrix[1][7]);
  ck_assert_double_eq(result.matrix[1][8], s21_result.matrix[1][8]);
  ck_assert_double_eq(result.matrix[1][9], s21_result.matrix[1][9]);
  ck_assert_double_eq(result.matrix[1][10], s21_result.matrix[1][10]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[2][4], s21_result.matrix[2][4]);
  ck_assert_double_eq(result.matrix[2][5], s21_result.matrix[2][5]);
  ck_assert_double_eq(result.matrix[2][6], s21_result.matrix[2][6]);
  ck_assert_double_eq(result.matrix[2][7], s21_result.matrix[2][7]);
  ck_assert_double_eq(result.matrix[2][8], s21_result.matrix[2][8]);
  ck_assert_double_eq(result.matrix[2][9], s21_result.matrix[2][9]);
  ck_assert_double_eq(result.matrix[2][10], s21_result.matrix[2][10]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[3][4], s21_result.matrix[3][4]);
  ck_assert_double_eq(result.matrix[3][5], s21_result.matrix[3][5]);
  ck_assert_double_eq(result.matrix[3][6], s21_result.matrix[3][6]);
  ck_assert_double_eq(result.matrix[3][7], s21_result.matrix[3][7]);
  ck_assert_double_eq(result.matrix[3][8], s21_result.matrix[3][8]);
  ck_assert_double_eq(result.matrix[3][9], s21_result.matrix[3][9]);
  ck_assert_double_eq(result.matrix[3][10], s21_result.matrix[3][10]);
  ck_assert_double_eq(result.matrix[4][0], s21_result.matrix[4][0]);
  ck_assert_double_eq(result.matrix[4][1], s21_result.matrix[4][1]);
  ck_assert_double_eq(result.matrix[4][2], s21_result.matrix[4][2]);
  ck_assert_double_eq(result.matrix[4][3], s21_result.matrix[4][3]);
  ck_assert_double_eq(result.matrix[4][4], s21_result.matrix[4][4]);
  ck_assert_double_eq(result.matrix[4][5], s21_result.matrix[4][5]);
  ck_assert_double_eq(result.matrix[4][6], s21_result.matrix[4][6]);
  ck_assert_double_eq(result.matrix[4][7], s21_result.matrix[4][7]);
  ck_assert_double_eq(result.matrix[4][8], s21_result.matrix[4][8]);
  ck_assert_double_eq(result.matrix[4][9], s21_result.matrix[4][9]);
  ck_assert_double_eq(result.matrix[4][10], s21_result.matrix[4][10]);
  ck_assert_double_eq(result.matrix[5][0], s21_result.matrix[5][0]);
  ck_assert_double_eq(result.matrix[5][1], s21_result.matrix[5][1]);
  ck_assert_double_eq(result.matrix[5][2], s21_result.matrix[5][2]);
  ck_assert_double_eq(result.matrix[5][3], s21_result.matrix[5][3]);
  ck_assert_double_eq(result.matrix[5][4], s21_result.matrix[5][4]);
  ck_assert_double_eq(result.matrix[5][5], s21_result.matrix[5][5]);
  ck_assert_double_eq(result.matrix[5][6], s21_result.matrix[5][6]);
  ck_assert_double_eq(result.matrix[5][7], s21_result.matrix[5][7]);
  ck_assert_double_eq(result.matrix[5][8], s21_result.matrix[5][8]);
  ck_assert_double_eq(result.matrix[5][9], s21_result.matrix[5][9]);
  ck_assert_double_eq(result.matrix[5][10], s21_result.matrix[5][10]);
  ck_assert_double_eq(result.matrix[6][0], s21_result.matrix[6][0]);
  ck_assert_double_eq(result.matrix[6][1], s21_result.matrix[6][1]);
  ck_assert_double_eq(result.matrix[6][2], s21_result.matrix[6][2]);
  ck_assert_double_eq(result.matrix[6][3], s21_result.matrix[6][3]);
  ck_assert_double_eq(result.matrix[6][4], s21_result.matrix[6][4]);
  ck_assert_double_eq(result.matrix[6][5], s21_result.matrix[6][5]);
  ck_assert_double_eq(result.matrix[6][6], s21_result.matrix[6][6]);
  ck_assert_double_eq(result.matrix[6][7], s21_result.matrix[6][7]);
  ck_assert_double_eq(result.matrix[6][8], s21_result.matrix[6][8]);
  ck_assert_double_eq(result.matrix[6][9], s21_result.matrix[6][9]);
  ck_assert_double_eq(result.matrix[6][10], s21_result.matrix[6][10]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_number_12_by_py) {
  int res, s21_res;
  double number = -245.94864843839463;
  matrix_t A, result, s21_result;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &result);
  A.matrix[0][0] = -274;
  A.matrix[0][1] = 535;
  A.matrix[0][2] = 402;
  A.matrix[0][3] = 154;
  A.matrix[0][4] = 915;
  A.matrix[1][0] = 58;
  A.matrix[1][1] = -111;
  A.matrix[1][2] = -695;
  A.matrix[1][3] = -688;
  A.matrix[1][4] = -55;
  A.matrix[2][0] = -535;
  A.matrix[2][1] = 282;
  A.matrix[2][2] = 701;
  A.matrix[2][3] = 429;
  A.matrix[2][4] = -569;
  A.matrix[3][0] = 850;
  A.matrix[3][1] = -218;
  A.matrix[3][2] = 355;
  A.matrix[3][3] = 636;
  A.matrix[3][4] = 885;
  A.matrix[4][0] = 411;
  A.matrix[4][1] = 7;
  A.matrix[4][2] = 275;
  A.matrix[4][3] = 145;
  A.matrix[4][4] = 908;
  result.matrix[0][0] = 67389.92967212012;
  result.matrix[0][1] = -131582.52691454111;
  result.matrix[0][2] = -98871.35667223464;
  result.matrix[0][3] = -37876.09185951277;
  result.matrix[0][4] = -225043.0133211311;
  result.matrix[1][0] = -14265.021609426889;
  result.matrix[1][1] = 27300.299976661805;
  result.matrix[1][2] = 170934.31066468428;
  result.matrix[1][3] = 169212.6701256155;
  result.matrix[1][4] = 13527.175664111704;
  result.matrix[2][0] = 131582.52691454111;
  result.matrix[2][1] = -69357.51885962728;
  result.matrix[2][2] = -172410.00255531463;
  result.matrix[2][3] = -105511.9701800713;
  result.matrix[2][4] = 139944.78096144655;
  result.matrix[3][0] = -209056.35117263542;
  result.matrix[3][1] = 53616.805359570026;
  result.matrix[3][2] = -87311.77019563009;
  result.matrix[3][3] = -156423.34040681898;
  result.matrix[3][4] = -217664.55386797924;
  result.matrix[4][0] = -101084.8945081802;
  result.matrix[4][1] = -1721.6405390687623;
  result.matrix[4][2] = -67635.87832055852;
  result.matrix[4][3] = -35662.55402356722;
  result.matrix[4][4] = -223321.37278206233;

  res = OK;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[1][4], s21_result.matrix[1][4]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[2][4], s21_result.matrix[2][4]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[3][4], s21_result.matrix[3][4]);
  ck_assert_double_eq(result.matrix[4][0], s21_result.matrix[4][0]);
  ck_assert_double_eq(result.matrix[4][1], s21_result.matrix[4][1]);
  ck_assert_double_eq(result.matrix[4][2], s21_result.matrix[4][2]);
  ck_assert_double_eq(result.matrix[4][3], s21_result.matrix[4][3]);
  ck_assert_double_eq(result.matrix[4][4], s21_result.matrix[4][4]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_number_13_by_py) {
  int res, s21_res;
  double number = -180.67434068469652;
  matrix_t A, result, s21_result;
  s21_create_matrix(4, 7, &A);
  s21_create_matrix(4, 7, &result);
  A.matrix[0][0] = -607;
  A.matrix[0][1] = 55;
  A.matrix[0][2] = -481;
  A.matrix[0][3] = 605;
  A.matrix[0][4] = 756;
  A.matrix[0][5] = 32;
  A.matrix[0][6] = -140;
  A.matrix[1][0] = 572;
  A.matrix[1][1] = -950;
  A.matrix[1][2] = -285;
  A.matrix[1][3] = -973;
  A.matrix[1][4] = 906;
  A.matrix[1][5] = 330;
  A.matrix[1][6] = 662;
  A.matrix[2][0] = 911;
  A.matrix[2][1] = 47;
  A.matrix[2][2] = 433;
  A.matrix[2][3] = -70;
  A.matrix[2][4] = -680;
  A.matrix[2][5] = -435;
  A.matrix[2][6] = -16;
  A.matrix[3][0] = 825;
  A.matrix[3][1] = 909;
  A.matrix[3][2] = -280;
  A.matrix[3][3] = 788;
  A.matrix[3][4] = -321;
  A.matrix[3][5] = 465;
  A.matrix[3][6] = 387;
  result.matrix[0][0] = 109669.32479561078;
  result.matrix[0][1] = -9937.088737658309;
  result.matrix[0][2] = 86904.35786933903;
  result.matrix[0][3] = -109307.97611424139;
  result.matrix[0][4] = -136589.80155763056;
  result.matrix[0][5] = -5781.578901910289;
  result.matrix[0][6] = 25294.407695857513;
  result.matrix[1][0] = -103345.7228716464;
  result.matrix[1][1] = 171640.6236504617;
  result.matrix[1][2] = 51492.187095138506;
  result.matrix[1][3] = 175796.1334862097;
  result.matrix[1][4] = -163690.95266033505;
  result.matrix[1][5] = -59622.53242594985;
  result.matrix[1][6] = -119606.4135332691;
  result.matrix[2][0] = -164594.32436375853;
  result.matrix[2][1] = -8491.694012180737;
  result.matrix[2][2] = -78231.9895164736;
  result.matrix[2][3] = 12647.203847928757;
  result.matrix[2][4] = 122858.55166559364;
  result.matrix[2][5] = 78593.33819784298;
  result.matrix[2][6] = 2890.7894509551443;
  result.matrix[3][0] = -149056.33106487463;
  result.matrix[3][1] = -164232.97568238914;
  result.matrix[3][2] = 50588.81539171503;
  result.matrix[3][3] = -142371.38045954084;
  result.matrix[3][4] = 57996.46335978758;
  result.matrix[3][5] = -84013.56841838387;
  result.matrix[3][6] = -69920.96984497755;

  res = OK;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);
  ck_assert_double_eq(result.matrix[0][5], s21_result.matrix[0][5]);
  ck_assert_double_eq(result.matrix[0][6], s21_result.matrix[0][6]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[1][4], s21_result.matrix[1][4]);
  ck_assert_double_eq(result.matrix[1][5], s21_result.matrix[1][5]);
  ck_assert_double_eq(result.matrix[1][6], s21_result.matrix[1][6]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[2][4], s21_result.matrix[2][4]);
  ck_assert_double_eq(result.matrix[2][5], s21_result.matrix[2][5]);
  ck_assert_double_eq(result.matrix[2][6], s21_result.matrix[2][6]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[3][4], s21_result.matrix[3][4]);
  ck_assert_double_eq(result.matrix[3][5], s21_result.matrix[3][5]);
  ck_assert_double_eq(result.matrix[3][6], s21_result.matrix[3][6]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(mult_number_14_by_py) {
  int res, s21_res;
  double number = -833.2185862961104;
  matrix_t A, result, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &result);
  A.matrix[0][0] = -544;
  A.matrix[0][1] = -134;
  A.matrix[1][0] = -862;
  A.matrix[1][1] = -674;
  result.matrix[0][0] = 453270.91094508406;
  result.matrix[0][1] = 111651.29056367879;
  result.matrix[1][0] = 718234.4213872472;
  result.matrix[1][1] = 561589.3271635784;

  res = OK;
  s21_res = s21_mult_number(&A, number, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(remove_matrix_1) {
  int rows, columns;
  matrix_t s21_mtrx = {.rows = 0, .columns = 0, .matrix = 0};

  rows = 2;
  columns = 3;

  s21_create_matrix(rows, columns, &s21_mtrx);

  s21_remove_matrix(&s21_mtrx);
  ck_assert_int_eq(OK, s21_mtrx.rows);
  ck_assert_int_eq(OK, s21_mtrx.columns);
}
END_TEST

START_TEST(remove_matrix_2) {
  int rows, columns;
  matrix_t s21_mtrx = {.rows = 0, .columns = 0, .matrix = 0};

  rows = -2;
  columns = 3;

  s21_create_matrix(rows, columns, &s21_mtrx);

  s21_remove_matrix(&s21_mtrx);
  ck_assert_int_eq(OK, s21_mtrx.rows);
  ck_assert_int_eq(OK, s21_mtrx.columns);
}
END_TEST

START_TEST(remove_matrix_3) {
  int rows, columns;
  matrix_t s21_mtrx = {.rows = 0, .columns = 0, .matrix = 0};

  rows = 2;
  columns = -3;

  s21_create_matrix(rows, columns, &s21_mtrx);

  s21_remove_matrix(&s21_mtrx);
  ck_assert_int_eq(OK, s21_mtrx.rows);
  ck_assert_int_eq(OK, s21_mtrx.columns);
}
END_TEST

START_TEST(remove_matrix_4) {
  int rows, columns;
  matrix_t s21_mtrx = {.rows = 0, .columns = 0, .matrix = NULL};

  rows = 0;
  columns = 0;

  s21_create_matrix(rows, columns, &s21_mtrx);

  s21_remove_matrix(&s21_mtrx);
  ck_assert_int_eq(OK, s21_mtrx.rows);
  ck_assert_int_eq(OK, s21_mtrx.columns);
}
END_TEST

START_TEST(sub_matrix_1) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  s21_create_matrix(2, 2, &result);

  A.matrix[0][0] = 2;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 8;
  A.matrix[1][1] = 10;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;
  result.matrix[0][0] = 1;
  result.matrix[0][1] = 2;
  result.matrix[1][0] = 4;
  result.matrix[1][1] = 5;
  result.rows = 2;
  result.columns = 2;

  res = OK;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sub_matrix_2) {
  int res, s21_res;
  matrix_t A, B = {.matrix = 0, .rows = 0, .columns = 0}, s21_result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;

  res = INCORRECT_MATRIX;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(sub_matrix_3) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(3, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sub_matrix_4) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = S21_INFINITY;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sub_matrix_5) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = S21_INFINITY_NEGATIVE;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sub_matrix_6) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 0.0 / 0.0L;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sub_matrix_7) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = -1.4e+308;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 1.4e+308;

  res = CALCULATION_ERROR;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sub_matrix_8) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 0;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 1.7976931348623158e+308;

  res = OK;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sub_matrix_9) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 1.7976931348623158e+308;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 0;

  res = OK;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sub_matrix_10_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(6, 2, &A);
  s21_create_matrix(6, 2, &B);
  s21_create_matrix(6, 2, &result);
  A.matrix[0][0] = 143;
  A.matrix[0][1] = -749;
  A.matrix[1][0] = 35;
  A.matrix[1][1] = 854;
  A.matrix[2][0] = 492;
  A.matrix[2][1] = 673;
  A.matrix[3][0] = 216;
  A.matrix[3][1] = 508;
  A.matrix[4][0] = -410;
  A.matrix[4][1] = 161;
  A.matrix[5][0] = -46;
  A.matrix[5][1] = -18;
  B.matrix[0][0] = -176;
  B.matrix[0][1] = 511;
  B.matrix[1][0] = 505;
  B.matrix[1][1] = -379;
  B.matrix[2][0] = 589;
  B.matrix[2][1] = 90;
  B.matrix[3][0] = 554;
  B.matrix[3][1] = 772;
  B.matrix[4][0] = -997;
  B.matrix[4][1] = 195;
  B.matrix[5][0] = 125;
  B.matrix[5][1] = 947;
  result.matrix[0][0] = 319;
  result.matrix[0][1] = -1260;
  result.matrix[1][0] = -470;
  result.matrix[1][1] = 1233;
  result.matrix[2][0] = -97;
  result.matrix[2][1] = 583;
  result.matrix[3][0] = -338;
  result.matrix[3][1] = -264;
  result.matrix[4][0] = 587;
  result.matrix[4][1] = -34;
  result.matrix[5][0] = -171;
  result.matrix[5][1] = -965;

  res = OK;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[4][0], s21_result.matrix[4][0]);
  ck_assert_double_eq(result.matrix[4][1], s21_result.matrix[4][1]);
  ck_assert_double_eq(result.matrix[5][0], s21_result.matrix[5][0]);
  ck_assert_double_eq(result.matrix[5][1], s21_result.matrix[5][1]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sub_matrix_11_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(10, 10, &A);
  s21_create_matrix(10, 10, &B);
  s21_create_matrix(10, 10, &result);
  A.matrix[0][0] = -152;
  A.matrix[0][1] = 807;
  A.matrix[0][2] = 649;
  A.matrix[0][3] = -687;
  A.matrix[0][4] = -190;
  A.matrix[0][5] = 417;
  A.matrix[0][6] = 806;
  A.matrix[0][7] = -722;
  A.matrix[0][8] = -70;
  A.matrix[0][9] = -33;
  A.matrix[1][0] = -701;
  A.matrix[1][1] = -774;
  A.matrix[1][2] = -801;
  A.matrix[1][3] = -184;
  A.matrix[1][4] = 915;
  A.matrix[1][5] = -407;
  A.matrix[1][6] = -429;
  A.matrix[1][7] = 368;
  A.matrix[1][8] = -696;
  A.matrix[1][9] = 759;
  A.matrix[2][0] = 874;
  A.matrix[2][1] = 381;
  A.matrix[2][2] = -684;
  A.matrix[2][3] = 429;
  A.matrix[2][4] = 715;
  A.matrix[2][5] = -883;
  A.matrix[2][6] = 949;
  A.matrix[2][7] = -826;
  A.matrix[2][8] = -997;
  A.matrix[2][9] = 224;
  A.matrix[3][0] = -796;
  A.matrix[3][1] = -551;
  A.matrix[3][2] = 90;
  A.matrix[3][3] = -225;
  A.matrix[3][4] = -392;
  A.matrix[3][5] = -746;
  A.matrix[3][6] = 142;
  A.matrix[3][7] = 436;
  A.matrix[3][8] = -938;
  A.matrix[3][9] = -659;
  A.matrix[4][0] = -843;
  A.matrix[4][1] = 704;
  A.matrix[4][2] = 82;
  A.matrix[4][3] = -236;
  A.matrix[4][4] = -283;
  A.matrix[4][5] = 211;
  A.matrix[4][6] = 125;
  A.matrix[4][7] = 770;
  A.matrix[4][8] = -324;
  A.matrix[4][9] = 410;
  A.matrix[5][0] = 442;
  A.matrix[5][1] = 505;
  A.matrix[5][2] = 652;
  A.matrix[5][3] = -830;
  A.matrix[5][4] = -889;
  A.matrix[5][5] = 466;
  A.matrix[5][6] = -79;
  A.matrix[5][7] = -928;
  A.matrix[5][8] = 625;
  A.matrix[5][9] = 835;
  A.matrix[6][0] = 247;
  A.matrix[6][1] = -621;
  A.matrix[6][2] = -898;
  A.matrix[6][3] = -868;
  A.matrix[6][4] = 990;
  A.matrix[6][5] = -948;
  A.matrix[6][6] = -982;
  A.matrix[6][7] = 507;
  A.matrix[6][8] = 51;
  A.matrix[6][9] = 848;
  A.matrix[7][0] = 787;
  A.matrix[7][1] = -373;
  A.matrix[7][2] = -910;
  A.matrix[7][3] = 425;
  A.matrix[7][4] = 417;
  A.matrix[7][5] = -17;
  A.matrix[7][6] = -707;
  A.matrix[7][7] = -213;
  A.matrix[7][8] = -100;
  A.matrix[7][9] = -96;
  A.matrix[8][0] = -714;
  A.matrix[8][1] = -435;
  A.matrix[8][2] = -376;
  A.matrix[8][3] = 766;
  A.matrix[8][4] = 682;
  A.matrix[8][5] = -98;
  A.matrix[8][6] = -249;
  A.matrix[8][7] = -994;
  A.matrix[8][8] = 471;
  A.matrix[8][9] = -947;
  A.matrix[9][0] = 993;
  A.matrix[9][1] = 350;
  A.matrix[9][2] = -594;
  A.matrix[9][3] = 467;
  A.matrix[9][4] = -154;
  A.matrix[9][5] = -468;
  A.matrix[9][6] = -810;
  A.matrix[9][7] = 31;
  A.matrix[9][8] = 622;
  A.matrix[9][9] = -342;
  B.matrix[0][0] = 7;
  B.matrix[0][1] = -252;
  B.matrix[0][2] = -615;
  B.matrix[0][3] = -141;
  B.matrix[0][4] = 189;
  B.matrix[0][5] = 849;
  B.matrix[0][6] = -94;
  B.matrix[0][7] = 993;
  B.matrix[0][8] = 620;
  B.matrix[0][9] = 467;
  B.matrix[1][0] = -261;
  B.matrix[1][1] = 589;
  B.matrix[1][2] = -4;
  B.matrix[1][3] = -890;
  B.matrix[1][4] = 686;
  B.matrix[1][5] = -914;
  B.matrix[1][6] = -727;
  B.matrix[1][7] = 444;
  B.matrix[1][8] = -58;
  B.matrix[1][9] = 467;
  B.matrix[2][0] = 361;
  B.matrix[2][1] = -156;
  B.matrix[2][2] = 794;
  B.matrix[2][3] = -762;
  B.matrix[2][4] = -603;
  B.matrix[2][5] = -992;
  B.matrix[2][6] = -66;
  B.matrix[2][7] = 27;
  B.matrix[2][8] = -28;
  B.matrix[2][9] = 379;
  B.matrix[3][0] = 102;
  B.matrix[3][1] = -438;
  B.matrix[3][2] = 782;
  B.matrix[3][3] = -811;
  B.matrix[3][4] = -982;
  B.matrix[3][5] = -619;
  B.matrix[3][6] = -653;
  B.matrix[3][7] = 163;
  B.matrix[3][8] = -769;
  B.matrix[3][9] = -325;
  B.matrix[4][0] = 855;
  B.matrix[4][1] = -550;
  B.matrix[4][2] = -320;
  B.matrix[4][3] = 534;
  B.matrix[4][4] = -83;
  B.matrix[4][5] = -804;
  B.matrix[4][6] = 820;
  B.matrix[4][7] = 235;
  B.matrix[4][8] = -799;
  B.matrix[4][9] = -300;
  B.matrix[5][0] = -836;
  B.matrix[5][1] = -931;
  B.matrix[5][2] = 551;
  B.matrix[5][3] = 602;
  B.matrix[5][4] = 506;
  B.matrix[5][5] = 115;
  B.matrix[5][6] = -499;
  B.matrix[5][7] = -514;
  B.matrix[5][8] = 451;
  B.matrix[5][9] = 646;
  B.matrix[6][0] = 880;
  B.matrix[6][1] = -896;
  B.matrix[6][2] = -538;
  B.matrix[6][3] = -205;
  B.matrix[6][4] = -66;
  B.matrix[6][5] = 757;
  B.matrix[6][6] = 917;
  B.matrix[6][7] = -65;
  B.matrix[6][8] = -437;
  B.matrix[6][9] = -588;
  B.matrix[7][0] = 261;
  B.matrix[7][1] = -86;
  B.matrix[7][2] = -829;
  B.matrix[7][3] = 808;
  B.matrix[7][4] = 440;
  B.matrix[7][5] = -362;
  B.matrix[7][6] = -350;
  B.matrix[7][7] = 596;
  B.matrix[7][8] = -813;
  B.matrix[7][9] = 422;
  B.matrix[8][0] = 12;
  B.matrix[8][1] = -617;
  B.matrix[8][2] = 88;
  B.matrix[8][3] = 221;
  B.matrix[8][4] = 346;
  B.matrix[8][5] = -358;
  B.matrix[8][6] = -783;
  B.matrix[8][7] = -539;
  B.matrix[8][8] = -171;
  B.matrix[8][9] = 43;
  B.matrix[9][0] = 362;
  B.matrix[9][1] = -718;
  B.matrix[9][2] = -875;
  B.matrix[9][3] = 620;
  B.matrix[9][4] = -805;
  B.matrix[9][5] = -628;
  B.matrix[9][6] = -223;
  B.matrix[9][7] = 211;
  B.matrix[9][8] = 750;
  B.matrix[9][9] = -861;
  result.matrix[0][0] = -159;
  result.matrix[0][1] = 1059;
  result.matrix[0][2] = 1264;
  result.matrix[0][3] = -546;
  result.matrix[0][4] = -379;
  result.matrix[0][5] = -432;
  result.matrix[0][6] = 900;
  result.matrix[0][7] = -1715;
  result.matrix[0][8] = -690;
  result.matrix[0][9] = -500;
  result.matrix[1][0] = -440;
  result.matrix[1][1] = -1363;
  result.matrix[1][2] = -797;
  result.matrix[1][3] = 706;
  result.matrix[1][4] = 229;
  result.matrix[1][5] = 507;
  result.matrix[1][6] = 298;
  result.matrix[1][7] = -76;
  result.matrix[1][8] = -638;
  result.matrix[1][9] = 292;
  result.matrix[2][0] = 513;
  result.matrix[2][1] = 537;
  result.matrix[2][2] = -1478;
  result.matrix[2][3] = 1191;
  result.matrix[2][4] = 1318;
  result.matrix[2][5] = 109;
  result.matrix[2][6] = 1015;
  result.matrix[2][7] = -853;
  result.matrix[2][8] = -969;
  result.matrix[2][9] = -155;
  result.matrix[3][0] = -898;
  result.matrix[3][1] = -113;
  result.matrix[3][2] = -692;
  result.matrix[3][3] = 586;
  result.matrix[3][4] = 590;
  result.matrix[3][5] = -127;
  result.matrix[3][6] = 795;
  result.matrix[3][7] = 273;
  result.matrix[3][8] = -169;
  result.matrix[3][9] = -334;
  result.matrix[4][0] = -1698;
  result.matrix[4][1] = 1254;
  result.matrix[4][2] = 402;
  result.matrix[4][3] = -770;
  result.matrix[4][4] = -200;
  result.matrix[4][5] = 1015;
  result.matrix[4][6] = -695;
  result.matrix[4][7] = 535;
  result.matrix[4][8] = 475;
  result.matrix[4][9] = 710;
  result.matrix[5][0] = 1278;
  result.matrix[5][1] = 1436;
  result.matrix[5][2] = 101;
  result.matrix[5][3] = -1432;
  result.matrix[5][4] = -1395;
  result.matrix[5][5] = 351;
  result.matrix[5][6] = 420;
  result.matrix[5][7] = -414;
  result.matrix[5][8] = 174;
  result.matrix[5][9] = 189;
  result.matrix[6][0] = -633;
  result.matrix[6][1] = 275;
  result.matrix[6][2] = -360;
  result.matrix[6][3] = -663;
  result.matrix[6][4] = 1056;
  result.matrix[6][5] = -1705;
  result.matrix[6][6] = -1899;
  result.matrix[6][7] = 572;
  result.matrix[6][8] = 488;
  result.matrix[6][9] = 1436;
  result.matrix[7][0] = 526;
  result.matrix[7][1] = -287;
  result.matrix[7][2] = -81;
  result.matrix[7][3] = -383;
  result.matrix[7][4] = -23;
  result.matrix[7][5] = 345;
  result.matrix[7][6] = -357;
  result.matrix[7][7] = -809;
  result.matrix[7][8] = 713;
  result.matrix[7][9] = -518;
  result.matrix[8][0] = -726;
  result.matrix[8][1] = 182;
  result.matrix[8][2] = -464;
  result.matrix[8][3] = 545;
  result.matrix[8][4] = 336;
  result.matrix[8][5] = 260;
  result.matrix[8][6] = 534;
  result.matrix[8][7] = -455;
  result.matrix[8][8] = 642;
  result.matrix[8][9] = -990;
  result.matrix[9][0] = 631;
  result.matrix[9][1] = 1068;
  result.matrix[9][2] = 281;
  result.matrix[9][3] = -153;
  result.matrix[9][4] = 651;
  result.matrix[9][5] = 160;
  result.matrix[9][6] = -587;
  result.matrix[9][7] = -180;
  result.matrix[9][8] = -128;
  result.matrix[9][9] = 519;

  res = OK;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);
  ck_assert_double_eq(result.matrix[0][5], s21_result.matrix[0][5]);
  ck_assert_double_eq(result.matrix[0][6], s21_result.matrix[0][6]);
  ck_assert_double_eq(result.matrix[0][7], s21_result.matrix[0][7]);
  ck_assert_double_eq(result.matrix[0][8], s21_result.matrix[0][8]);
  ck_assert_double_eq(result.matrix[0][9], s21_result.matrix[0][9]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[1][4], s21_result.matrix[1][4]);
  ck_assert_double_eq(result.matrix[1][5], s21_result.matrix[1][5]);
  ck_assert_double_eq(result.matrix[1][6], s21_result.matrix[1][6]);
  ck_assert_double_eq(result.matrix[1][7], s21_result.matrix[1][7]);
  ck_assert_double_eq(result.matrix[1][8], s21_result.matrix[1][8]);
  ck_assert_double_eq(result.matrix[1][9], s21_result.matrix[1][9]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[2][4], s21_result.matrix[2][4]);
  ck_assert_double_eq(result.matrix[2][5], s21_result.matrix[2][5]);
  ck_assert_double_eq(result.matrix[2][6], s21_result.matrix[2][6]);
  ck_assert_double_eq(result.matrix[2][7], s21_result.matrix[2][7]);
  ck_assert_double_eq(result.matrix[2][8], s21_result.matrix[2][8]);
  ck_assert_double_eq(result.matrix[2][9], s21_result.matrix[2][9]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[3][4], s21_result.matrix[3][4]);
  ck_assert_double_eq(result.matrix[3][5], s21_result.matrix[3][5]);
  ck_assert_double_eq(result.matrix[3][6], s21_result.matrix[3][6]);
  ck_assert_double_eq(result.matrix[3][7], s21_result.matrix[3][7]);
  ck_assert_double_eq(result.matrix[3][8], s21_result.matrix[3][8]);
  ck_assert_double_eq(result.matrix[3][9], s21_result.matrix[3][9]);
  ck_assert_double_eq(result.matrix[4][0], s21_result.matrix[4][0]);
  ck_assert_double_eq(result.matrix[4][1], s21_result.matrix[4][1]);
  ck_assert_double_eq(result.matrix[4][2], s21_result.matrix[4][2]);
  ck_assert_double_eq(result.matrix[4][3], s21_result.matrix[4][3]);
  ck_assert_double_eq(result.matrix[4][4], s21_result.matrix[4][4]);
  ck_assert_double_eq(result.matrix[4][5], s21_result.matrix[4][5]);
  ck_assert_double_eq(result.matrix[4][6], s21_result.matrix[4][6]);
  ck_assert_double_eq(result.matrix[4][7], s21_result.matrix[4][7]);
  ck_assert_double_eq(result.matrix[4][8], s21_result.matrix[4][8]);
  ck_assert_double_eq(result.matrix[4][9], s21_result.matrix[4][9]);
  ck_assert_double_eq(result.matrix[5][0], s21_result.matrix[5][0]);
  ck_assert_double_eq(result.matrix[5][1], s21_result.matrix[5][1]);
  ck_assert_double_eq(result.matrix[5][2], s21_result.matrix[5][2]);
  ck_assert_double_eq(result.matrix[5][3], s21_result.matrix[5][3]);
  ck_assert_double_eq(result.matrix[5][4], s21_result.matrix[5][4]);
  ck_assert_double_eq(result.matrix[5][5], s21_result.matrix[5][5]);
  ck_assert_double_eq(result.matrix[5][6], s21_result.matrix[5][6]);
  ck_assert_double_eq(result.matrix[5][7], s21_result.matrix[5][7]);
  ck_assert_double_eq(result.matrix[5][8], s21_result.matrix[5][8]);
  ck_assert_double_eq(result.matrix[5][9], s21_result.matrix[5][9]);
  ck_assert_double_eq(result.matrix[6][0], s21_result.matrix[6][0]);
  ck_assert_double_eq(result.matrix[6][1], s21_result.matrix[6][1]);
  ck_assert_double_eq(result.matrix[6][2], s21_result.matrix[6][2]);
  ck_assert_double_eq(result.matrix[6][3], s21_result.matrix[6][3]);
  ck_assert_double_eq(result.matrix[6][4], s21_result.matrix[6][4]);
  ck_assert_double_eq(result.matrix[6][5], s21_result.matrix[6][5]);
  ck_assert_double_eq(result.matrix[6][6], s21_result.matrix[6][6]);
  ck_assert_double_eq(result.matrix[6][7], s21_result.matrix[6][7]);
  ck_assert_double_eq(result.matrix[6][8], s21_result.matrix[6][8]);
  ck_assert_double_eq(result.matrix[6][9], s21_result.matrix[6][9]);
  ck_assert_double_eq(result.matrix[7][0], s21_result.matrix[7][0]);
  ck_assert_double_eq(result.matrix[7][1], s21_result.matrix[7][1]);
  ck_assert_double_eq(result.matrix[7][2], s21_result.matrix[7][2]);
  ck_assert_double_eq(result.matrix[7][3], s21_result.matrix[7][3]);
  ck_assert_double_eq(result.matrix[7][4], s21_result.matrix[7][4]);
  ck_assert_double_eq(result.matrix[7][5], s21_result.matrix[7][5]);
  ck_assert_double_eq(result.matrix[7][6], s21_result.matrix[7][6]);
  ck_assert_double_eq(result.matrix[7][7], s21_result.matrix[7][7]);
  ck_assert_double_eq(result.matrix[7][8], s21_result.matrix[7][8]);
  ck_assert_double_eq(result.matrix[7][9], s21_result.matrix[7][9]);
  ck_assert_double_eq(result.matrix[8][0], s21_result.matrix[8][0]);
  ck_assert_double_eq(result.matrix[8][1], s21_result.matrix[8][1]);
  ck_assert_double_eq(result.matrix[8][2], s21_result.matrix[8][2]);
  ck_assert_double_eq(result.matrix[8][3], s21_result.matrix[8][3]);
  ck_assert_double_eq(result.matrix[8][4], s21_result.matrix[8][4]);
  ck_assert_double_eq(result.matrix[8][5], s21_result.matrix[8][5]);
  ck_assert_double_eq(result.matrix[8][6], s21_result.matrix[8][6]);
  ck_assert_double_eq(result.matrix[8][7], s21_result.matrix[8][7]);
  ck_assert_double_eq(result.matrix[8][8], s21_result.matrix[8][8]);
  ck_assert_double_eq(result.matrix[8][9], s21_result.matrix[8][9]);
  ck_assert_double_eq(result.matrix[9][0], s21_result.matrix[9][0]);
  ck_assert_double_eq(result.matrix[9][1], s21_result.matrix[9][1]);
  ck_assert_double_eq(result.matrix[9][2], s21_result.matrix[9][2]);
  ck_assert_double_eq(result.matrix[9][3], s21_result.matrix[9][3]);
  ck_assert_double_eq(result.matrix[9][4], s21_result.matrix[9][4]);
  ck_assert_double_eq(result.matrix[9][5], s21_result.matrix[9][5]);
  ck_assert_double_eq(result.matrix[9][6], s21_result.matrix[9][6]);
  ck_assert_double_eq(result.matrix[9][7], s21_result.matrix[9][7]);
  ck_assert_double_eq(result.matrix[9][8], s21_result.matrix[9][8]);
  ck_assert_double_eq(result.matrix[9][9], s21_result.matrix[9][9]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sub_matrix_12_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(4, 1, &A);
  s21_create_matrix(4, 1, &B);
  s21_create_matrix(4, 1, &result);
  A.matrix[0][0] = -633;
  A.matrix[1][0] = -658;
  A.matrix[2][0] = 120;
  A.matrix[3][0] = 422;
  B.matrix[0][0] = 259;
  B.matrix[1][0] = -921;
  B.matrix[2][0] = 36;
  B.matrix[3][0] = 84;
  result.matrix[0][0] = -892;
  result.matrix[1][0] = 263;
  result.matrix[2][0] = 84;
  result.matrix[3][0] = 338;

  res = OK;
  s21_res = s21_sub_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_1) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  s21_create_matrix(2, 2, &result);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;
  result.matrix[0][0] = 2;
  result.matrix[0][1] = 4;
  result.matrix[1][0] = 8;
  result.matrix[1][1] = 10;
  result.rows = 2;
  result.columns = 2;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_2) {
  int res, s21_res;
  matrix_t A, B = {.matrix = 0, .rows = 0, .columns = 0}, s21_result;
  s21_create_matrix(2, 2, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;

  res = INCORRECT_MATRIX;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(sum_matrix_3) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(3, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sum_matrix_4) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = S21_INFINITY;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sum_matrix_5) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = S21_INFINITY_NEGATIVE;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sum_matrix_6) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 0.0 / 0.0L;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 5;

  res = CALCULATION_ERROR;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sum_matrix_7) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 1.4e+308;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 1.4e+308;

  res = CALCULATION_ERROR;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sum_matrix_8) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 1.7976931348623158e+308;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 0;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_9) {
  int res, s21_res;
  matrix_t A, B, s21_result;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = -1.7976931348623158e+308;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = 0;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_10_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(6, 6, &A);
  s21_create_matrix(6, 6, &B);
  s21_create_matrix(6, 6, &result);
  A.matrix[0][0] = -606;
  A.matrix[0][1] = 642;
  A.matrix[0][2] = 904;
  A.matrix[0][3] = 134;
  A.matrix[0][4] = 624;
  A.matrix[0][5] = -900;
  A.matrix[1][0] = -210;
  A.matrix[1][1] = -616;
  A.matrix[1][2] = -490;
  A.matrix[1][3] = -877;
  A.matrix[1][4] = -461;
  A.matrix[1][5] = 207;
  A.matrix[2][0] = 19;
  A.matrix[2][1] = 523;
  A.matrix[2][2] = -723;
  A.matrix[2][3] = -583;
  A.matrix[2][4] = -367;
  A.matrix[2][5] = -351;
  A.matrix[3][0] = 574;
  A.matrix[3][1] = 162;
  A.matrix[3][2] = 975;
  A.matrix[3][3] = 4;
  A.matrix[3][4] = 49;
  A.matrix[3][5] = 417;
  A.matrix[4][0] = -286;
  A.matrix[4][1] = 151;
  A.matrix[4][2] = -876;
  A.matrix[4][3] = 689;
  A.matrix[4][4] = 254;
  A.matrix[4][5] = 535;
  A.matrix[5][0] = 191;
  A.matrix[5][1] = -2;
  A.matrix[5][2] = -30;
  A.matrix[5][3] = 26;
  A.matrix[5][4] = 234;
  A.matrix[5][5] = -246;
  B.matrix[0][0] = -971;
  B.matrix[0][1] = 204;
  B.matrix[0][2] = 564;
  B.matrix[0][3] = -860;
  B.matrix[0][4] = -636;
  B.matrix[0][5] = -432;
  B.matrix[1][0] = -40;
  B.matrix[1][1] = -576;
  B.matrix[1][2] = -406;
  B.matrix[1][3] = 416;
  B.matrix[1][4] = -900;
  B.matrix[1][5] = 631;
  B.matrix[2][0] = -424;
  B.matrix[2][1] = 26;
  B.matrix[2][2] = -84;
  B.matrix[2][3] = -798;
  B.matrix[2][4] = 795;
  B.matrix[2][5] = -533;
  B.matrix[3][0] = 479;
  B.matrix[3][1] = -695;
  B.matrix[3][2] = -43;
  B.matrix[3][3] = 62;
  B.matrix[3][4] = 284;
  B.matrix[3][5] = 431;
  B.matrix[4][0] = 819;
  B.matrix[4][1] = -498;
  B.matrix[4][2] = 602;
  B.matrix[4][3] = -61;
  B.matrix[4][4] = -359;
  B.matrix[4][5] = 146;
  B.matrix[5][0] = -316;
  B.matrix[5][1] = 196;
  B.matrix[5][2] = -63;
  B.matrix[5][3] = -238;
  B.matrix[5][4] = 615;
  B.matrix[5][5] = -964;
  result.matrix[0][0] = -1577;
  result.matrix[0][1] = 846;
  result.matrix[0][2] = 1468;
  result.matrix[0][3] = -726;
  result.matrix[0][4] = -12;
  result.matrix[0][5] = -1332;
  result.matrix[1][0] = -250;
  result.matrix[1][1] = -1192;
  result.matrix[1][2] = -896;
  result.matrix[1][3] = -461;
  result.matrix[1][4] = -1361;
  result.matrix[1][5] = 838;
  result.matrix[2][0] = -405;
  result.matrix[2][1] = 549;
  result.matrix[2][2] = -807;
  result.matrix[2][3] = -1381;
  result.matrix[2][4] = 428;
  result.matrix[2][5] = -884;
  result.matrix[3][0] = 1053;
  result.matrix[3][1] = -533;
  result.matrix[3][2] = 932;
  result.matrix[3][3] = 66;
  result.matrix[3][4] = 333;
  result.matrix[3][5] = 848;
  result.matrix[4][0] = 533;
  result.matrix[4][1] = -347;
  result.matrix[4][2] = -274;
  result.matrix[4][3] = 628;
  result.matrix[4][4] = -105;
  result.matrix[4][5] = 681;
  result.matrix[5][0] = -125;
  result.matrix[5][1] = 194;
  result.matrix[5][2] = -93;
  result.matrix[5][3] = -212;
  result.matrix[5][4] = 849;
  result.matrix[5][5] = -1210;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);
  ck_assert_double_eq(result.matrix[0][5], s21_result.matrix[0][5]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[1][4], s21_result.matrix[1][4]);
  ck_assert_double_eq(result.matrix[1][5], s21_result.matrix[1][5]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[2][4], s21_result.matrix[2][4]);
  ck_assert_double_eq(result.matrix[2][5], s21_result.matrix[2][5]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[3][4], s21_result.matrix[3][4]);
  ck_assert_double_eq(result.matrix[3][5], s21_result.matrix[3][5]);
  ck_assert_double_eq(result.matrix[4][0], s21_result.matrix[4][0]);
  ck_assert_double_eq(result.matrix[4][1], s21_result.matrix[4][1]);
  ck_assert_double_eq(result.matrix[4][2], s21_result.matrix[4][2]);
  ck_assert_double_eq(result.matrix[4][3], s21_result.matrix[4][3]);
  ck_assert_double_eq(result.matrix[4][4], s21_result.matrix[4][4]);
  ck_assert_double_eq(result.matrix[4][5], s21_result.matrix[4][5]);
  ck_assert_double_eq(result.matrix[5][0], s21_result.matrix[5][0]);
  ck_assert_double_eq(result.matrix[5][1], s21_result.matrix[5][1]);
  ck_assert_double_eq(result.matrix[5][2], s21_result.matrix[5][2]);
  ck_assert_double_eq(result.matrix[5][3], s21_result.matrix[5][3]);
  ck_assert_double_eq(result.matrix[5][4], s21_result.matrix[5][4]);
  ck_assert_double_eq(result.matrix[5][5], s21_result.matrix[5][5]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_11_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(8, 4, &A);
  s21_create_matrix(8, 4, &B);
  s21_create_matrix(8, 4, &result);
  A.matrix[0][0] = 107;
  A.matrix[0][1] = -956;
  A.matrix[0][2] = -824;
  A.matrix[0][3] = 65;
  A.matrix[1][0] = 270;
  A.matrix[1][1] = -978;
  A.matrix[1][2] = -235;
  A.matrix[1][3] = 613;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = -695;
  A.matrix[2][2] = 164;
  A.matrix[2][3] = 584;
  A.matrix[3][0] = -354;
  A.matrix[3][1] = 130;
  A.matrix[3][2] = -928;
  A.matrix[3][3] = -632;
  A.matrix[4][0] = 359;
  A.matrix[4][1] = -34;
  A.matrix[4][2] = -20;
  A.matrix[4][3] = -772;
  A.matrix[5][0] = -176;
  A.matrix[5][1] = 257;
  A.matrix[5][2] = -487;
  A.matrix[5][3] = -744;
  A.matrix[6][0] = 710;
  A.matrix[6][1] = -371;
  A.matrix[6][2] = 893;
  A.matrix[6][3] = 212;
  A.matrix[7][0] = 922;
  A.matrix[7][1] = 834;
  A.matrix[7][2] = 578;
  A.matrix[7][3] = 903;
  B.matrix[0][0] = -408;
  B.matrix[0][1] = -317;
  B.matrix[0][2] = 417;
  B.matrix[0][3] = 967;
  B.matrix[1][0] = -3;
  B.matrix[1][1] = 375;
  B.matrix[1][2] = -680;
  B.matrix[1][3] = 856;
  B.matrix[2][0] = 445;
  B.matrix[2][1] = -72;
  B.matrix[2][2] = 212;
  B.matrix[2][3] = -894;
  B.matrix[3][0] = -821;
  B.matrix[3][1] = -942;
  B.matrix[3][2] = 880;
  B.matrix[3][3] = 36;
  B.matrix[4][0] = 95;
  B.matrix[4][1] = -365;
  B.matrix[4][2] = 817;
  B.matrix[4][3] = 291;
  B.matrix[5][0] = -765;
  B.matrix[5][1] = -183;
  B.matrix[5][2] = -72;
  B.matrix[5][3] = 699;
  B.matrix[6][0] = -757;
  B.matrix[6][1] = 136;
  B.matrix[6][2] = 167;
  B.matrix[6][3] = -741;
  B.matrix[7][0] = 147;
  B.matrix[7][1] = 777;
  B.matrix[7][2] = -341;
  B.matrix[7][3] = -876;
  result.matrix[0][0] = -301;
  result.matrix[0][1] = -1273;
  result.matrix[0][2] = -407;
  result.matrix[0][3] = 1032;
  result.matrix[1][0] = 267;
  result.matrix[1][1] = -603;
  result.matrix[1][2] = -915;
  result.matrix[1][3] = 1469;
  result.matrix[2][0] = 447;
  result.matrix[2][1] = -767;
  result.matrix[2][2] = 376;
  result.matrix[2][3] = -310;
  result.matrix[3][0] = -1175;
  result.matrix[3][1] = -812;
  result.matrix[3][2] = -48;
  result.matrix[3][3] = -596;
  result.matrix[4][0] = 454;
  result.matrix[4][1] = -399;
  result.matrix[4][2] = 797;
  result.matrix[4][3] = -481;
  result.matrix[5][0] = -941;
  result.matrix[5][1] = 74;
  result.matrix[5][2] = -559;
  result.matrix[5][3] = -45;
  result.matrix[6][0] = -47;
  result.matrix[6][1] = -235;
  result.matrix[6][2] = 1060;
  result.matrix[6][3] = -529;
  result.matrix[7][0] = 1069;
  result.matrix[7][1] = 1611;
  result.matrix[7][2] = 237;
  result.matrix[7][3] = 27;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[4][0], s21_result.matrix[4][0]);
  ck_assert_double_eq(result.matrix[4][1], s21_result.matrix[4][1]);
  ck_assert_double_eq(result.matrix[4][2], s21_result.matrix[4][2]);
  ck_assert_double_eq(result.matrix[4][3], s21_result.matrix[4][3]);
  ck_assert_double_eq(result.matrix[5][0], s21_result.matrix[5][0]);
  ck_assert_double_eq(result.matrix[5][1], s21_result.matrix[5][1]);
  ck_assert_double_eq(result.matrix[5][2], s21_result.matrix[5][2]);
  ck_assert_double_eq(result.matrix[5][3], s21_result.matrix[5][3]);
  ck_assert_double_eq(result.matrix[6][0], s21_result.matrix[6][0]);
  ck_assert_double_eq(result.matrix[6][1], s21_result.matrix[6][1]);
  ck_assert_double_eq(result.matrix[6][2], s21_result.matrix[6][2]);
  ck_assert_double_eq(result.matrix[6][3], s21_result.matrix[6][3]);
  ck_assert_double_eq(result.matrix[7][0], s21_result.matrix[7][0]);
  ck_assert_double_eq(result.matrix[7][1], s21_result.matrix[7][1]);
  ck_assert_double_eq(result.matrix[7][2], s21_result.matrix[7][2]);
  ck_assert_double_eq(result.matrix[7][3], s21_result.matrix[7][3]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_12_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(4, 12, &A);
  s21_create_matrix(4, 12, &B);
  s21_create_matrix(4, 12, &result);
  A.matrix[0][0] = -523;
  A.matrix[0][1] = -232;
  A.matrix[0][2] = -152;
  A.matrix[0][3] = -925;
  A.matrix[0][4] = -747;
  A.matrix[0][5] = -261;
  A.matrix[0][6] = 583;
  A.matrix[0][7] = -525;
  A.matrix[0][8] = -873;
  A.matrix[0][9] = 630;
  A.matrix[0][10] = 234;
  A.matrix[0][11] = -241;
  A.matrix[1][0] = 566;
  A.matrix[1][1] = 699;
  A.matrix[1][2] = -827;
  A.matrix[1][3] = 14;
  A.matrix[1][4] = -970;
  A.matrix[1][5] = 37;
  A.matrix[1][6] = -623;
  A.matrix[1][7] = -728;
  A.matrix[1][8] = 295;
  A.matrix[1][9] = 225;
  A.matrix[1][10] = -446;
  A.matrix[1][11] = -797;
  A.matrix[2][0] = 978;
  A.matrix[2][1] = -652;
  A.matrix[2][2] = -207;
  A.matrix[2][3] = -175;
  A.matrix[2][4] = -733;
  A.matrix[2][5] = -794;
  A.matrix[2][6] = -802;
  A.matrix[2][7] = -962;
  A.matrix[2][8] = -302;
  A.matrix[2][9] = -71;
  A.matrix[2][10] = -38;
  A.matrix[2][11] = 98;
  A.matrix[3][0] = 248;
  A.matrix[3][1] = 369;
  A.matrix[3][2] = 325;
  A.matrix[3][3] = -644;
  A.matrix[3][4] = 605;
  A.matrix[3][5] = 349;
  A.matrix[3][6] = -268;
  A.matrix[3][7] = -542;
  A.matrix[3][8] = -467;
  A.matrix[3][9] = -202;
  A.matrix[3][10] = 528;
  A.matrix[3][11] = -211;
  B.matrix[0][0] = -455;
  B.matrix[0][1] = 546;
  B.matrix[0][2] = 832;
  B.matrix[0][3] = -263;
  B.matrix[0][4] = -158;
  B.matrix[0][5] = 697;
  B.matrix[0][6] = -425;
  B.matrix[0][7] = 696;
  B.matrix[0][8] = -685;
  B.matrix[0][9] = 731;
  B.matrix[0][10] = -215;
  B.matrix[0][11] = 554;
  B.matrix[1][0] = 408;
  B.matrix[1][1] = -601;
  B.matrix[1][2] = -333;
  B.matrix[1][3] = -948;
  B.matrix[1][4] = -421;
  B.matrix[1][5] = -507;
  B.matrix[1][6] = 964;
  B.matrix[1][7] = 875;
  B.matrix[1][8] = -134;
  B.matrix[1][9] = 576;
  B.matrix[1][10] = 704;
  B.matrix[1][11] = -97;
  B.matrix[2][0] = 881;
  B.matrix[2][1] = 686;
  B.matrix[2][2] = -937;
  B.matrix[2][3] = 580;
  B.matrix[2][4] = -746;
  B.matrix[2][5] = -402;
  B.matrix[2][6] = -235;
  B.matrix[2][7] = 53;
  B.matrix[2][8] = 707;
  B.matrix[2][9] = 342;
  B.matrix[2][10] = 760;
  B.matrix[2][11] = -35;
  B.matrix[3][0] = 194;
  B.matrix[3][1] = -389;
  B.matrix[3][2] = -179;
  B.matrix[3][3] = -954;
  B.matrix[3][4] = 47;
  B.matrix[3][5] = -736;
  B.matrix[3][6] = 699;
  B.matrix[3][7] = 574;
  B.matrix[3][8] = 1;
  B.matrix[3][9] = 42;
  B.matrix[3][10] = 38;
  B.matrix[3][11] = -289;
  result.matrix[0][0] = -978;
  result.matrix[0][1] = 314;
  result.matrix[0][2] = 680;
  result.matrix[0][3] = -1188;
  result.matrix[0][4] = -905;
  result.matrix[0][5] = 436;
  result.matrix[0][6] = 158;
  result.matrix[0][7] = 171;
  result.matrix[0][8] = -1558;
  result.matrix[0][9] = 1361;
  result.matrix[0][10] = 19;
  result.matrix[0][11] = 313;
  result.matrix[1][0] = 974;
  result.matrix[1][1] = 98;
  result.matrix[1][2] = -1160;
  result.matrix[1][3] = -934;
  result.matrix[1][4] = -1391;
  result.matrix[1][5] = -470;
  result.matrix[1][6] = 341;
  result.matrix[1][7] = 147;
  result.matrix[1][8] = 161;
  result.matrix[1][9] = 801;
  result.matrix[1][10] = 258;
  result.matrix[1][11] = -894;
  result.matrix[2][0] = 1859;
  result.matrix[2][1] = 34;
  result.matrix[2][2] = -1144;
  result.matrix[2][3] = 405;
  result.matrix[2][4] = -1479;
  result.matrix[2][5] = -1196;
  result.matrix[2][6] = -1037;
  result.matrix[2][7] = -909;
  result.matrix[2][8] = 405;
  result.matrix[2][9] = 271;
  result.matrix[2][10] = 722;
  result.matrix[2][11] = 63;
  result.matrix[3][0] = 442;
  result.matrix[3][1] = -20;
  result.matrix[3][2] = 146;
  result.matrix[3][3] = -1598;
  result.matrix[3][4] = 652;
  result.matrix[3][5] = -387;
  result.matrix[3][6] = 431;
  result.matrix[3][7] = 32;
  result.matrix[3][8] = -466;
  result.matrix[3][9] = -160;
  result.matrix[3][10] = 566;
  result.matrix[3][11] = -500;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);
  ck_assert_double_eq(result.matrix[0][5], s21_result.matrix[0][5]);
  ck_assert_double_eq(result.matrix[0][6], s21_result.matrix[0][6]);
  ck_assert_double_eq(result.matrix[0][7], s21_result.matrix[0][7]);
  ck_assert_double_eq(result.matrix[0][8], s21_result.matrix[0][8]);
  ck_assert_double_eq(result.matrix[0][9], s21_result.matrix[0][9]);
  ck_assert_double_eq(result.matrix[0][10], s21_result.matrix[0][10]);
  ck_assert_double_eq(result.matrix[0][11], s21_result.matrix[0][11]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[1][4], s21_result.matrix[1][4]);
  ck_assert_double_eq(result.matrix[1][5], s21_result.matrix[1][5]);
  ck_assert_double_eq(result.matrix[1][6], s21_result.matrix[1][6]);
  ck_assert_double_eq(result.matrix[1][7], s21_result.matrix[1][7]);
  ck_assert_double_eq(result.matrix[1][8], s21_result.matrix[1][8]);
  ck_assert_double_eq(result.matrix[1][9], s21_result.matrix[1][9]);
  ck_assert_double_eq(result.matrix[1][10], s21_result.matrix[1][10]);
  ck_assert_double_eq(result.matrix[1][11], s21_result.matrix[1][11]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[2][4], s21_result.matrix[2][4]);
  ck_assert_double_eq(result.matrix[2][5], s21_result.matrix[2][5]);
  ck_assert_double_eq(result.matrix[2][6], s21_result.matrix[2][6]);
  ck_assert_double_eq(result.matrix[2][7], s21_result.matrix[2][7]);
  ck_assert_double_eq(result.matrix[2][8], s21_result.matrix[2][8]);
  ck_assert_double_eq(result.matrix[2][9], s21_result.matrix[2][9]);
  ck_assert_double_eq(result.matrix[2][10], s21_result.matrix[2][10]);
  ck_assert_double_eq(result.matrix[2][11], s21_result.matrix[2][11]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[3][4], s21_result.matrix[3][4]);
  ck_assert_double_eq(result.matrix[3][5], s21_result.matrix[3][5]);
  ck_assert_double_eq(result.matrix[3][6], s21_result.matrix[3][6]);
  ck_assert_double_eq(result.matrix[3][7], s21_result.matrix[3][7]);
  ck_assert_double_eq(result.matrix[3][8], s21_result.matrix[3][8]);
  ck_assert_double_eq(result.matrix[3][9], s21_result.matrix[3][9]);
  ck_assert_double_eq(result.matrix[3][10], s21_result.matrix[3][10]);
  ck_assert_double_eq(result.matrix[3][11], s21_result.matrix[3][11]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_13_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(1, 5, &A);
  s21_create_matrix(1, 5, &B);
  s21_create_matrix(1, 5, &result);
  A.matrix[0][0] = -235;
  A.matrix[0][1] = 64;
  A.matrix[0][2] = -46;
  A.matrix[0][3] = -952;
  A.matrix[0][4] = -646;
  B.matrix[0][0] = 772;
  B.matrix[0][1] = 794;
  B.matrix[0][2] = 239;
  B.matrix[0][3] = 409;
  B.matrix[0][4] = -964;
  result.matrix[0][0] = 537;
  result.matrix[0][1] = 858;
  result.matrix[0][2] = 193;
  result.matrix[0][3] = -543;
  result.matrix[0][4] = -1610;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_14_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(4, 5, &A);
  s21_create_matrix(4, 5, &B);
  s21_create_matrix(4, 5, &result);
  A.matrix[0][0] = 817;
  A.matrix[0][1] = -947;
  A.matrix[0][2] = 283;
  A.matrix[0][3] = -901;
  A.matrix[0][4] = 21;
  A.matrix[1][0] = -26;
  A.matrix[1][1] = 455;
  A.matrix[1][2] = 654;
  A.matrix[1][3] = 712;
  A.matrix[1][4] = -216;
  A.matrix[2][0] = 166;
  A.matrix[2][1] = -462;
  A.matrix[2][2] = -532;
  A.matrix[2][3] = 602;
  A.matrix[2][4] = -970;
  A.matrix[3][0] = -864;
  A.matrix[3][1] = -821;
  A.matrix[3][2] = 909;
  A.matrix[3][3] = -20;
  A.matrix[3][4] = 161;
  B.matrix[0][0] = 98;
  B.matrix[0][1] = -62;
  B.matrix[0][2] = 63;
  B.matrix[0][3] = -376;
  B.matrix[0][4] = -619;
  B.matrix[1][0] = 121;
  B.matrix[1][1] = 908;
  B.matrix[1][2] = -806;
  B.matrix[1][3] = -117;
  B.matrix[1][4] = 435;
  B.matrix[2][0] = -495;
  B.matrix[2][1] = -560;
  B.matrix[2][2] = -570;
  B.matrix[2][3] = 643;
  B.matrix[2][4] = 648;
  B.matrix[3][0] = 739;
  B.matrix[3][1] = -304;
  B.matrix[3][2] = 982;
  B.matrix[3][3] = -110;
  B.matrix[3][4] = -456;
  result.matrix[0][0] = 915;
  result.matrix[0][1] = -1009;
  result.matrix[0][2] = 346;
  result.matrix[0][3] = -1277;
  result.matrix[0][4] = -598;
  result.matrix[1][0] = 95;
  result.matrix[1][1] = 1363;
  result.matrix[1][2] = -152;
  result.matrix[1][3] = 595;
  result.matrix[1][4] = 219;
  result.matrix[2][0] = -329;
  result.matrix[2][1] = -1022;
  result.matrix[2][2] = -1102;
  result.matrix[2][3] = 1245;
  result.matrix[2][4] = -322;
  result.matrix[3][0] = -125;
  result.matrix[3][1] = -1125;
  result.matrix[3][2] = 1891;
  result.matrix[3][3] = -130;
  result.matrix[3][4] = -295;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[0][3], s21_result.matrix[0][3]);
  ck_assert_double_eq(result.matrix[0][4], s21_result.matrix[0][4]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[1][3], s21_result.matrix[1][3]);
  ck_assert_double_eq(result.matrix[1][4], s21_result.matrix[1][4]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);
  ck_assert_double_eq(result.matrix[2][3], s21_result.matrix[2][3]);
  ck_assert_double_eq(result.matrix[2][4], s21_result.matrix[2][4]);
  ck_assert_double_eq(result.matrix[3][0], s21_result.matrix[3][0]);
  ck_assert_double_eq(result.matrix[3][1], s21_result.matrix[3][1]);
  ck_assert_double_eq(result.matrix[3][2], s21_result.matrix[3][2]);
  ck_assert_double_eq(result.matrix[3][3], s21_result.matrix[3][3]);
  ck_assert_double_eq(result.matrix[3][4], s21_result.matrix[3][4]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(sum_matrix_15_by_py) {
  int res, s21_res;
  matrix_t A, B, result, s21_result;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_create_matrix(3, 3, &result);
  A.matrix[0][0] = 649;
  A.matrix[0][1] = 657;
  A.matrix[0][2] = -325;
  A.matrix[1][0] = 452;
  A.matrix[1][1] = -816;
  A.matrix[1][2] = 329;
  A.matrix[2][0] = -357;
  A.matrix[2][1] = -966;
  A.matrix[2][2] = -115;
  B.matrix[0][0] = 924;
  B.matrix[0][1] = 615;
  B.matrix[0][2] = 801;
  B.matrix[1][0] = 477;
  B.matrix[1][1] = -783;
  B.matrix[1][2] = -546;
  B.matrix[2][0] = 873;
  B.matrix[2][1] = -945;
  B.matrix[2][2] = 943;
  result.matrix[0][0] = 1573;
  result.matrix[0][1] = 1272;
  result.matrix[0][2] = 476;
  result.matrix[1][0] = 929;
  result.matrix[1][1] = -1599;
  result.matrix[1][2] = -217;
  result.matrix[2][0] = 516;
  result.matrix[2][1] = -1911;
  result.matrix[2][2] = 828;

  res = OK;
  s21_res = s21_sum_matrix(&A, &B, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  ck_assert_double_eq(result.matrix[2][0], s21_result.matrix[2][0]);
  ck_assert_double_eq(result.matrix[2][1], s21_result.matrix[2][1]);
  ck_assert_double_eq(result.matrix[2][2], s21_result.matrix[2][2]);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(transpose_1) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &result);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  result.matrix[0][0] = 1;
  result.matrix[0][1] = 2;
  result.matrix[0][2] = 3;
  result.matrix[1][0] = 4;
  result.matrix[1][1] = 5;
  result.matrix[1][2] = 6;

  res = OK;
  s21_res = s21_transpose(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  ck_assert_double_eq(result.matrix[0][1], s21_result.matrix[0][1]);
  ck_assert_double_eq(result.matrix[0][2], s21_result.matrix[0][2]);
  ck_assert_double_eq(result.matrix[1][0], s21_result.matrix[1][0]);
  ck_assert_double_eq(result.matrix[1][1], s21_result.matrix[1][1]);
  ck_assert_double_eq(result.matrix[1][2], s21_result.matrix[1][2]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

START_TEST(transpose_2) {
  int res, s21_res;
  matrix_t A, s21_result;
  s21_create_matrix(3, 2, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = S21_NAN;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;

  res = CALCULATION_ERROR;
  s21_res = s21_transpose(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(transpose_3) {
  int res, s21_res;
  matrix_t A, s21_result;
  s21_create_matrix(3, 2, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = S21_INFINITY;

  res = CALCULATION_ERROR;
  s21_res = s21_transpose(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(transpose_4) {
  int res, s21_res;
  matrix_t A, s21_result;
  s21_create_matrix(3, 2, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = S21_INFINITY_NEGATIVE;

  res = CALCULATION_ERROR;
  s21_res = s21_transpose(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(transpose_5) {
  int res, s21_res;
  matrix_t A = {.rows = 0, .columns = 3}, s21_result;

  res = INCORRECT_MATRIX;
  s21_res = s21_transpose(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
}
END_TEST

START_TEST(transpose_6) {
  int res, s21_res;
  matrix_t A, result, s21_result;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &result);

  A.matrix[0][0] = 3;
  result.matrix[0][0] = 3;

  res = OK;
  s21_res = s21_transpose(&A, &s21_result);

  ck_assert_int_eq(res, s21_res);
  ck_assert_int_eq(result.rows, s21_result.rows);
  ck_assert_int_eq(result.columns, s21_result.columns);
  ck_assert_double_eq(result.matrix[0][0], s21_result.matrix[0][0]);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&s21_result);
}
END_TEST

Suite *suite_transpose() {
  Suite *s = suite_create("s21_transpose");
  TCase *tc = tcase_create("transpose_tc");
  tcase_add_test(tc, transpose_1);
  tcase_add_test(tc, transpose_2);
  tcase_add_test(tc, transpose_3);
  tcase_add_test(tc, transpose_4);
  tcase_add_test(tc, transpose_5);
  tcase_add_test(tc, transpose_6);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_sum_matrix() {
  Suite *s = suite_create("s21_sum_matrix");
  TCase *tc = tcase_create("sum_matrix_tc");
  tcase_add_test(tc, sum_matrix_1);
  tcase_add_test(tc, sum_matrix_2);
  tcase_add_test(tc, sum_matrix_3);
  tcase_add_test(tc, sum_matrix_4);
  tcase_add_test(tc, sum_matrix_5);
  tcase_add_test(tc, sum_matrix_6);
  tcase_add_test(tc, sum_matrix_7);
  tcase_add_test(tc, sum_matrix_8);
  tcase_add_test(tc, sum_matrix_9);
  tcase_add_test(tc, sum_matrix_10_by_py);
  tcase_add_test(tc, sum_matrix_11_by_py);
  tcase_add_test(tc, sum_matrix_12_by_py);
  tcase_add_test(tc, sum_matrix_13_by_py);
  tcase_add_test(tc, sum_matrix_14_by_py);
  tcase_add_test(tc, sum_matrix_15_by_py);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_sub_matrix() {
  Suite *s = suite_create("s21_sub_matrix");
  TCase *tc = tcase_create("sub_matrix_tc");
  tcase_add_test(tc, sub_matrix_1);
  tcase_add_test(tc, sub_matrix_2);
  tcase_add_test(tc, sub_matrix_3);
  tcase_add_test(tc, sub_matrix_4);
  tcase_add_test(tc, sub_matrix_5);
  tcase_add_test(tc, sub_matrix_6);
  tcase_add_test(tc, sub_matrix_7);
  tcase_add_test(tc, sub_matrix_8);
  tcase_add_test(tc, sub_matrix_9);
  tcase_add_test(tc, sub_matrix_10_by_py);
  tcase_add_test(tc, sub_matrix_11_by_py);
  tcase_add_test(tc, sub_matrix_12_by_py);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_remove_matrix() {
  Suite *s = suite_create("s21_remove_matrix");
  TCase *tc = tcase_create("remove_matrix_tc");
  tcase_add_test(tc, remove_matrix_1);
  tcase_add_test(tc, remove_matrix_2);
  tcase_add_test(tc, remove_matrix_3);
  tcase_add_test(tc, remove_matrix_4);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_mult_number() {
  Suite *s = suite_create("s21_mult_number");
  TCase *tc = tcase_create("mult_number_tc");
  tcase_add_test(tc, mult_number_1);
  tcase_add_test(tc, mult_number_2);
  tcase_add_test(tc, mult_number_3);
  tcase_add_test(tc, mult_number_4);
  tcase_add_test(tc, mult_number_5);
  tcase_add_test(tc, mult_number_6);
  tcase_add_test(tc, mult_number_7);
  tcase_add_test(tc, mult_number_8);
  tcase_add_test(tc, mult_number_9);
  tcase_add_test(tc, mult_number_10);
  tcase_add_test(tc, mult_number_11_by_py);
  tcase_add_test(tc, mult_number_12_by_py);
  tcase_add_test(tc, mult_number_13_by_py);
  tcase_add_test(tc, mult_number_14_by_py);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_mult_matrix() {
  Suite *s = suite_create("s21_mult_matrix");
  TCase *tc = tcase_create("mult_matrix_tc");
  tcase_add_test(tc, mult_matrix_1);
  tcase_add_test(tc, mult_matrix_2);
  tcase_add_test(tc, mult_matrix_3);
  tcase_add_test(tc, mult_matrix_4);
  tcase_add_test(tc, mult_matrix_5);
  tcase_add_test(tc, mult_matrix_5_5);
  tcase_add_test(tc, mult_matrix_6);
  tcase_add_test(tc, mult_matrix_7);
  tcase_add_test(tc, mult_matrix_8);
  tcase_add_test(tc, mult_matrix_9);
  tcase_add_test(tc, mult_matrix_10);
  tcase_add_test(tc, mult_matrix_11_by_py);
  tcase_add_test(tc, mult_matrix_12_by_py);
  tcase_add_test(tc, mult_matrix_13_by_py);
  tcase_add_test(tc, mult_matrix_14_by_py);
  tcase_add_test(tc, mult_matrix_15_by_py);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_inverse_matrix() {
  Suite *s = suite_create("s21_inverse_matrix");
  TCase *tc = tcase_create("inverse_matrix_tc");
  tcase_add_test(tc, inverse_matrix_1);
  tcase_add_test(tc, inverse_matrix_2);
  tcase_add_test(tc, inverse_matrix_3);
  tcase_add_test(tc, inverse_matrix_4);
  tcase_add_test(tc, inverse_matrix_5);
  tcase_add_test(tc, inverse);
  tcase_add_test(tc, s21_inverse_1);
  tcase_add_test(tc, determinant);
  tcase_add_test(tc, test_incorrect);
  tcase_add_test(tc, inverse_matrix_6_by_py);
  tcase_add_test(tc, inverse_matrix_7_by_py);
  tcase_add_test(tc, inverse_matrix_8_by_py);
  tcase_add_test(tc, inverse_matrix_9_by_py);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_eq_matrix() {
  Suite *s = suite_create("s21_eq_matrix");
  TCase *tc = tcase_create("eq_matrix_tc");
  tcase_add_test(tc, eq_matrix_1);
  tcase_add_test(tc, eq_matrix_2);
  tcase_add_test(tc, eq_matrix_3);
  tcase_add_test(tc, eq_matrix_4);
  tcase_add_test(tc, eq_matrix_5);
  tcase_add_test(tc, eq_matrix_6);
  tcase_add_test(tc, eq_matrix_7);
  tcase_add_test(tc, eq_matrix_8);

  tcase_add_test(tc, zero_matrix);
  tcase_add_test(tc, zero_matrix_1);
  tcase_add_test(tc, casual_matrix_1);
  tcase_add_test(tc, casual_matrix_2);
  tcase_add_test(tc, casual_matrix_3);
  tcase_add_test(tc, casual_matrix_4);
  tcase_add_test(tc, casual_matrix_5);
  tcase_add_test(tc, casual_matrix_6);
  tcase_add_test(tc, casual_matrix_7);
  tcase_add_test(tc, casual_matrix_8);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_determinant() {
  Suite *s = suite_create("s21_determinant");
  TCase *tc = tcase_create("determinant_tc");
  tcase_add_test(tc, determinant_1);
  tcase_add_test(tc, determinant_2);
  tcase_add_test(tc, determinant_3);
  tcase_add_test(tc, determinant_4);
  tcase_add_test(tc, determinant_5);
  tcase_add_test(tc, determinant_6);
  tcase_add_test(tc, determinant_7);
  tcase_add_test(tc, determinant_8);
  tcase_add_test(tc, determinant_9);
  tcase_add_test(tc, determinant_10);
  tcase_add_test(tc, determinant_11_by_py);
  tcase_add_test(tc, determinant_12_by_py);
  tcase_add_test(tc, determinant_13_by_py);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_create_matrix() {
  Suite *s = suite_create("s21_create_matrix");
  TCase *tc = tcase_create("create_matrix_tc");
  tcase_add_test(tc, create_matrix_1);
  tcase_add_test(tc, create_matrix_2);
  tcase_add_test(tc, create_matrix_3);
  tcase_add_test(tc, create_matrix_4);
  tcase_add_test(tc, create_matrix_5);
  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_calc_complements() {
  Suite *s = suite_create("s21_calc_complements");
  TCase *tc = tcase_create("calc_complements_tc");
  tcase_add_test(tc, calc_complements_1);
  tcase_add_test(tc, calc_complements_2);
  tcase_add_test(tc, calc_complements_3);
  tcase_add_test(tc, calc_complements_4);
  tcase_add_test(tc, calc_complements_5);
  tcase_add_test(tc, calc_complements_5_5);
  tcase_add_test(tc, calc_complements_6);
  tcase_add_test(tc, calc_complements_7);
  suite_add_tcase(s, tc);
  return s;
}

int main() {
  Suite *tests[] = {
      suite_create_matrix(), suite_remove_matrix(),  suite_eq_matrix(),
      suite_sum_matrix(),    suite_sub_matrix(),     suite_mult_number(),
      suite_mult_matrix(),   suite_transpose(),      suite_calc_complements(),
      suite_determinant(),   suite_inverse_matrix(), NULL,
  };

  for (Suite **current = tests; *current; current++) {
    SRunner *sr = srunner_create(*current);
    srunner_set_fork_status(sr, CK_NOFORK);
    srunner_run_all(sr, CK_NORMAL);
    srunner_free(sr);
  }

  return 0;
}
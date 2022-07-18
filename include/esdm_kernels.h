/*
    ESDM-PAV Analytical Kernels
    Copyright (C) 2022 CMCC Foundation

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __ESDM_READ_STREAM_H
#define __ESDM_READ_STREAM_H

#include <esdm.h>

#define ESDM_SEPARATOR ","

#define ESDM_FUNCTION_NOP "nop"
#define ESDM_FUNCTION_STREAM "stream"

#define ESDM_FUNCTION_MAX "max"
#define ESDM_FUNCTION_MIN "min"
#define ESDM_FUNCTION_AVG "avg"
#define ESDM_FUNCTION_SUM "sum"
#define ESDM_FUNCTION_STD "std"
#define ESDM_FUNCTION_VAR "var"
#define ESDM_FUNCTION_STAT "stat"

#define ESDM_FUNCTION_OUTLIER "outlier"

#define ESDM_FUNCTION_SUM_SCALAR "sum_scalar"
#define ESDM_FUNCTION_MUL_SCALAR "mul_scalar"

#define ESDM_FUNCTION_ABS "abs"
#define ESDM_FUNCTION_SQR "sqr"
#define ESDM_FUNCTION_SQRT "sqrt"
#define ESDM_FUNCTION_CEIL "ceil"
#define ESDM_FUNCTION_FLOOR "floor"
#define ESDM_FUNCTION_ROUND "round"
#define ESDM_FUNCTION_INT "int"
#define ESDM_FUNCTION_NINT "nint"

#define ESDM_FUNCTION_POW "pow"
#define ESDM_FUNCTION_EXP "exp"
#define ESDM_FUNCTION_LOG "log"
#define ESDM_FUNCTION_LOG10 "log10"

#define ESDM_FUNCTION_SIN "sin"
#define ESDM_FUNCTION_COS "cos"
#define ESDM_FUNCTION_TAN "tan"
#define ESDM_FUNCTION_ASIN "asin"
#define ESDM_FUNCTION_ACOS "acos"
#define ESDM_FUNCTION_ATAN "atan"
#define ESDM_FUNCTION_SINH "sinh"
#define ESDM_FUNCTION_COSH "cosh"
#define ESDM_FUNCTION_TANH "tanh"

#define ESDM_FUNCTION_RECI "reci"
#define ESDM_FUNCTION_NOT "not"

typedef struct _esdm_stream_data_t {
	char *operation;
	char *args;
	void *buff;
	char valid;
	double value1;
	double value2;
	uint64_t number;
	void *fill_value;
} esdm_stream_data_t;

int esdm_is_a_reduce_func(const char *operation, const char *args);
void *esdm_stream_func(esdm_dataspace_t * space, void *buff, void *user_ptr, void *esdm_fill_value);
void esdm_reduce_func(esdm_dataspace_t * space, void *user_ptr, void *stream_func_out);

#endif				//__ESDM_READ_STREAM_H

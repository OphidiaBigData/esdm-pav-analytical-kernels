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

#include <stdlib.h>
#include <math.h>

#include "esdm_kernels.h"

#define UNUSED(x) {(void)(x);}

#define ESDM_FUNCTION_OP_N 3
#define ESDM_FUNCTION_OP_SET '1'
#define ESDM_FUNCTION_OP_LESS_THAN '<'
#define ESDM_FUNCTION_OP_MORE_THAN '>'

typedef struct _esdm_stream_data_out_t {
	double value1;
	double value2;
	double value3;
	uint64_t number;
} esdm_stream_data_out_t;

int esdm_is_a_reduce_func(const char *operation, const char *args)
{
	if (!operation)
		return 0;

	if (!strcmp(operation, ESDM_FUNCTION_MAX))
		return 1;
	if (!strcmp(operation, ESDM_FUNCTION_MIN))
		return 1;
	if (!strcmp(operation, ESDM_FUNCTION_AVG))
		return 1;
	if (!strcmp(operation, ESDM_FUNCTION_SUM))
		return 1;
	if (!strcmp(operation, ESDM_FUNCTION_STD))
		return 1;
	if (!strcmp(operation, ESDM_FUNCTION_VAR))
		return 1;
	if (!strcmp(operation, ESDM_FUNCTION_OUTLIER))
		return 1;
	if (!strcmp(operation, ESDM_FUNCTION_STAT)) {
		int i, option = 0;
		for (i = 0; i < ESDM_FUNCTION_OP_N; ++i)
			if (args && args[0]) {
				if (args[0] == ESDM_FUNCTION_OP_SET)
					option++;
				args++;
			} else
				break;
		return option;
	}

	return 0;
}

void *esdm_stream_func(esdm_dataspace_t * space, void *buff, void *user_ptr, void *esdm_fill_value)
{
	UNUSED(esdm_fill_value);

	if (!space || !buff || !user_ptr)
		return NULL;

	esdm_stream_data_t *stream_data = (esdm_stream_data_t *) user_ptr;
	if (!stream_data->operation)
		return NULL;

	char *args = stream_data->args ? strdup(stream_data->args) : NULL;	// Copy for strtok

	int64_t i, idx, ndims = esdm_dataspace_get_dims(space);
	int64_t const *s = esdm_dataspace_get_size(space);
	//int64_t const *si = esdm_dataspace_get_offset(space);
	int64_t ci[ndims], ei[ndims];
	for (i = 0; i < ndims; ++i) {
		ci[i] = 0;	// + si[i]
		ei[i] = s[i];	// + si[i]
	}

	uint64_t k = 1, n = esdm_dataspace_element_count(space);
	esdm_type_t type = esdm_dataspace_get_type(space);
	void *fill_value = stream_data->fill_value;
	esdm_stream_data_out_t *tmp = NULL;

	if (!strcmp(stream_data->operation, ESDM_FUNCTION_NOP) || !strcmp(stream_data->operation, ESDM_FUNCTION_STREAM)) {

		// TODO: copy only the data related to the dataspace
		memcpy(stream_data->buff, buff, esdm_dataspace_total_bytes(space));

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_MAX)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 0;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v < a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v < a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v < a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v < a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v < a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v < a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_MIN)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 0;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v > a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v > a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v > a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v > a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v > a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if ((!fill_value || (a[idx] != fv)) && (!tmp->number || (v > a[idx]))) {
					v = a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_AVG) || !strcmp(stream_data->operation, ESDM_FUNCTION_SUM)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->value1 = 0;
		tmp->number = 0;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, fv = fill_value ? *(char *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, fv = fill_value ? *(short *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, fv = fill_value ? *(int *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, fv = fill_value ? *(long long *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, fv = fill_value ? *(float *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, fv = fill_value ? *(double *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_STD) || !strcmp(stream_data->operation, ESDM_FUNCTION_VAR)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->value1 = 0;
		tmp->value2 = 0;
		tmp->number = 0;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, fv = fill_value ? *(char *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->value2 += a[idx] * a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, fv = fill_value ? *(short *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->value2 += a[idx] * a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, fv = fill_value ? *(int *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->value2 += a[idx] * a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, fv = fill_value ? *(long long *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->value2 += a[idx] * a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, fv = fill_value ? *(float *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->value2 += a[idx] * a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, fv = fill_value ? *(double *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					tmp->value1 += a[idx];
					tmp->value2 += a[idx] * a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_STAT)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->value3 = 0;
		tmp->number = 0;

		char *save_pointer = NULL, *arg = args ? strtok_r(args, ESDM_SEPARATOR, &save_pointer) : NULL;
		char i, option = 0, tbyte = 1;
		for (i = 0; i < ESDM_FUNCTION_OP_N; ++i)
			if (arg && arg[0]) {
				if (arg[0] == ESDM_FUNCTION_OP_SET)
					option |= tbyte;
				arg++;
				tbyte <<= 1;
			} else
				break;
		if (!option) {

			// No operation is executed in this case

		} else if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v1 = 0, v2 = 0, fv = fill_value ? *(char *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					if ((option & 1) && (!tmp->number || (v1 < a[idx])))	// Min
						v1 = a[idx];
					if ((option & 2) && (!tmp->number || (v2 > a[idx])))	// Max
						v2 = a[idx];
					if (option & 4)	// Avg
						tmp->value3 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v1;
			tmp->value2 = v2;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v1 = 0, v2 = 0, fv = fill_value ? *(short *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					if ((option & 1) && (!tmp->number || (v1 < a[idx])))	// Min
						v1 = a[idx];
					if ((option & 2) && (!tmp->number || (v2 > a[idx])))	// Max
						v2 = a[idx];
					if (option & 4)	// Avg
						tmp->value3 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v1;
			tmp->value2 = v2;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v1 = 0, v2 = 0, fv = fill_value ? *(int *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					if ((option & 1) && (!tmp->number || (v1 < a[idx])))	// Min
						v1 = a[idx];
					if ((option & 2) && (!tmp->number || (v2 > a[idx])))	// Max
						v2 = a[idx];
					if (option & 4)	// Avg
						tmp->value3 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v1;
			tmp->value2 = v2;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v1 = 0, v2 = 0, fv = fill_value ? *(long long *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					if ((option & 1) && (!tmp->number || (v1 < a[idx])))	// Min
						v1 = a[idx];
					if ((option & 2) && (!tmp->number || (v2 > a[idx])))	// Max
						v2 = a[idx];
					if (option & 4)	// Avg
						tmp->value3 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v1;
			tmp->value2 = v2;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v1 = 0, v2 = 0, fv = fill_value ? *(float *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					if ((option & 1) && (!tmp->number || (v1 < a[idx])))	// Min
						v1 = a[idx];
					if ((option & 2) && (!tmp->number || (v2 > a[idx])))	// Max
						v2 = a[idx];
					if (option & 4)	// Avg
						tmp->value3 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v1;
			tmp->value2 = v2;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v1 = 0, v2 = 0, fv = fill_value ? *(double *) fill_value : 0;
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv)) {
					if ((option & 1) && (!tmp->number || (v1 < a[idx])))	// Min
						v1 = a[idx];
					if ((option & 2) && (!tmp->number || (v2 > a[idx])))	// Max
						v2 = a[idx];
					if (option & 4)	// Avg
						tmp->value3 += a[idx];
					tmp->number++;
				}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v1;
			tmp->value2 = v2;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_OUTLIER)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->value1 = 0;
		tmp->number = 1;	// Use only to avoid errors during the reduction phase

		char *save_pointer = NULL, *arg = args ? strtok_r(args, ESDM_SEPARATOR, &save_pointer) : NULL;
		char thresh_type = ESDM_FUNCTION_OP_MORE_THAN;
		if (arg) {
			if (!isdigit(arg[0])) {
				if (arg[0] == ESDM_FUNCTION_OP_LESS_THAN)
					thresh_type = ESDM_FUNCTION_OP_LESS_THAN;
				arg++;
			}
		}
		if (!arg || !arg[0]) {

			// No element is considered outlier in case the threshold is not given

		} else if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, fv = fill_value ? *(char *) fill_value : 0, v = strtol(arg, NULL, 10);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv))
					switch (thresh_type) {
						case ESDM_FUNCTION_OP_LESS_THAN:
							if (v > a[idx])
								tmp->value1++;
							break;
						case ESDM_FUNCTION_OP_MORE_THAN:
						default:
							if (v < a[idx])
								tmp->value1++;
					}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, fv = fill_value ? *(short *) fill_value : 0, v = strtol(arg, NULL, 10);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv))
					switch (thresh_type) {
						case ESDM_FUNCTION_OP_LESS_THAN:
							if (v > a[idx])
								tmp->value1++;
							break;
						case ESDM_FUNCTION_OP_MORE_THAN:
						default:
							if (v < a[idx])
								tmp->value1++;
					}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, fv = fill_value ? *(int *) fill_value : 0, v = strtol(arg, NULL, 10);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv))
					switch (thresh_type) {
						case ESDM_FUNCTION_OP_LESS_THAN:
							if (v > a[idx])
								tmp->value1++;
							break;
						case ESDM_FUNCTION_OP_MORE_THAN:
						default:
							if (v < a[idx])
								tmp->value1++;
					}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, fv = fill_value ? *(long long *) fill_value : 0, v = strtoll(arg, NULL, 10);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv))
					switch (thresh_type) {
						case ESDM_FUNCTION_OP_LESS_THAN:
							if (v > a[idx])
								tmp->value1++;
							break;
						case ESDM_FUNCTION_OP_MORE_THAN:
						default:
							if (v < a[idx])
								tmp->value1++;
					}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, fv = fill_value ? *(float *) fill_value : 0, v = strtof(arg, NULL);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv))
					switch (thresh_type) {
						case ESDM_FUNCTION_OP_LESS_THAN:
							if (v > a[idx])
								tmp->value1++;
							break;
						case ESDM_FUNCTION_OP_MORE_THAN:
						default:
							if (v < a[idx])
								tmp->value1++;
					}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, fv = fill_value ? *(double *) fill_value : 0, v = strtod(arg, NULL);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				if (!fill_value || (a[idx] != fv))
					switch (thresh_type) {
						case ESDM_FUNCTION_OP_LESS_THAN:
							if (v > a[idx])
								tmp->value1++;
							break;
						case ESDM_FUNCTION_OP_MORE_THAN:
						default:
							if (v < a[idx])
								tmp->value1++;
					}
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_SUM_SCALAR)) {

		if (!args) {
			// TODO: copy only the data related to the dataspace
			memcpy(stream_data->buff, buff, esdm_dataspace_total_bytes(space));
			return NULL;
		}

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		char *save_pointer = NULL, *arg = args ? strtok_r(args, ESDM_SEPARATOR, &save_pointer) : NULL;

		if (type == SMD_DTYPE_INT8) {

			char scalar = arg ? strtol(arg, NULL, 10) : 0;
			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] + scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short scalar = arg ? strtol(arg, NULL, 10) : 0;
			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] + scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int scalar = arg ? strtol(arg, NULL, 10) : 0;
			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] + scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long scalar = arg ? strtoll(arg, NULL, 10) : 0;
			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] + scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float scalar = arg ? strtof(arg, NULL) : 0;
			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] + scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double scalar = arg ? strtod(arg, NULL) : 0;
			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] + scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_MUL_SCALAR)) {

		if (!args) {
			// TODO: copy only the data related to the dataspace
			memcpy(stream_data->buff, buff, esdm_dataspace_total_bytes(space));
			return NULL;
		}

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		char *save_pointer = NULL, *arg = args ? strtok_r(args, ESDM_SEPARATOR, &save_pointer) : NULL;

		if (type == SMD_DTYPE_INT8) {

			char scalar = arg ? strtol(arg, NULL, 10) : 1;
			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short scalar = arg ? strtol(arg, NULL, 10) : 1;
			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int scalar = arg ? strtol(arg, NULL, 10) : 1;
			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long scalar = arg ? strtoll(arg, NULL, 10) : 1;
			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float scalar = arg ? strtof(arg, NULL) : 1;
			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double scalar = arg ? strtod(arg, NULL) : 1;
			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * scalar : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_ABS)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? abs(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? abs(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? abs(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? abs(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? fabs(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? fabs(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_SQRT)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sqrt(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sqrt(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sqrt(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sqrt(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sqrt(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sqrt(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_SQR)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? a[idx] * a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

		// TODO: to be optimized for integer values
	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_CEIL)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? ceil(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? ceil(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? ceil(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? ceil(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? ceil(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? ceil(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

		// TODO: to be optimized for integer values
	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_FLOOR) || !strcmp(stream_data->operation, ESDM_FUNCTION_INT)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

		// TODO: to be optimized for integer values
	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_ROUND) || !strcmp(stream_data->operation, ESDM_FUNCTION_NINT)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx] + 0.5) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx] + 0.5) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx] + 0.5) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx] + 0.5) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx] + 0.5) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? floor(a[idx] + 0.5) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_POW)) {

		if (!args) {
			// TODO: copy only the data related to the dataspace
			memcpy(stream_data->buff, buff, esdm_dataspace_total_bytes(space));
			return NULL;
		}

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		char *save_pointer = NULL, *arg = args ? strtok_r(args, ESDM_SEPARATOR, &save_pointer) : NULL;

		if (type == SMD_DTYPE_INT8) {

			char scalar = arg ? strtol(arg, NULL, 10) : 1;
			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? pow(a[idx], scalar) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short scalar = arg ? strtol(arg, NULL, 10) : 1;
			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? pow(a[idx], scalar) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int scalar = arg ? strtol(arg, NULL, 10) : 1;
			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? pow(a[idx], scalar) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long scalar = arg ? strtoll(arg, NULL, 10) : 1;
			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? pow(a[idx], scalar) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float scalar = arg ? strtof(arg, NULL) : 1;
			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? pow(a[idx], scalar) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double scalar = arg ? strtod(arg, NULL) : 1;
			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? pow(a[idx], scalar) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_EXP)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? exp(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? exp(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? exp(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? exp(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? exp(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? exp(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_LOG)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_LOG10)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log10(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log10(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log10(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log10(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log10(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? log10(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_SIN)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_COS)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_TAN)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_ASIN)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? asin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? asin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? asin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? asin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? asin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? asin(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_ACOS)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? acos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? acos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? acos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? acos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? acos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? acos(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_ATAN)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? atan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? atan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? atan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? atan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? atan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? atan(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_SINH)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sinh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sinh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sinh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sinh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sinh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? sinh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_COSH)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cosh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cosh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cosh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cosh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cosh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? cosh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_TANH)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tanh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tanh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tanh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tanh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tanh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? tanh(a[idx]) : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

		// TODO: to be optimized for integer values
	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_RECI)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? 1.0 / a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? 1.0 / a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? 1.0 / a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? 1.0 / a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? 1.0 / a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? 1.0 / a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_NOT)) {

		tmp = (esdm_stream_data_out_t *) malloc(sizeof(esdm_stream_data_out_t));
		tmp->number = 1;

		if (type == SMD_DTYPE_INT8) {

			char *a = (char *) buff, v = 0, fv = fill_value ? *(char *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? !a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT16) {

			short *a = (short *) buff, v = 0, fv = fill_value ? *(short *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? !a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT32) {

			int *a = (int *) buff, v = 0, fv = fill_value ? *(int *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? !a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_INT64) {

			long long *a = (long long *) buff, v = 0, fv = fill_value ? *(long long *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? !a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_FLOAT) {

			float *a = (float *) buff, v = 0, fv = fill_value ? *(float *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? !a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else if (type == SMD_DTYPE_DOUBLE) {

			double *a = (double *) buff, v = 0, fv = fill_value ? *(double *) fill_value : 0;
			size_t step = sizeof(v);
			for (k = 0; k < n; k++) {
				idx = 0;
				for (i = 0; i < ndims; i++)
					idx = idx * s[i] + ci[i];
				v = !fill_value || (a[idx] != fv) ? !a[idx] : fv;
				memcpy(stream_data->buff + idx * step, &v, step);
				for (i = ndims - 1; i >= 0; i--) {
					ci[i]++;
					if (ci[i] < ei[i])
						break;
					ci[i] = 0;	// si[i];
				}
			}
			tmp->value1 = v;

		} else {
			free(tmp);
			if (args)
				free(args);
			return NULL;
		}

	}

	if (args)
		free(args);

	return tmp;
}

void esdm_reduce_func(esdm_dataspace_t * space, void *user_ptr, void *stream_func_out)
{
	esdm_stream_data_out_t *tmp = (esdm_stream_data_out_t *) stream_func_out;

	do {

		if (!space || !user_ptr || (tmp && !tmp->number))
			break;

		esdm_type_t type = esdm_dataspace_get_type(space);
		esdm_stream_data_t *stream_data = (esdm_stream_data_t *) user_ptr;
		if (!stream_data->operation)
			break;

		if (!strcmp(stream_data->operation, ESDM_FUNCTION_MAX)) {

			if (!tmp)
				break;

			if (type == SMD_DTYPE_INT8) {

				char v = (char) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(char *) stream_data->buff;
					if (pre < v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_INT16) {

				short v = (short) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(short *) stream_data->buff;
					if (pre < v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_INT32) {

				int v = (int) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(int *) stream_data->buff;
					if (pre < v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_INT64) {

				long long v = (long long) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(long long *) stream_data->buff;
					if (pre < v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_FLOAT) {

				float v = (float) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(float *) stream_data->buff;
					if (pre < v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_DOUBLE) {

				double v = (double) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(double *) stream_data->buff;
					if (pre < v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			}

		} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_MIN)) {

			if (!tmp)
				break;

			if (type == SMD_DTYPE_INT8) {

				char v = (char) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(char *) stream_data->buff;
					if (pre > v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_INT16) {

				short v = (short) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(short *) stream_data->buff;
					if (pre > v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_INT32) {

				int v = (int) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(int *) stream_data->buff;
					if (pre > v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_INT64) {

				long long v = (long long) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(long long *) stream_data->buff;
					if (pre > v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_FLOAT) {

				float v = (float) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(float *) stream_data->buff;
					if (pre > v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			} else if (type == SMD_DTYPE_DOUBLE) {

				double v = (double) tmp->value1, pre;
				if (stream_data->valid) {
					pre = *(double *) stream_data->buff;
					if (pre > v)
						memcpy(stream_data->buff, &v, sizeof(v));
				} else {
					memcpy(stream_data->buff, &v, sizeof(v));
					stream_data->valid = 1;
				}

			}

		} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_AVG) || !strcmp(stream_data->operation, ESDM_FUNCTION_SUM) || !strcmp(stream_data->operation, ESDM_FUNCTION_OUTLIER)) {

			if (!tmp)
				break;

			if (!stream_data->valid) {
				stream_data->valid = 1;
				stream_data->value1 = 0;
				stream_data->number = 0;
			}
			stream_data->value1 += tmp->value1;
			if (!strcmp(stream_data->operation, ESDM_FUNCTION_AVG))
				stream_data->number += tmp->number;
			else
				stream_data->number = 1;

			if (type == SMD_DTYPE_INT8) {

				char v = (char) (stream_data->value1 / stream_data->number);
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_INT16) {

				short v = (short) (stream_data->value1 / stream_data->number);
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_INT32) {

				int v = (int) (stream_data->value1 / stream_data->number);
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_INT64) {

				long long v = (long long) (stream_data->value1 / stream_data->number);
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_FLOAT) {

				float v = (float) (stream_data->value1 / stream_data->number);
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_DOUBLE) {

				double v = (double) (stream_data->value1 / stream_data->number);
				memcpy(stream_data->buff, &v, sizeof(v));

			}

		} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_STD) || !strcmp(stream_data->operation, ESDM_FUNCTION_VAR)) {

			if (!tmp)
				break;

			if (!stream_data->valid) {
				stream_data->valid = 1;
				stream_data->value1 = 0;
				stream_data->value2 = 0;
				stream_data->number = 0;
			}
			stream_data->value1 += tmp->value1;
			stream_data->value2 += tmp->value2;
			stream_data->number += tmp->number;

			double result = (stream_data->value2 - stream_data->value1 * stream_data->value1 / stream_data->number) / stream_data->number;
			if (stream_data->number > 1)
				result *= stream_data->number / (stream_data->number - 1.0);
			if (!strcmp(stream_data->operation, ESDM_FUNCTION_STD))
				result = sqrt(result);

			if (type == SMD_DTYPE_INT8) {

				char v = (char) result;
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_INT16) {

				short v = (short) result;
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_INT32) {

				int v = (int) result;
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_INT64) {

				long long v = (long long) result;
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_FLOAT) {

				float v = (float) result;
				memcpy(stream_data->buff, &v, sizeof(v));

			} else if (type == SMD_DTYPE_DOUBLE) {

				double v = (double) result;
				memcpy(stream_data->buff, &v, sizeof(v));

			}

		} else if (!strcmp(stream_data->operation, ESDM_FUNCTION_STAT)) {

			if (!tmp)
				break;

			char i, option = 0, tbyte = 1, *arg = stream_data->args;
			size_t offset = 0;
			for (i = 0; i < ESDM_FUNCTION_OP_N; ++i)
				if (arg && arg[0]) {
					if (arg[0] == ESDM_FUNCTION_OP_SET)
						option |= tbyte;
					arg++;
					tbyte <<= 1;
				} else
					break;

			if (option & 4) {
				if (!stream_data->valid) {
					stream_data->value1 = 0;
					stream_data->number = 0;
				}
				stream_data->value1 += tmp->value3;
				stream_data->number += tmp->number;
			}

			if (type == SMD_DTYPE_INT8) {

				char v, pre;
				if (option & 1) {
					v = (char) tmp->value1;
					if (stream_data->valid) {
						pre = *(char *) stream_data->buff;
						if (pre > v)
							memcpy(stream_data->buff, &v, sizeof(v));
					} else
						memcpy(stream_data->buff, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 2) {
					v = (char) tmp->value2;
					if (stream_data->valid) {
						pre = *((char *) (stream_data->buff + offset));
						if (pre < v)
							memcpy(stream_data->buff + offset, &v, sizeof(v));
					} else
						memcpy(stream_data->buff + offset, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 4) {
					v = (char) (stream_data->value1 / stream_data->number);
					memcpy(stream_data->buff + offset, &v, sizeof(v));
					//offset += sizeof(v);  // Useless
				}
				stream_data->valid = 1;

			} else if (type == SMD_DTYPE_INT16) {

				short v, pre;
				if (option & 1) {
					v = (short) tmp->value1;
					if (stream_data->valid) {
						pre = *(short *) stream_data->buff;
						if (pre > v)
							memcpy(stream_data->buff, &v, sizeof(v));
					} else
						memcpy(stream_data->buff, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 2) {
					v = (short) tmp->value2;
					if (stream_data->valid) {
						pre = *((short *) (stream_data->buff + offset));
						if (pre < v)
							memcpy(stream_data->buff + offset, &v, sizeof(v));
					} else
						memcpy(stream_data->buff + offset, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 4) {
					v = (short) (stream_data->value1 / stream_data->number);
					memcpy(stream_data->buff + offset, &v, sizeof(v));
					//offset += sizeof(v);  // Useless
				}
				stream_data->valid = 1;

			} else if (type == SMD_DTYPE_INT32) {

				char v, pre;
				if (option & 1) {
					v = (int) tmp->value1;
					if (stream_data->valid) {
						pre = *(int *) stream_data->buff;
						if (pre > v)
							memcpy(stream_data->buff, &v, sizeof(v));
					} else
						memcpy(stream_data->buff, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 2) {
					v = (int) tmp->value2;
					if (stream_data->valid) {
						pre = *((int *) (stream_data->buff + offset));
						if (pre < v)
							memcpy(stream_data->buff + offset, &v, sizeof(v));
					} else
						memcpy(stream_data->buff + offset, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 4) {
					v = (int) (stream_data->value1 / stream_data->number);
					memcpy(stream_data->buff + offset, &v, sizeof(v));
					//offset += sizeof(v);  // Useless
				}
				stream_data->valid = 1;

			} else if (type == SMD_DTYPE_INT64) {

				long long v, pre;
				if (option & 1) {
					v = (long long) tmp->value1;
					if (stream_data->valid) {
						pre = *(long long *) stream_data->buff;
						if (pre > v)
							memcpy(stream_data->buff, &v, sizeof(v));
					} else
						memcpy(stream_data->buff, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 2) {
					v = (long long) tmp->value2;
					if (stream_data->valid) {
						pre = *((long long *) (stream_data->buff + offset));
						if (pre < v)
							memcpy(stream_data->buff + offset, &v, sizeof(v));
					} else
						memcpy(stream_data->buff + offset, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 4) {
					v = (long long) (stream_data->value1 / stream_data->number);
					memcpy(stream_data->buff + offset, &v, sizeof(v));
					//offset += sizeof(v);  // Useless
				}
				stream_data->valid = 1;

			} else if (type == SMD_DTYPE_FLOAT) {

				float v, pre;
				if (option & 1) {
					v = (float) tmp->value1;
					if (stream_data->valid) {
						pre = *(float *) stream_data->buff;
						if (pre > v)
							memcpy(stream_data->buff, &v, sizeof(v));
					} else
						memcpy(stream_data->buff, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 2) {
					v = (float) tmp->value2;
					if (stream_data->valid) {
						pre = *((float *) (stream_data->buff + offset));
						if (pre < v)
							memcpy(stream_data->buff + offset, &v, sizeof(v));
					} else
						memcpy(stream_data->buff + offset, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 4) {
					v = (float) (stream_data->value1 / stream_data->number);
					memcpy(stream_data->buff + offset, &v, sizeof(v));
					//offset += sizeof(v);  // Useless
				}
				stream_data->valid = 1;

			} else if (type == SMD_DTYPE_DOUBLE) {

				double v, pre;
				if (option & 1) {
					v = (double) tmp->value1;
					if (stream_data->valid) {
						pre = *(double *) stream_data->buff;
						if (pre > v)
							memcpy(stream_data->buff, &v, sizeof(v));
					} else
						memcpy(stream_data->buff, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 2) {
					v = (double) tmp->value2;
					if (stream_data->valid) {
						pre = *((double *) (stream_data->buff + offset));
						if (pre < v)
							memcpy(stream_data->buff + offset, &v, sizeof(v));
					} else
						memcpy(stream_data->buff + offset, &v, sizeof(v));
					offset += sizeof(v);
				}
				if (option & 4) {
					v = (double) (stream_data->value1 / stream_data->number);
					memcpy(stream_data->buff + offset, &v, sizeof(v));
					//offset += sizeof(v);  // Useless
				}
				stream_data->valid = 1;

			}

		}

	} while (0);

	if (stream_func_out)
		free(stream_func_out);
}

/*
 * -------------------------------------------------------------------------
 * validate.c -- Input validation function
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2006  Roger Liu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------------
 */

#include "swt_common.h"
#include "dwt.h"


str_error_notification strErrNoti[] ={
	{POSITIVE_INTEGER_ONLY, "Input integer must be positive!\n"},
	{LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION, "Length Parameter is not valid for input vector!\n"},
	{SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION, "Size Parameter is not valid for input matrix!\n"},
	{OPT_CHAR_NOT_VALID, "Optional Charactor Parameter is not valid!\n"},
	{EXTENSION_OPT_NOT_VALID, "Extension Method is not valid!\n"},
	{WAVELET_NAME_NOT_VALID, "Wavelet Name is not valid!\n"},
	{DECOMPOSITION_LEVEL_NOT_VALID, "Input signal is not valid for selected decompostion level and wavelets!\n"},
	{MULTI_DECOM_LEVEL_LESS_THAN_TWO ,"Decomposition level must be no less than 2!\n"},
	{UNKNOWN_INPUT_ERR, "Unrecognized Input and Output Numbers or parameter not valid for the algorithm!\nPlease refer to help pages!\n"},
	{WRONG_LHS, "Output parameter quantity is inappropriate!\n"}
};
int errorNum = sizeof(strErrNoti)/sizeof(str_error_notification);

/*-------------------------------------------
 * Dimension Checking
 *-----------------------------------------*/
void void_check (int number, int *type)
{
  int row, col;
  GetMatrixdims(number,&row,&col);
  if ((row==0) && (col==0))
    *type = 1;
  else
    *type = 0;
  //return;
}

void scalar_check (int number, int *type)
{
  int row, col;
  GetMatrixdims(number,&row,&col);
  if ((row==1) && (col==1))
    *type = 1;
  else
    *type = 0;
  //return;
}

void vector_check (int number, int *type)
{
  int row, col;
  GetMatrixdims(number,&row,&col);
  if ((row==1) && (col>1))
    *type = 1;
  else if ((row>1) && (col==1))
    *type = 1;
  else
    *type = 0;
  //return;
}

void matrix_check (int number, int *type)
{
  int row, col;
  GetMatrixdims(number,&row,&col);
  if ((row>1) && (col>1))
    *type = 1;
  else
    *type = 0;
  //return;
}

void real_or_complex (int number, int *type)
{
  int il, lw;
  lw = number + Top - Rhs;
  il = iadr(*Lstk(lw));
  *type = *istk(il+3); // 0 if real, 1 if complex
  //return;
}

int sci_matrix_vector_real (int number)
{
  int typ1, typ2;
  vector_check(number, &typ1);
  real_or_complex(number, &typ2);
  if (typ1 && (!typ2) && (GetType(number)==sci_matrix))
    return 1;
  else
    return 0;
}

int sci_matrix_vector_complex (int number)
{
  int typ1, typ2;
  vector_check(number, &typ1);
  real_or_complex(number, &typ2);
  if (typ1 && typ2 && (GetType(number)==sci_matrix))
    return 1;
  else
    return 0;
}

int sci_matrix_matrix_complex (int number)
{
  int typ1, typ2;
  matrix_check(number, &typ1);
  real_or_complex(number, &typ2);
  if (typ1 && typ2 && (GetType(number)==sci_matrix))
    return 1;
  else
    return 0;
}


int sci_mlist_check (int number)
{
	if (GetType(number)==sci_mlist)
		return 1;
	else
		return 0;
}

//int sci_mlist_real (int number)
//{
//  int typ1, typ2;
//  typ1 = sci_mlist_check (number);
//  real_or_complex(number, &typ2);
//  if (typ1 && (!typ2))
//    return 1;
//  else
//    return 0;
//}

int sci_matrix_scalar_real (int number)
{
  int typ1, typ2;
  scalar_check(number, &typ1);
  real_or_complex(number, &typ2);
  if (typ1 && (!typ2) && (GetType(number)==sci_matrix))
    return 1;
  else
    return 0;
}

int sci_matrix_void (int number)
{
  int type;
  void_check(number, &type);
  if (type && (GetType(number)==sci_matrix))
    return 1;
  else
    return 0;
}

int sci_matrix_matrix_real (int number)
{
  int typ1, typ2;
  matrix_check(number, &typ1);
  real_or_complex(number, &typ2);
  if (typ1 && (!typ2) && (GetType(number)==sci_matrix))
    return 1;
  else
    return 0;
}

int sci_strings_scalar (int number)
{
  int typ1;
  scalar_check(number, &typ1);
  if (typ1 && (GetType(number)==sci_strings))
    return 1;
  else
    return 0;
}

int sci_strings_vector (int number)
{
  int typ1;
  vector_check(number, &typ1);
  if (typ1 && (GetType(number)==sci_strings))
    return 1;
  else
    return 0;
}




int scalar_string_check(char *l, char c)
{
  if (*l==c)
    return 1;
  else 
    return 0;
}

int length_check(int number, int leng)
{
  int row, col;
  GetMatrixdims(number,&row,&col);
  if (row*col == leng)
    return 1;
  else
    return 0;
}

int vector_length_check(int number1, int number2)
{
  int row1, col1, row2, col2;
  GetMatrixdims(number1, &row1, &col1);
  GetMatrixdims(number2, &row2, &col2);
  if ((row1*col1)==(row2*col2))
    return 1;
  else
    return 0;
}

void vector_length_compare(int number, int leng, int *res)
{
  int row, col;
  GetMatrixdims(number,&row,&col);
  if (row*col == leng)
    *res = 0;
  else if (row*col > leng)
    *res = 1;
  else
    *res = -1;
  return;
}

void matrix_length_compare(int number, int rowLeng, int colLeng, 
			   int *resRow, int *resCol)
{
  int row, col;
  GetMatrixdims(number,&row,&col);
  if (row == rowLeng)
    *resRow = 0;
  else if (row > rowLeng)
    *resRow = 1;
  else
    *resRow = -1;
  if (col == colLeng)
    *resCol = 0;
  else if (col > colLeng)
    *resCol = 1;
  else
    *resCol = -1;
  return;
}

int 
matrix_col_length_check(int number, int leng)
{
  int row, col;
  GetMatrixdims(number, &row, &col);
  if (col==leng)
    return 1;
  else 
    return 0;
}

int 
matrix_row_length_check(int number, int leng)
{
  int row, col;
  GetMatrixdims(number, &row, &col);
  if (row==leng)
    return 1;
  else 
    return 0;
}

int
matrix_length_check (int number1, int number2)
{
  int row1, col1, row2, col2;
  GetMatrixdims(number1, &row1, &col1);
  GetMatrixdims(number2, &row2, &col2);
  if ((row1==row2) && (col1==col2))
    return 1;
  else
    return 0;
}

void extension_check(char *mode, int *type)
{
  int count;
  *type = 0;
  for (count=0;count<extensionIdentityNum;count++)
    {
      if (strcmp(mode,ei[count].extMethodName) == 0)
	{
	  *type = 1;
	  break;
	}
    }
  return;
}

void wavelet_family_check(char *wname, int wf, int *type)
{
  int count;
  *type = 0;
  for(count=0;count<waveletIdentityNum;count++)
    {
      if ((strcmp(wname,wi[count].wname) == 0) && 
	  (wi[count].family==wf))
	{
	  *type = 1;
	  break;
	}
    }
  return;
}


/*-------------------------------------------
 * Validation Error Notification
 *-----------------------------------------*/

void
validate_print (int errCode)
{
  int count;
  for(count=0;count<errorNum;count++)
  {
    if(strErrNoti[count].errorNumber==errCode)
	{
		sciprint(strErrNoti[count].message);
		break;
	}
  }
  return;
}
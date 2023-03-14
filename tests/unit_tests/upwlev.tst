// -------------------------------------------------------------------------
// SWT - Scilab wavelet toolbox
// Copyright (C) 2010-2014  Holger Nahrstaedt
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//-------------------------------------------------------------------------
//
//  <-- NO CHECK ERROR OUTPUT -->


// dwt1d  Test

loadmatfile("-mat",fullfile(get_swt_path(),"tests/unit_tests/Data.mat"));

// upwlev

[c,l]=wavedec(s1,3,'sym10');
[cA1,cD1]=dwt(s1,'sym10');
[cA2,cD2]=dwt(cA1,'sym10');
[cA3,cD3]=dwt(cA2,'sym10');
[c2,l2,ca3]=upwlev(c,l,'sym10');
a3=c(1:length(cA3));
aa2=idwt(cA3,cD3,'sym10',length(cA2));
cc2=[aa2,c(length(cA3)+length(cD3)+1:$)];
ll2=[l(3),l(3:$)];

assert_checkalmostequal ( a3 , ca3 , %eps );
assert_checkalmostequal ( cc2 , c2 , %eps );
assert_checkalmostequal ( ll2 , l2 , %eps );

[Lo_R,Hi_R]=wfilters('sym10','r');
[c2,l2,ca3]=upwlev(c,l,Lo_R,Hi_R);
assert_checkalmostequal ( a3 , ca3 , %eps );
assert_checkalmostequal ( cc2 , c2 , %eps );
assert_checkalmostequal ( ll2 , l2 , %eps );


clear x1;
clear x2;
clear s1;
clear d1;
clear d2;
clear Lo_R;
clear Hi_R;

clear c;
clear l;
clear cA1;
clear cA2;
clear cA3;
clear cD1;
clear cD2;
clear cD3;
clear a3;
clear aa2;
clear ll2;
clear l2;
clear cc2;
clear a3;
clear c2;
clear l2;
clear ca3;

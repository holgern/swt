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

// upcoef
 a=rand(1,50,'normal');
 [Lo_R,Hi_R]=wfilters('sym10','r');
x1=idwt(a,[],'sym10');
x2=idwt(a,[],'sym10',25);
x22=idwt(x1,[],'sym10');
x3=idwt(x22,[],'sym10');
x33=idwt(x22,[],'sym10',25);
y1=upcoef('a',a,'sym10');
y2=upcoef('a',a,Lo_R,Hi_R);
y3=upcoef('a',a,'sym10',1,25);
y4=upcoef('a',a,Lo_R,Hi_R,1,25);
y5=upcoef('a',a,'sym10',1);
y6=upcoef('a',a,Lo_R,Hi_R,1);
y7=upcoef('a',a,'sym10',3);
y8=upcoef('a',a,Lo_R,Hi_R,3);
y9=upcoef('a',a,'sym10',3,25);
y10=upcoef('a',a,Lo_R,Hi_R,3,25);

assert_checkalmostequal ( x1 , y1(19:$-18) , %eps );
assert_checkalmostequal ( x1 , y2(19:$-18) , %eps );
assert_checkalmostequal ( x2 , y3 , %eps );
assert_checkalmostequal ( x2 , y4 , %eps );
assert_checkalmostequal ( x1 , y5(19:$-18) , %eps );
assert_checkalmostequal ( x1 , y6(19:$-18) , %eps );
assert_checkalmostequal ( x3 , y7(127:$-126) , %eps );
assert_checkalmostequal ( x3 , y8(127:$-126) , %eps );
assert_checkalmostequal ( x33 , y9 , %eps );
assert_checkalmostequal ( x33 , y10 , %eps );



d=rand(1,50,'normal');
x1=idwt([],d,'sym10');
x2=idwt([],d,'sym10',25);
x22=idwt(x1,[],'sym10');
x3=idwt(x22,[],'sym10');
x33=idwt(x22,[],'sym10',25);
y1=upcoef('d',d,'sym10');
y2=upcoef('d',d,Lo_R,Hi_R);
y3=upcoef('d',d,'sym10',1,25);
y4=upcoef('d',d,Lo_R,Hi_R,1,25);
y5=upcoef('d',d,'sym10',1);
y6=upcoef('d',d,Lo_R,Hi_R,1);
y7=upcoef('d',d,'sym10',3);
y8=upcoef('d',d,Lo_R,Hi_R,3);
y9=upcoef('d',d,'sym10',3,25);
y10=upcoef('d',d,Lo_R,Hi_R,3,25);
assert_checkalmostequal ( x1 , y1(19:$-18) , %eps );
assert_checkalmostequal ( x1 , y2(19:$-18) , %eps );
assert_checkalmostequal ( x2 , y3 , %eps );
assert_checkalmostequal ( x2 , y4 , %eps );
assert_checkalmostequal ( x1 , y5(19:$-18) , %eps );
assert_checkalmostequal ( x1 , y6(19:$-18) , %eps );
assert_checkalmostequal ( x3 , y7(127:$-126) , %eps );
assert_checkalmostequal ( x3 , y8(127:$-126) , %eps );
assert_checkalmostequal ( x33 , y9 , %eps );
assert_checkalmostequal ( x33 , y10 , %eps );


//size of upcoef

y1=upcoef('a',1,'sym8',1);

assert_checkalmostequal ( size(y1) , [1,16] , %eps );


clear a;
clear d;
clear x1;
clear x2;
clear x3;
clear x22;
clear x33;
clear y1;
clear y2;
clear y3;
clear y4;
clear y5;
clear y6;
clear y7;
clear y8;
clear y9;
clear y10;

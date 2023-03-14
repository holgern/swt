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
// wrcoef

[c,l]=wavedec(s1,3,'sym10');
a3=wrcoef('a',c,l,'sym10');
aa3=wrcoef('a',c,l,'sym10',3);
a2=wrcoef('a',c,l,'sym10',2);
a1=wrcoef('a',c,l,'sym10',1);
d3=wrcoef('d',c,l,'sym10');
dd3=wrcoef('d',c,l,'sym10',3);
d2=wrcoef('d',c,l,'sym10',2);
d1=wrcoef('d',c,l,'sym10',1);
[cA1,cD1]=dwt(s1,'sym10');
[cA2,cD2]=dwt(cA1,'sym10');
[cA3,cD3]=dwt(cA2,'sym10');
A3=waverec([cA3,zeros(1,length(cA3)+length(cD2)+length(cD1))],l,'sym10');
A2=waverec([cA3,cD3,zeros(1,length(cD2)+length(cD1))],l,'sym10');
A1=waverec([cA3,cD3,cD2,zeros(1,length(cD1))],l,'sym10');
D3=waverec([zeros(1,length(cA3)),cD3,zeros(1,length(cD2)+length(cD1))],l,'sym10');
D2=waverec([zeros(1,length(cA3)+length(cD3)),cD2,zeros(1,length(cD1))],l,'sym10');
D1=waverec([zeros(1,length(cA3)+length(cD3)+length(cD2)),cD1],l,'sym10');
assert_checkalmostequal ( a3 , A3 , %eps );
assert_checkalmostequal ( aa3 , A3 , %eps );
assert_checkalmostequal ( a2 , A2 , %eps );
assert_checkalmostequal ( a1 , A1 , %eps );
assert_checkalmostequal ( d3 , D3 , %eps );
assert_checkalmostequal ( dd3 , D3 , %eps );
assert_checkalmostequal ( d2 , D2 , %eps );
assert_checkalmostequal ( d1 , D1 , %eps );


[Lo_R,Hi_R]=wfilters('sym10','r');
a3=wrcoef('a',c,l,Lo_R,Hi_R);
aa3=wrcoef('a',c,l,Lo_R,Hi_R,3);
a2=wrcoef('a',c,l,Lo_R,Hi_R,2);
a1=wrcoef('a',c,l,Lo_R,Hi_R,1);
d3=wrcoef('d',c,l,Lo_R,Hi_R);
dd3=wrcoef('d',c,l,Lo_R,Hi_R,3);
d2=wrcoef('d',c,l,Lo_R,Hi_R,2);
d1=wrcoef('d',c,l,Lo_R,Hi_R,1);
[cA1,cD1]=dwt(s1,'sym10');
[cA2,cD2]=dwt(cA1,'sym10');
[cA3,cD3]=dwt(cA2,'sym10');
A3=waverec([cA3,zeros(1,length(cA3)+length(cD2)+length(cD1))],l,'sym10');
A2=waverec([cA3,cD3,zeros(1,length(cD2)+length(cD1))],l,'sym10');
A1=waverec([cA3,cD3,cD2,zeros(1,length(cD1))],l,'sym10');
D3=waverec([zeros(1,length(cA3)),cD3,zeros(1,length(cD2)+length(cD1))],l,'sym10');
D2=waverec([zeros(1,length(cA3)+length(cD3)),cD2,zeros(1,length(cD1))],l,'sym10');
D1=waverec([zeros(1,length(cA3)+length(cD3)+length(cD2)),cD1],l,'sym10');
assert_checkalmostequal ( a3 , A3 , %eps );
assert_checkalmostequal ( aa3 , A3 , %eps );
assert_checkalmostequal ( a2 , A2 , %eps );
assert_checkalmostequal ( a1 , A1 , %eps );
assert_checkalmostequal ( d3 , D3 , %eps );
assert_checkalmostequal ( dd3 , D3 , %eps );
assert_checkalmostequal ( d2 , D2 , %eps );
assert_checkalmostequal ( d1 , D1 , %eps );

clear cA1;
clear cA2;
clear cA3;
clear cD1;
clear cD2;
clear cD3;
clear A1;
clear A2;
clear A3;
clear D1;
clear D2;
clear D3;
clear a3;
clear aa3;
clear a2;
clear a1;
clear dd;
clear d3;
clear d2;
clear d1;
clear c;
clear l;
clear x1;
clear x2;
clear s1;
clear d1;
clear d2;

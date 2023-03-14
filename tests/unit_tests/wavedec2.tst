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


// dwt2d  Test
version = getversion("scilab");
if (version(1)<6) then
sz=stacksize();
stacksize(1e7);
end;
loadmatfile("-mat",fullfile(get_swt_path(),"tests/unit_tests/Data.mat"));
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;

// wavedec2
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('haar','d');
[c1,s1]=wavedec2(a,3,'haar');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'haar');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'haar');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'haar');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];

assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );

clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// db1
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db1','d');
[c1,s1]=wavedec2(a,3,'db1');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db1');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db1');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db1');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// db2
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db2','d');
[c1,s1]=wavedec2(a,3,'db2');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db2');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db2');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db2');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// db3
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db3','d');
[c1,s1]=wavedec2(a,3,'db3');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db3');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db3');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db3');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// db4
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db4','d');
[c1,s1]=wavedec2(a,3,'db4');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db4');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db4');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db4');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// db5
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db5','d');
[c1,s1]=wavedec2(a,3,'db5');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db5');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db5');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db5');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// db6
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db6','d');
[c1,s1]=wavedec2(a,3,'db6');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db6');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db6');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db6');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// db7
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db7','d');
[c1,s1]=wavedec2(a,3,'db7');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db7');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db7');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db7');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// db8
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db8','d');
[c1,s1]=wavedec2(a,3,'db8');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db8');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db8');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db8');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// db9
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db9','d');
[c1,s1]=wavedec2(a,3,'db9');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db9');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db9');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db9');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// db10
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('db10','d');
[c1,s1]=wavedec2(a,3,'db10');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'db10');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'db10');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'db10');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// coif1
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('coif1','d');
[c1,s1]=wavedec2(a,3,'coif1');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'coif1');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'coif1');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'coif1');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// coif2
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('coif2','d');
[c1,s1]=wavedec2(a,3,'coif2');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'coif2');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'coif2');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'coif2');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// coif3
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('coif3','d');
[c1,s1]=wavedec2(a,3,'coif3');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'coif3');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'coif3');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'coif3');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// coif4
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('coif4','d');
[c1,s1]=wavedec2(a,3,'coif4');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'coif4');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'coif4');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'coif4');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// coif5
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('coif5','d');
[c1,s1]=wavedec2(a,3,'coif5');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'coif5');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'coif5');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'coif5');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// sym4
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('sym4','d');
[c1,s1]=wavedec2(a,3,'sym4');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'sym4');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'sym4');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'sym4');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// sym5
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('sym5','d');
[c1,s1]=wavedec2(a,3,'sym5');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'sym5');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'sym5');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'sym5');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// sym6
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('sym6','d');
[c1,s1]=wavedec2(a,3,'sym6');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'sym6');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'sym6');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'sym6');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// sym7
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('sym7','d');
[c1,s1]=wavedec2(a,3,'sym7');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'sym7');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'sym7');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'sym7');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;

// sym8
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('sym8','d');
[c1,s1]=wavedec2(a,3,'sym8');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'sym8');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'sym8');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'sym8');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// sym9
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('sym9','d');
[c1,s1]=wavedec2(a,3,'sym9');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'sym9');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'sym9');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'sym9');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;

// sym10
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('sym10','d');
[c1,s1]=wavedec2(a,3,'sym10');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'sym10');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'sym10');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'sym10');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// bior1.1
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior1.1','d');
[c1,s1]=wavedec2(a,3,'bior1.1');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior1.1');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior1.1');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior1.1');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// bior1.3
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior1.3','d');
[c1,s1]=wavedec2(a,3,'bior1.3');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior1.3');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior1.3');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior1.3');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// bior1.5
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior1.5','d');
[c1,s1]=wavedec2(a,3,'bior1.5');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior1.5');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior1.5');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior1.5');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// bior2.2
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior2.2','d');
[c1,s1]=wavedec2(a,3,'bior2.2');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior2.2');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior2.2');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior2.2');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// bior2.4
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior2.4','d');
[c1,s1]=wavedec2(a,3,'bior2.4');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior2.4');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior2.4');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior2.4');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// bior2.6
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior2.6','d');
[c1,s1]=wavedec2(a,3,'bior2.6');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior2.6');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior2.6');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior2.6');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// bior2.8
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior2.8','d');
[c1,s1]=wavedec2(a,3,'bior2.8');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior2.8');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior2.8');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior2.8');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;


// bior3.1
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior3.1','d');
[c1,s1]=wavedec2(a,3,'bior3.1');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior3.1');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior3.1');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior3.1');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// bior3.3
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior3.3','d');
[c1,s1]=wavedec2(a,3,'bior3.3');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior3.3');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior3.3');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior3.3');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;

// bior3.5
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior3.5','d');
[c1,s1]=wavedec2(a,3,'bior3.5');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior3.5');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior3.5');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior3.5');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// bior3.7
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior3.7','d');
[c1,s1]=wavedec2(a,3,'bior3.7');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior3.7');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior3.7');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior3.7');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;

// Bior3.9
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[c1,s1]=wavedec2(a,3,'bior3.9');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior3.9');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior3.9');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior3.9');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );


// Bior4.4
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior4.4','d');
[c1,s1]=wavedec2(a,3,'bior4.4');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior4.4');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior4.4');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior4.4');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;



// per
dwtmode('per','nodisp');
a=rand(500,501,'normal');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[c1,s1]=wavedec2(a,3,'bior3.9');
[c2,s2]=wavedec2(a,3,Lo_D,Hi_D);
[cA1,cH1,cV1,cD1]=dwt2(a,'bior3.9');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'bior3.9');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'bior3.9');
cc=[matrix(cA3,1,length(cA3)) matrix(cH3,1,length(cH3)) matrix(cV3,1,length(cV3)) matrix(cD3,1,length(cD3))];
cc=[cc matrix(cH2,1,length(cH2)) matrix(cV2,1,length(cV2)) matrix(cD2,1,length(cD2))];
cc=[cc matrix(cH1,1,length(cH1)) matrix(cV1,1,length(cV1)) matrix(cD1,1,length(cD1))];
ss=[size(cA3);size(cH3);size(cH2);size(cH1);size(a)];
assert_checkalmostequal ( c1 , cc , %eps );
assert_checkalmostequal ( c2 , cc , %eps );
clear a;
clear Lo_D;
clear Hi_D;
clear c1;
clear c2;
clear s1;
clear s2;
clear cA1;
clear cH1;
clear cV1;
clear cD1;
clear cA2;
clear cH2;
clear cV2;
clear cD2;
clear cA3;
clear cH3;
clear cV3;
clear cD3;
clear cc;
clear ss;

dwtmode('symh','nodisp');
if (version(1)<6) then
stacksize(sz(1));
clear sz;
end;

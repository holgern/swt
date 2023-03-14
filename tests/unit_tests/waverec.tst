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

// waverec
//haar
[c,l]=wavedec(s1,3,'haar');
[cA1,cD1]=dwt(s1,'haar');
[cA2,cD2]=dwt(cA1,'haar');
[cA3,cD3]=dwt(cA2,'haar');
ca2=idwt(cA3,cD3,'haar',length(cA2));
ca1=idwt(ca2,cD2,'haar',length(cA1));
a0=idwt(ca1,cD1,'haar',length(s1));
s0=waverec(c,l,'haar');

assert_checkalmostequal ( a0 , s0 , %eps );


// db1
for N=1:36
  wname="db"+sprintf("%d",N);

  [c,l]=wavedec(s1,3,wname);
  [cA1,cD1]=dwt(s1,wname);
  [cA2,cD2]=dwt(cA1,wname);
  [cA3,cD3]=dwt(cA2,wname);
  ca2=idwt(cA3,cD3,wname,length(cA2));
  ca1=idwt(ca2,cD2,wname,length(cA1));
  a0=idwt(ca1,cD1,wname,length(s1));
  s0=waverec(c,l,wname);
  assert_checkalmostequal ( a0 , s0 , %eps );
end;


// coif1
for N=1:17
  wname="coif"+sprintf("%d",N);
  [c,l]=wavedec(s1,3,wname);
  [cA1,cD1]=dwt(s1,wname);
  [cA2,cD2]=dwt(cA1,wname);
  [cA3,cD3]=dwt(cA2,wname);
  ca2=idwt(cA3,cD3,wname,length(cA2));
  ca1=idwt(ca2,cD2,wname,length(cA1));
  a0=idwt(ca1,cD1,wname,length(s1));
  s0=waverec(c,l,wname);
  assert_checkalmostequal ( a0 , s0 , %eps );
end;


// sym4
for N=2:20
wname="sym"+sprintf("%d",N);


[c,l]=wavedec(s1,3,wname);
[cA1,cD1]=dwt(s1,wname);
[cA2,cD2]=dwt(cA1,wname);
[cA3,cD3]=dwt(cA2,wname);
ca2=idwt(cA3,cD3,wname,length(cA2));
ca1=idwt(ca2,cD2,wname,length(cA1));
a0=idwt(ca1,cD1,wname,length(s1));
s0=waverec(c,l,wname);
assert_checkalmostequal ( a0 , s0 , %eps );
end;

bior_fam=["bior1.1","bior1.3", "bior1.5", "bior2.2", "bior2.4", "bior2.6",...
"bior2.8", "bior3.1", "bior3.3", "bior3.5", "bior3.7",...
"bior3.9", "bior4.4", "bior5.5", "bior6.8"];
// bior1.1
for i=1:max(size(bior_fam))
[c,l]=wavedec(s1,3,bior_fam(i));
[cA1,cD1]=dwt(s1,bior_fam(i));
[cA2,cD2]=dwt(cA1,bior_fam(i));
[cA3,cD3]=dwt(cA2,bior_fam(i));
ca2=idwt(cA3,cD3,bior_fam(i),length(cA2));
ca1=idwt(ca2,cD2,bior_fam(i),length(cA1));
a0=idwt(ca1,cD1,bior_fam(i),length(s1));
s0=waverec(c,l,bior_fam(i));
assert_checkalmostequal ( a0 , s0 , %eps );
end;


// type 2
for i=1:max(size(bior_fam))
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(bior_fam(i));
[c,l]=wavedec(s1,3,bior_fam(i));
s0=waverec(c,l,Lo_R,Hi_R);
[cA1,cD1]=dwt(s1,bior_fam(i));
[cA2,cD2]=dwt(cA1,bior_fam(i));
[cA3,cD3]=dwt(cA2,bior_fam(i));
ca2=idwt(cA3,cD3,bior_fam(i),length(cA2));
ca1=idwt(ca2,cD2,bior_fam(i),length(cA1));
a0=idwt(ca1,cD1,bior_fam(i),length(s1));
assert_checkalmostequal ( a0 , s0 , %eps );
Lo_R=rand(1,length(Lo_D),'normal');
Hi_R=rand(1,length(Lo_D),'normal');
ca2=idwt(cA3,cD3,Lo_R,Hi_R,length(cA2));
ca1=idwt(ca2,cD2,Lo_R,Hi_R,length(cA1));
a0=idwt(ca1,cD1,Lo_R,Hi_R,length(s1));
s0=waverec(c,l,Lo_R,Hi_R);
assert_checkalmostequal ( a0 , s0 , %eps );
end;



// type 3
for i=1:max(size(bior_fam))
dwtmode("symh",'nodisp');
[c,l]=wavedec(s1,3,bior_fam(i));
a0=waverec(c,l,bior_fam(i));
dwtmode("symw",'nodisp');
a1=waverec(c,l,bior_fam(i));
dwtmode("asymh",'nodisp');
a2=waverec(c,l,bior_fam(i));
dwtmode("asymw",'nodisp');
a3=waverec(c,l,bior_fam(i));
dwtmode("zpd",'nodisp');
a4=waverec(c,l,bior_fam(i));
dwtmode("sp0",'nodisp');
a5=waverec(c,l,bior_fam(i));
dwtmode("sp1",'nodisp');
a6=waverec(c,l,bior_fam(i));
dwtmode("ppd",'nodisp');
a7=waverec(c,l,bior_fam(i));
dwtmode("per",'nodisp');
a8=waverec(c,l,bior_fam(i));
dwtmode("symh",'nodisp');

assert_checkalmostequal ( a0 , a1 , %eps );
assert_checkalmostequal ( a0 , a2 , %eps );
assert_checkalmostequal ( a0 , a3 , %eps );
assert_checkalmostequal ( a0 , a4 , %eps );
assert_checkalmostequal ( a0 , a5 , %eps );
assert_checkalmostequal ( a0 , a6 , %eps );
assert_checkalmostequal ( a0 , a7 , %eps );
assert_checkalmostequal ( a0 , a8 , %eps );
end;

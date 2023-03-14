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

// dwt
s1=s1(:)';
// type 1 input
// haar
[cA,cD]=dwt(x1,'haar');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('haar');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA ,  %eps, %eps*10);
assert_checkalmostequal ( cdd , cD ,  %eps, %eps*10);

[cA,cD]=dwt(x2,'haar');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA ,  %eps, %eps*10);
assert_checkalmostequal ( cdd , cD ,  %eps, %eps*10);
[cA,cD]=dwt(s1,'haar');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA ,  %eps, %eps*10);
assert_checkalmostequal ( cdd , cD ,  %eps, %eps*10);

// db1 - db36
for N=1:30
  wname="db"+sprintf("%d",N);
  [cA,cD]=dwt(x1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10 );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10 );
  [cA,cD]=dwt(x2,wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10 );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10 );
  [cA,cD]=dwt(s1,wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10 );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
end;


// coif1 - coif17
for N=1:10
  wname="coif"+sprintf("%d",N);
  [cA,cD]=dwt(x1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10  );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10  );
  [cA,cD]=dwt(x2,wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10   );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10  );
  [cA,cD]=dwt(s1,wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10  );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10  );
end;



// sym4 - sym10
for N=4:10
  wname="sym"+sprintf("%d",N);
  [cA,cD]=dwt(x1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10 );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10);;
  [cA,cD]=dwt(x2,wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10 );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10 );
  [cA,cD]=dwt(s1,wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10);;
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10 );
end;

// bior1.1
wnames=['bior1.1','bior1.3','bior1.5','bior2.2','bior2.4','bior2.6','bior2.8','bior3.1','bior3.3','bior3.5','bior3.7','bior3.9'];
for N=1:12
  wname=wnames(N);

  [cA,cD]=dwt(x1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10 );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10 );
  [cA,cD]=dwt(x2,wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10 );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10 );
  [cA,cD]=dwt(s1,wname);
  caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
  cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
  assert_checkalmostequal ( caa , cA , %eps, %eps*10 );
  assert_checkalmostequal ( cdd , cD , %eps, %eps*10 );
end;


// type 2 input
Lo_D=rand(1,20,'normal');
Hi_D=rand(1,20,'normal');
[cA,cD]=dwt(x1,Lo_D,Hi_D);
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);;
assert_checkalmostequal ( cdd , cD , %eps, %eps*100);;
[cA,cD]=dwt(x2,Lo_D,Hi_D);
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);;
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);;
[cA,cD]=dwt(x2,Lo_D,Hi_D);
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);;
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);;


// type 3 input
// symw
// haar
[cA,cD]=dwt(x1,'db10','mode','symw');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'symw',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symw',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(x2,'db10','mode','symw');
caa=dyaddown(wkeep(conv(wextend(1,'symw',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symw',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(s1,'db10','mode','symw');
caa=dyaddown(wkeep(conv(wextend(1,'symw',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symw',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);


// asymh
[cA,cD]=dwt(x1,'db10','mode','asymh');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'asymh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(x2,'db10','mode','asymh');
caa=dyaddown(wkeep(conv(wextend(1,'asymh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(s1,'db10','mode','asymh');
caa=dyaddown(wkeep(conv(wextend(1,'asymh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);


// asymw
[cA,cD]=dwt(x1,'db10','mode','asymw');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'asymw',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymw',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(x2,'db10','mode','asymw');
caa=dyaddown(wkeep(conv(wextend(1,'asymw',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymw',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(s1,'db10','mode','asymw');
caa=dyaddown(wkeep(conv(wextend(1,'asymw',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymw',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);

// zpd
[cA,cD]=dwt(x1,'db10','mode','zpd');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'zpd',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'zpd',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(x2,'db10','mode','zpd');
caa=dyaddown(wkeep(conv(wextend(1,'zpd',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'zpd',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(s1,'db10','mode','zpd');
caa=dyaddown(wkeep(conv(wextend(1,'zpd',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'zpd',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);

// sp0
[cA,cD]=dwt(x1,'db10','mode','sp0');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'sp0',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp0',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(x2,'db10','mode','sp0');
caa=dyaddown(wkeep(conv(wextend(1,'sp0',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp0',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(s1,'db10','mode','sp0');
caa=dyaddown(wkeep(conv(wextend(1,'sp0',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp0',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);


// sp1
[cA,cD]=dwt(x1,'db10','mode','sp1');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'sp1',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp1',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*100);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(x2,'db10','mode','sp1');
caa=dyaddown(wkeep(conv(wextend(1,'sp1',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp1',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(s1,'db10','mode','sp1');
caa=dyaddown(wkeep(conv(wextend(1,'sp1',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp1',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);


// ppd
[cA,cD]=dwt(x1,'db10','mode','ppd');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'ppd',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'ppd',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(x2,'db10','mode','ppd');
caa=dyaddown(wkeep(conv(wextend(1,'ppd',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'ppd',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(s1,'db10','mode','ppd');
caa=dyaddown(wkeep(conv(wextend(1,'ppd',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'ppd',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);


// per
[cA,cD]=dwt(x1,'db10','mode','per');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'per',x1,length(Lo_D),'b'),Lo_D),(length(x1))));
cdd=dyaddown(wkeep(conv(wextend(1,'per',x1,length(Lo_D),'b'),Hi_D),(length(x1))));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(x2,'db10','mode','per');
caa=dyaddown(wkeep(conv(wextend(1,'per',x2,length(Lo_D),'b'),Lo_D),(length(x2))));
cdd=dyaddown(wkeep(conv(wextend(1,'per',x2,length(Lo_D),'b'),Hi_D),(length(x2))));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);
[cA,cD]=dwt(s1,'db10','mode','per');
caa=dyaddown(wkeep(conv(wextend(1,'per',s1,length(Lo_D),'b'),Lo_D),(length(s1))));
cdd=dyaddown(wkeep(conv(wextend(1,'per',s1,length(Lo_D),'b'),Hi_D),(length(s1))));
assert_checkalmostequal ( caa , cA , %eps, %eps*10);
assert_checkalmostequal ( cdd , cD , %eps, %eps*10);

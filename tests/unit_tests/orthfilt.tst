// Copyright (C) 2010 - H. Nahrstaedt
//
// filter Test
// orthfilt


a=rand(1,50,'normal');
lo_r=a;
lo_d=wrev(lo_r);
hi_r=qmf(lo_r);
hi_d=wrev(hi_r);
[Lo_D,Hi_D,Lo_R,Hi_R]=orthfilt(a);

assert_checkalmostequal ( Hi_R , hi_r , %eps );
assert_checkalmostequal ( Lo_D , lo_d , %eps );
assert_checkalmostequal ( Hi_D , hi_d , %eps );
assert_checkalmostequal ( Lo_R , lo_r , %eps );
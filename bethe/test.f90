  Program test
  complex*16:: a(2,2),ai(2,2)
  a(1,1) = dcmplx(2.d0,0.d0)
  a(1,2) = dcmplx(0.d0,0.d0)
  a(2,2) = dcmplx(3.d0,0.d0)
  a(2,1) = dcmplx(0.d0,0.d0)
  call cmatinv(a,ai,2)
  write(6,*) ai(1,1), ai(1,2)
  write(6,*) ai(2,1), ai(2,2)
  
  end

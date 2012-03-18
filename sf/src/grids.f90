function linspace(start,stop,num,endpoint,mesh) result(array)
  real(8)          :: start,stop,step,array(num)
  integer          :: num,i
  logical,optional :: endpoint
  logical          :: endpoint_
  real(8),optional :: mesh
  if(num<0)call error("linspace: N<0, abort.")
  endpoint_=.true.;if(present(endpoint))endpoint_=endpoint
  if(endpoint_)then
     if(num==1)array(num)=start
     step = (stop-start)/real(num-1,8)
  else
     step = (stop-start)/real(num,8)
  end if
  forall(i=1:num)array(i)=start + real(i-1,8)*step
  if(present(mesh))mesh=step
end function linspace


function logspace(start,stop,num,endpoint,base) result(array)
  real(8)          :: start,stop,array(num)
  integer          :: num,i
  logical,optional :: endpoint
  logical          :: endpoint_
  real(8),optional :: base
  real(8)          :: base_
  real(8)          :: A,B
  base_= 10.d0;if(present(base))base_=base
  if(num<0)call error("logspace: N<0, abort.")
  A=start;if(start==0.d0)A=1.d-5
  B=stop;if(stop==0.d0)B=1.d-5
  A=log(A)/log(base) ; B=log(B)/log(base)
  array = linspace(A,B,num=num,endpoint=endpoint)
  array = base**array
end function logspace



function arange(start,num,endpoint) result(array)
  integer          :: start,array(num)
  integer          :: num,i
  logical,optional :: endpoint
  logical          :: endpoint_
  if(num<0)call error("numspace: N<0, abort.")
  endpoint_=.true.;if(present(endpoint))endpoint_=endpoint
  if(endpoint_)then
     forall(i=1:num)array(i)=start+i-1
  else
     forall(i=1:num-1)array(i)=start+i-1
  end if
end function arange

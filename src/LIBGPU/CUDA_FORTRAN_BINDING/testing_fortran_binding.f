      program cuda_fortran_run_module
      integer data1(128,16), data2(16), res(128,16)

*     Initialize the GPU, use the 1st (0)
      call cuInit(0)
      call cuDeviceGet(idev, 0)
      call cuCtxCreate(ictx, 0, idev)

*     Load module and function
      call cuModuleLoad(imod,
     .     '/data/fortran_cuda/test2.cubin')
      call cuModuleGetFunction(ifunc, imod, 'scale')

*     Fill arrays
      do i=1,16
         do j=1,128
            data1(j,i) = 128*i + j
         enddo
      enddo
      do i=1,16
         data2(i) = 16-i 
      enddo

*     Allocate device pointers
      call cuMemAlloc(idata1_ptr, 16*128*4)
      call cuMemAlloc(idata2_ptr, 16*4)

*     Copy data to device
      call cuMemcpyHtoD(idata1_ptr, data1, 16*128*4)
      call cuMemcpyHtoD(idata2_ptr, data2, 16*4)

*     Set function parameters
*     Function parameters size is 16, because we use 2 pointers in 64 bit
*     each being 8 bytes. For 32 platform, this will be 8.
      call cuParamSeti(ifunc, 0, idata1_ptr)
      call cuParamSeti(ifunc, 8, idata2_ptr)
      call cuParamSetSize(ifunc, 16)

*     Launch the calculation on the GPU, use the 'y' axis of the block
      call cuFuncSetBlockShape(ifunc, 128, 1, 1)
      call cuLaunchGrid(ifunc, 1, 16, 1)
*     Utility function to check what result returned the last call to
*     a CUDA driver function
*     call GetLastCUDAResult(ires)
*     print *,ires

*     Copy the results back
      call cuMemcpyDtoH(res, idata1_ptr, 16*128*4)

*     Release all resources
      call cuMemFree(idata1_ptr)
      call cuMemFree(idata2_ptr)

*     Verify results
      call verify_data(data1, data2, res)

      return
      end program

* Verify the results from the GPU
      subroutine verify_data(data1, data2, res)
      integer data1(128,16), data2(16), res(128,16)
      integer test
      
      test = 1
      do i=1,16
         do j=1,128
             if(res(j,i).ne.(data1(j,i)*data2(i)))then
                  print *,res(j,i),(data1(j,i)*data2(i))
                  test=0
             endif
         enddo
      enddo

      if(test.eq.0)then
          print *,'Results do not match'
      else
          print *,'Execution was OK'
      endif

      return
      end

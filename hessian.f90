module hessian
        implicit none

        contains

        subroutine hess(output, func, n, args, param, m,eps)
    !calculates the output matrix of a function using the centred difference
    !algorithm.
    !------------------------------------------------------
    !Parameters:
    !    func : python function with array like return value
    !        function of which the numerical output should be calculated
!
!        eps : float
!            difference parameter
!
!        args : array like
!            independent variables of function
!
!        param : array-like
!            Array of parameters of the function
!    ----------------------------------------------------
!    returns:
!        output : numpy array
!            the output matrix of the function
                implicit none
                integer, intent(in) :: n
                integer, intent(in) :: m
                real(kind=8), dimension(n,n), intent(out) :: output
                real(kind=8), dimension(n), intent(in) :: args
                real(kind=8), dimension(n) :: temp_args_up
                real(kind=8), dimension(n) :: temp_args_up_up
                real(kind=8), dimension(n) :: temp_args_down
                real(kind=8), dimension(n) :: temp_args_up_down
                real(kind=8), dimension(n) :: temp_args_down_up
                real(kind=8), dimension(n) :: temp_args_down_down
                real(kind=8), dimension(m), intent(in) :: param
                real(kind=8) , intent(in) :: eps
                external  func
                real(kind=8) :: func
                integer :: i,j
                do i=1,n
                        do j=1,i

                            temp_args_up_up = args
                            temp_args_down_down = args
                            temp_args_down_up = args
                            temp_args_up_down = args
                            temp_args_up_up(i) = temp_args_up_up(i)+eps*(1.0d0+abs(args(i)))
                            temp_args_up_up(j) = temp_args_up_up(j)+eps*(1.0d0+abs(args(j)))
                            temp_args_down_down(i) = temp_args_down_down(i) - eps*(1.0d0+abs(args(i)))
                            temp_args_down_down(j) = temp_args_down_down(j)-eps*(1.0d0+abs(args(j)))
                            temp_args_up_down(i) = temp_args_up_down(i)+eps*(1.0d0+abs(args(i)))
                            temp_args_up_down(j) = temp_args_up_down(j)-eps*(1.0d0+abs(args(j)))
                            temp_args_down_up(i) = temp_args_down_up(i)-eps*(1.0d0+abs(args(i)))
                            temp_args_down_up(j) = temp_args_down_up(j)+eps*(1.0d0+abs(args(j)))
                            temp_args_up = args
                            temp_args_down =args
                            if (i==j) then

                                temp_args_up(i) = temp_args_up(i)+ eps*(1.0d0+abs(args(i)))
                                temp_args_down(j) = temp_args_down(j) - eps*(1.0d0+abs(args(j)))

                                output(i,j) = -func(temp_args_down_down,n,m, param) - func(temp_args_up_up,n,m, param)&
                                        + 16.0d0*func(temp_args_up,n,m, param) + 16.0d0*func(temp_args_down, n,m,param) -&
                                        30.0d0*func(args, n,m,param)

                                output(i,j)=output(i,j)*((1/eps)**2.0d0/(1.0d0+abs(args(i)))/(1.0d0+abs(args(j)))/12.0d0)

                            else
                                output(i,j) = (func(temp_args_up_up,n,m, param) - func(temp_args_up_down,n,m, param) - &
                                            func(temp_args_down_up,n,m, param) + func(temp_args_down_down,n,m, param))
                                output(i,j)=output(i,j)*((1/eps)**2.0d0/(1.0d0+abs(args(i)))/(1.0d0+abs(args(j)))/4.0d0)
                                end if
                              output(j,i)=output(i,j)
                        end do
                        end do



        end subroutine hess


        end module hessian

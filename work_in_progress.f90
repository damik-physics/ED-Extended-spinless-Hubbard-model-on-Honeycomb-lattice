module work_in_progress 

    implicit none 

    contains 

    ! Function to calculate the energy of a state given its quantum numbers
        ! double precision function energy(k1, k2, filling, t, mass)
        !     integer :: k1, k2
        !     double precision :: filling, t, mass
        !     ! Energy calculation logic here
        !     energy = -t * (cos(2 * pi * k1 / sites) + cos(2 * pi * k2 / sites)) + mass * (filling - 0.5)
    ! end function energy
    function fact(n) result(res) !Calculates the factorial n!

        !----------------------------------!
        !            Factorial: n!         !
        !----------------------------------!

        implicit none

        integer, intent(in) :: n
        
        integer(kind=8) :: i, res

        if(n < 0) error stop 'Factorial is singular for negative integers'

        if(n == 0) then
            res = 1
        else
            res = 1
            do i = 2, n
                res = res * i
            end do
        end if

    end function fact

    function binomial(n, k) result(res) !Calculates the binomial coefficient C(n, k) = n! / (k! * (n-k)!)

        !---------------------------------------!
        !            Binomial coefficients      !
        !---------------------------------------!

        implicit none
        integer, intent(in) :: n, k
        
        double precision :: res

        res = fact(n) / (fact(k) * fact(n - k))

    end function binomial


    function getPosition(b1, b2, n) result(position)
        ! For two bit patterns b1, b2, it calculates how many combinations of smaller sizes come before them (in a combinatorial ordering).
        ! Returns a unique integer position based on these counts.
        ! Uses half of the total size (n/2) as the binomial coefficient parameter, because b1 and b2 represent halves of a system.
        implicit none
        integer, intent(in) :: b1, b2
        integer, intent(in) :: n  ! total size (previously N)

        integer :: i
        integer :: offset_b1, offset_b2
        integer :: position

        offset_b1 = 0
        ! Sum binomial coefficients C(n/2, i) for i=0 to popcnt(b1)-1
        ! Counts how many subsets of size less than popcnt(b1) exist for b1
        do i = 0, popcnt(b1) - 1
            offset_b1 = offset_b1 + binomial(n / 2, i)
        end do

        offset_b2 = 0
        ! Same for b2
        do i = 0, popcnt(b2) - 1
            offset_b2 = offset_b2 + binomial(n / 2, i)
        end do

        ! Combine offsets to get unique position/index (1-based)
        position = offset_b1 + offset_b2 + 1
    end function getPosition

end module work_in_progress
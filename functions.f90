module functions

    implicit none
    
    contains 


    function trim_name(file_name) result(trimmed_name) !Function to remove blank spaces from file names

        implicit none

        character(len=*), intent(in) :: file_name

        integer       :: ls1, ls2, i
        character*200 :: trimmed_name ! Declare a character variable to hold the trimmed name
        
        trimmed_name = ''
        ls1 = len_trim(file_name) ! Get the length of the input file name
        ls2 = 0 ! Initialize the length of the trimmed name

        do i = 1, ls1
            if(file_name(i:i) .ne. ' ') then ! Check if the character is not a blank space
            ls2 = ls2 + 1 
            trimmed_name(ls2:ls2) = file_name(i:i) ! Append non-blank characters to the trimmed name
            end if
        end do

        return
        
    end function trim_name

    function combs(sites, particles, magnetization) result(combinations) 
        !Calculates the number of possible combinations for a given magnetization and particle number

        !-------------------------------------------!
        !            Number of basis states         !
        !-------------------------------------------!

        implicit none

        integer, intent(in)          :: sites, particles
        double precision, intent(in) :: magnetization
        
        integer :: particles_up, particles_dn ! Number of particles with spin up and down
        integer(kind=8) :: combinations ! Number of combinations

        particles_up = int((particles + 2 * magnetization) / 2)
        particles_dn = particles - particles_up

        if(particles_up < 0 .or. particles_dn < 0) then
            combinations = 0
        else
            combinations = int(fact(sites) / (fact(particles_up) * fact(max(1,sites - particles_up))), 8) * int(fact(sites) / (fact(particles_dn) * fact(max(1, sites-particles_dn))), 8)
        end if

    end function combs

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

    function log2(x) result(res) ! Calculates the logarithm base 2 of x

        implicit none

        double precision, intent(in) :: x
        double precision :: res   

        res = log(x) / log(2.)

    end function log2 

end module functions
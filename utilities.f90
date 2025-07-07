module utilities 
    implicit none

    contains

    subroutine cleanup(inner)
        
        use variables

        implicit none
        logical, intent(in) :: inner
        ! This subroutine deallocates all dynamically allocated arrays used in the program.

        if(inner) then ! Deallocate all dynamically allocated arrays of inner loops
            if(allocated(rc))          deallocate(rc)
            if(allocated(ham))         deallocate(ham)
            if(allocated(ham_dc))      deallocate(ham_dc)
            if(allocated(evals))       deallocate(evals)
            if(allocated(ham_dp))      deallocate(ham_dp)
            if(allocated(norm))        deallocate(norm)
            if(allocated(hamDi))       deallocate(hamDi)
            if(allocated(rcOff))       deallocate(rcOff)    
            if(allocated(abasis))      deallocate(abasis)
            if(allocated(bbasis))      deallocate(bbasis)
            if(allocated(hamOff))      deallocate(hamOff)
            if(allocated(mombasis))    deallocate(mombasis)
            if(allocated(parities))    deallocate(parities)
            if(allocated(eigstate))    deallocate(eigstate)
            if(allocated(energies))    deallocate(energies)
            if(allocated(hamOff_dp))   deallocate(hamOff_dp)
            if(allocated(hamOff_dc))   deallocate(hamOff_dc)
            if(allocated(occ))         deallocate(occ)
            if(allocated(dplcts))      deallocate(dplcts)
            if(allocated(eigstate_dc)) deallocate(eigstate_dc)
            if(allocated(basis))       deallocate(basis)
        else ! Deallocate all dynamically allocated arrays of outer loops
            if(allocated(xy))          deallocate(xy)
            if(allocated(xtransl))     deallocate(xtransl)
            if(allocated(ytransl))     deallocate(ytransl)
            if(allocated(alattice))    deallocate(alattice)
            if(allocated(blattice))    deallocate(blattice)
            if(allocated(phases))      deallocate(phases)
            if(allocated(bsites))      deallocate(bsites)
            if(allocated(hexsites))    deallocate(hexsites)
            if(allocated(asitesbonds)) deallocate(asitesbonds)
            if(allocated(bsitesbonds)) deallocate(bsitesbonds)
            if(allocated(latticevecs)) deallocate(latticevecs)
        end if 


    end subroutine cleanup

end module utilities
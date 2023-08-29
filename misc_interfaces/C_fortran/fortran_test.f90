program main

  use iso_c_binding, only : C_CHAR, C_NULL_CHAR
  implicit none

  interface
    subroutine allocate_colvars ( string ) bind ( C, name = "allocate_Colvars" )
      use iso_c_binding, only : C_CHAR
      character ( kind = C_CHAR ) :: string ( * )
    end subroutine allocate_colvars

    subroutine call_colvars () bind ( C, name = "call_proxy_member_static" )
      use iso_c_binding, only : C_CHAR
    end subroutine call_colvars
  end interface

  call allocate_colvars ( C_CHAR_"Hello from Fortran!" // C_NULL_CHAR )
  call call_colvars ()
end
module precision_standard

  integer,parameter :: i1b = selected_int_kind(2)
  integer,parameter :: i2b = selected_int_kind(4)
  integer,parameter :: i4b = selected_int_kind(9)
  integer,parameter :: i8b = selected_int_kind(18)
  integer,parameter :: sp = kind(0.0)
  integer,parameter :: dp = kind(0.0d0)

end module precision_standard

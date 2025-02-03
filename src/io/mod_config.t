!> Module that allows working with a configuration file
module mod_config

  implicit none
  private

  !> The double precision kind-parameter
  integer, parameter :: dp               = kind(0.0d0)

  integer, parameter :: CFG_num_types    = 4 !< Number of variable types
  integer, parameter :: CFG_integer_type = 1 !< Integer type
  integer, parameter :: CFG_real_type    = 2 !< Real number type
  integer, parameter :: CFG_string_type  = 3 !< String type
  integer, parameter :: CFG_logic_type   = 4 !< Boolean/logical type
  integer, parameter :: CFG_unknown_type = 0 !< Used before a variable is created

  integer, parameter :: CFG_name_len   = 80  !< Maximum length of variable names
  integer, parameter :: CFG_string_len = 200 !< Fixed length of string type

  !> Names of the types
  character(len=10), parameter :: CFG_type_names(0:CFG_num_types) = &
       [character(len=10) :: "storage", "integer", "real", "string ", "logical"]

  !> Maximum number of entries in a variable (if it's an array)
  integer, parameter :: CFG_max_array_size = 20

  !> The separator(s) for array-like variables (space, comma, ', ", and tab)
  character(len=*), parameter :: CFG_separators = " ,'"""//char(9)

  !> The separator for categories (stored in var_name)
  character(len=*), parameter :: CFG_category_separator = "%"

  !> The type of a configuration variable
  type CFG_var_t
     private
     !> Type of variable
     integer                       :: var_type
     !> Size of variable, 1 means scalar, > 1 means array
     integer                       :: var_size
     !> Whether the variable size is flexible
     logical                       :: dynamic_size
     !> Whether the variable's value has been requested
     logical                       :: used
     !> Data that has been read in for this variable
     character(len=CFG_string_len) :: stored_data
     !> Name of the variable
     character(len=CFG_name_len)   :: var_name
     !> Description of variable
     character(len=CFG_string_len) :: description

     ! These are the arrays used for storage. In the future, a "pointer" based
     ! approach could be used.
     real(dp), allocatable                      :: real_data(:)
     integer, allocatable                       :: int_data(:)
     character(len=CFG_string_len), allocatable :: char_data(:)
     logical, allocatable                       :: logic_data(:)
  end type CFG_var_t

  !> The configuration that contains all the variables
  type CFG_t
     integer                      :: num_vars = 0
     logical                      :: sorted = .false.
     type(CFG_var_t), allocatable :: vars(:)
  end type CFG_t

  !> Interface to add variables to the configuration
  interface CFG_add
     module procedure :: add_real, add_real_array
     module procedure :: add_int, add_int_array
     module procedure :: add_string, add_string_array
     module procedure :: add_logic, add_logic_array
  end interface CFG_add

  !> Interface to get variables from the configuration
  interface CFG_get
     module procedure :: get_real, get_real_array
     module procedure :: get_int, get_int_array
     module procedure :: get_logic, get_logic_array
     module procedure :: get_string, get_string_array
  end interface CFG_get

  !> Interface to get variables from the configuration
  interface CFG_add_get
     module procedure :: add_get_real, add_get_real_array
     module procedure :: add_get_int, add_get_int_array
     module procedure :: add_get_logic, add_get_logic_array
     module procedure :: add_get_string, add_get_string_array
  end interface CFG_add_get

  ! Public types
  public :: CFG_t
  public :: CFG_integer_type
  public :: CFG_real_type
  public :: CFG_string_type
  public :: CFG_logic_type
  public :: CFG_type_names

  ! Constants
  public :: CFG_name_len
  public :: CFG_string_len
  public :: CFG_max_array_size

  ! Public methods
  public :: CFG_add
  public :: CFG_get
  public :: CFG_add_get
  public :: CFG_get_size
  public :: CFG_get_type
  public :: CFG_check
  public :: CFG_sort
  public :: CFG_write
  public :: CFG_write_markdown
  public :: CFG_read_file
  public :: CFG_update_from_arguments

contains

  subroutine CFG_update_from_arguments(cfg)
    type(CFG_t),intent(inout)      :: cfg
    character(len=100)             :: cfg_name
    integer                        :: ix

    do ix = 1, command_argument_count()
       call get_command_argument(ix, cfg_name)
       call CFG_read_file(cfg, trim(cfg_name))
    end do
  end subroutine CFG_update_from_arguments

  !> This routine will be called if an error occurs in one of the subroutines of
  !> this module.
  subroutine handle_error(err_string)
    character(len=*), intent(in) :: err_string

    print *, "The following error occured in mod_config:"
    print *, trim(err_string)

    ! It is usually best to quit after an error, to make sure the error message
    ! is not overlooked in the program's output
    error stop
  end subroutine handle_error

  !> Return the index of the variable with name 'var_name', or -1 if not found.
  subroutine get_var_index(cfg, var_name, ix)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: var_name
    integer, intent(out)         :: ix
    integer                      :: i

    if (cfg%sorted) then
       call binary_search_variable(cfg, var_name, ix)
    else
       ! Linear search
       do i = 1, cfg%num_vars
          if (cfg%vars(i)%var_name == var_name) exit
       end do

       ! If not found, set i to -1
       if (i == cfg%num_vars + 1) i = -1
       ix = i
    end if

  end subroutine get_var_index

  !> Update the variables in the configartion with the values found in 'filename'
  subroutine CFG_read_file(cfg, filename)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: filename

    integer, parameter            :: my_unit = 123
    integer                       :: io_state, equal_sign_ix
    integer                       :: ix, line_number
    character(len=CFG_name_len)   :: var_name, category
    character(len=CFG_name_len)   :: line_fmt
    character(len=CFG_string_len) :: err_string
    character(len=CFG_string_len) :: line

    open(my_unit, FILE=trim(filename), status = "OLD", &
         action="READ", err=998, iostat=io_state)
    line_number = 0

    write(line_fmt, "(A,I0,A)") "(A", CFG_string_len, ")"

    ! Default category is empty
    category = ""

    do
       read(my_unit, FMT=trim(line_fmt), ERR=998, end=999) line
       line_number = line_number + 1

       call trim_comment(line, '#')

       ! Skip empty lines
       if (line == "") cycle

       ! Locate the '=' sign
       equal_sign_ix = scan(line, '=')

       ! if there is no '='-sign then a category is indicated
       if (equal_sign_ix == 0) then
          line = adjustl(line)

          ! The category name should appear like this: [category_name]
          ix = scan(line, ']')
          if (line(1:1) /= '[' .or. ix == 0) then
             write(err_string, *) "Cannot read line ", line_number, &
                  " from ", trim(filename)
             call handle_error(err_string)
          else
             category = line(2:ix-1)
             cycle
          end if
       end if

       var_name = line(1 : equal_sign_ix - 1) ! Set variable name

       ! If there is no indent, reset to no category
       if (var_name(1:1) /= " " .and. var_name(1:1) /= char(9)) then
          category = ""
       end if

       ! Remove leading blanks
       var_name = adjustl(var_name)

       ! Add category if it is defined
       if (category /= "") then
          var_name = trim(category) // CFG_category_separator // var_name
       end if

       line     = line(equal_sign_ix + 1:)    ! Set line to the values behind the '=' sign
       line     = adjustl(line)               ! Remove leading blanks

       ! Find variable corresponding to name in file
       call get_var_index(cfg, var_name, ix)

       if (ix <= 0) then
          ! Variable still needs to be created, for now store data as a string
          call prepare_store_var(cfg, trim(var_name), CFG_unknown_type, 1, &
               "Not yet created", ix, .false.)
          cfg%vars(ix)%stored_data = line
       else
          cfg%vars(ix)%stored_data = line
          call read_variable(cfg%vars(ix))
       end if
    end do

998 write(err_string, *) "io_state = ", io_state, " while reading from ", &
         trim(filename), " at line ", line_number
    call handle_error("CFG_read_file:" // err_string)

    ! Routine ends here if the end of "filename" is reached
999 close(my_unit, iostat=io_state)

  end subroutine CFG_read_file

  subroutine read_variable(var)
    type(CFG_var_t), intent(inout)            :: var
    integer                                   :: n, n_entries
    integer                                   :: ix_start(CFG_max_array_size)
    integer                                   :: ix_end(CFG_max_array_size)

    ! Get the start and end positions of the line content, and the number of entries
    call get_fields_string(var%stored_data, CFG_separators, &
         CFG_max_array_size, n_entries, ix_start, ix_end)

    if (var%var_size /= n_entries) then
       if (.not. var%dynamic_size) then
          call handle_error("read_variable: variable [" // &
               & trim(var%var_name) // "] has the wrong size")
       else
          var%var_size = n_entries
          call resize_storage(var)
       end if
    end if

    do n = 1, n_entries
       select case (var%var_type)
       case (CFG_integer_type)
          read(var%stored_data(ix_start(n):ix_end(n)), *) var%int_data(n)
       case (CFG_real_type)
          read(var%stored_data(ix_start(n):ix_end(n)), *) var%real_data(n)
       case (CFG_string_type)
          var%char_data(n) = trim(var%stored_data(ix_start(n):ix_end(n)))
       case (CFG_logic_type)
          read(var%stored_data(ix_start(n):ix_end(n)), *) var%logic_data(n)
       end select
    end do
  end subroutine read_variable

  subroutine trim_comment(line, comment_chars)
    character(len=*), intent(inout) :: line
    character(len=*), intent(in)    :: comment_chars
    integer                         :: n
    character                       :: current_char, need_char

    ! Strip comments, but only outside quoted strings (so that var = '#yolo' is
    ! valid when # is a comment char)
    need_char = ""

    do n = 1, len(line)
       current_char = line(n:n)

       if (need_char == "") then
          if (current_char == "'") then
             need_char = "'"    ! Open string
          else if (current_char == '"') then
             need_char = '"'    ! Open string
          else if (index(current_char, comment_chars) /= 0) then
             line = line(1:n-1) ! Trim line up to comment character
             exit
          end if
       else if (current_char == need_char) then
          need_char = ""        ! Close string
       end if

    end do

  end subroutine trim_comment

  subroutine CFG_check(cfg)
    type(CFG_t), intent(in)       :: cfg
    integer                       :: n
    character(len=CFG_string_len) :: err_string

    do n = 1, cfg%num_vars
       if (cfg%vars(n)%var_type == CFG_unknown_type) then
          write(err_string, *) "CFG_check: unknown variable ", &
               trim(cfg%vars(n)%var_name), " in a config file"
          call handle_error(err_string)
       end if
    end do
  end subroutine CFG_check

  !> This routine writes the current configuration to a file with descriptions
  subroutine CFG_write(cfg_in, filename, hide_unused)
    use iso_fortran_env
    type(CFG_t), intent(in)       :: cfg_in
    character(len=*), intent(in)  :: filename
    logical, intent(in), optional :: hide_unused

    type(CFG_t)                   :: cfg
    integer                       :: i, j, io_state, myUnit
    logical                       :: hide_not_used
    character(len=CFG_name_len)   :: name_format, var_name
    character(len=CFG_name_len)   :: category, prev_category
    character(len=CFG_string_len) :: err_string

    hide_not_used = .false.
    if (present(hide_unused)) hide_not_used = hide_unused

    ! Always print a sorted configuration
    cfg = cfg_in
    if (.not. cfg%sorted) call CFG_sort(cfg)

    write(name_format, FMT="(A,I0,A)") "(A,A", CFG_name_len, ",A)"

    if (filename == "stdout") then
       myUnit = output_unit
    else
       myUnit = 333
       open(myUnit, FILE=filename, ACTION="WRITE", ERR=999, IOSTAT=io_state)
    end if

    category      = ""
    prev_category = ""

    do i = 1, cfg%num_vars
       if (.not. cfg%vars(i)%used .and. hide_not_used) cycle
       if (cfg%vars(i)%var_type == CFG_unknown_type) cycle

       ! Write category when it changes
       call split_category(cfg%vars(i), category, var_name)

       if (category /= prev_category .and. category /= '') then
          write(myUnit, ERR=998, FMT="(A)") '[' // trim(category) // ']'
          prev_category = category
       end if

       ! Indent if inside category
       if (category /= "") then
          write(myUnit, ERR=998, FMT="(A,A,A)") "    # ", &
               trim(cfg%vars(i)%description), ":"
          write(myUnit, ADVANCE="NO", ERR=998, FMT="(A)") &
               "    " // trim(var_name) // " ="
       else
          write(myUnit, ERR=998, FMT="(A,A,A)") "# ", &
               trim(cfg%vars(i)%description), ":"
          write(myUnit, ADVANCE="NO", ERR=998, FMT="(A)") &
               trim(var_name) // " ="
       end if

       select case(cfg%vars(i)%var_type)
       case (CFG_integer_type)
          do j = 1, cfg%vars(i)%var_size
             write(myUnit, ADVANCE="NO", ERR=998, FMT="(A,I0)") &
                  " ", cfg%vars(i)%int_data(j)
          end do
       case (CFG_real_type)
          do j = 1, cfg%vars(i)%var_size
             write(myUnit, ADVANCE="NO", ERR=998, FMT="(A,E11.4)") &
                  " ", cfg%vars(i)%real_data(j)
          end do
       case (CFG_string_type)
          do j = 1, cfg%vars(i)%var_size
             write(myUnit, ADVANCE="NO", ERR=998, FMT="(A)") &
                  " '" // trim(cfg%vars(i)%char_data(j)) // "'"
          end do
       case (CFG_logic_type)
          do j = 1, cfg%vars(i)%var_size
             write(myUnit, ADVANCE="NO", ERR=998, FMT="(A,L1)") &
                  " ", cfg%vars(i)%logic_data(j)
          end do
       end select
       write(myUnit, ERR=998, FMT="(A)") ""
       write(myUnit, ERR=998, FMT="(A)") ""
    end do

    if (myUnit /= output_unit) close(myUnit, ERR=999, IOSTAT=io_state)
    call CFG_check(cfg_in)
    return

998 continue
    write(err_string, *) "CFG_write error: io_state = ", io_state, &
         " while writing ", trim(var_name), " to ", filename
    call handle_error(err_string)

999 continue ! If there was an error, the routine will end here
    write(err_string, *) "CFG_write error: io_state = ", io_state, &
         " while writing to ", filename
    call handle_error(err_string)

  end subroutine CFG_write

  !> This routine writes the current configuration to a markdown file
  subroutine CFG_write_markdown(cfg_in, filename, hide_unused)
    use iso_fortran_env
    type(CFG_t), intent(in)       :: cfg_in
    character(len=*), intent(in)  :: filename
    logical, intent(in), optional :: hide_unused

    type(CFG_t)                   :: cfg
    integer                       :: i, j, io_state, myUnit
    logical                       :: hide_not_used
    character(len=CFG_name_len)   :: name_format, var_name
    character(len=CFG_name_len)   :: category, prev_category
    character(len=CFG_string_len) :: err_string

    hide_not_used = .false.
    if (present(hide_unused)) hide_not_used = hide_unused

    ! Always print a sorted configuration
    cfg = cfg_in
    if (.not. cfg%sorted) call CFG_sort(cfg)

    write(name_format, FMT="(A,I0,A)") "(A,A", CFG_name_len, ",A)"

    if (filename == "stdout") then
       myUnit = output_unit
    else
       myUnit = 333
       open(myUnit, FILE=filename, ACTION="WRITE", ERR=999, IOSTAT=io_state)
    end if

    category      = ""
    prev_category = "X"
    write(myUnit, ERR=998, FMT="(A)") "# Configuration file (markdown format)"
    write(myUnit, ERR=998, FMT="(A)") ""

    do i = 1, cfg%num_vars

       if (.not. cfg%vars(i)%used .and. hide_not_used) cycle
       if (cfg%vars(i)%var_type == CFG_unknown_type) cycle

       ! Write category when it changes
       call split_category(cfg%vars(i), category, var_name)

       if (category /= prev_category) then
          if (category == "") category = "No category"
          write(myUnit, ERR=998, FMT="(A)") '## ' // trim(category)
          write(myUnit, ERR=998, FMT="(A)") ""
          prev_category = category
       end if

       write(myUnit, ERR=998, FMT="(A)") "* " // trim(cfg%vars(i)%description)
       write(myUnit, ERR=998, FMT="(A)") ""
       write(myUnit, ADVANCE="NO", ERR=998, FMT="(A)") &
            '        ' // trim(var_name) // " ="

       select case(cfg%vars(i)%var_type)
       case (CFG_integer_type)
          do j = 1, cfg%vars(i)%var_size
             write(myUnit, ADVANCE="NO", ERR=998, FMT="(A,I0)") &
                  " ", cfg%vars(i)%int_data(j)
          end do
       case (CFG_real_type)
          do j = 1, cfg%vars(i)%var_size
             write(myUnit, ADVANCE="NO", ERR=998, FMT="(A,E11.4)") &
                  " ", cfg%vars(i)%real_data(j)
          end do
       case (CFG_string_type)
          do j = 1, cfg%vars(i)%var_size
             write(myUnit, ADVANCE="NO", ERR=998, FMT="(A)") &
                  " '" // trim(cfg%vars(i)%char_data(j)) // "'"
          end do
       case (CFG_logic_type)
          do j = 1, cfg%vars(i)%var_size
             write(myUnit, ADVANCE="NO", ERR=998, FMT="(A,L1)") &
                  " ", cfg%vars(i)%logic_data(j)
          end do
       end select
       write(myUnit, ERR=998, FMT="(A)") ""
       write(myUnit, ERR=998, FMT="(A)") ""
    end do

    if (myUnit /= output_unit) close(myUnit, ERR=999, IOSTAT=io_state)
    call CFG_check(cfg_in)
    return

998 continue
    write(err_string, *) "CFG_write_markdown error: io_state = ", io_state, &
         " while writing ", trim(var_name), " to ", filename
    call handle_error(err_string)

999 continue ! If there was an error, the routine will end here
    write(err_string, *) "CFG_write_markdown error: io_state = ", io_state, &
         " while writing to ", filename
    call handle_error(err_string)

  end subroutine CFG_write_markdown

  subroutine split_category(variable, category, var_name)
    type(CFG_var_t), intent(in)          :: variable
    character(CFG_name_len), intent(out) :: category
    character(CFG_name_len), intent(out) :: var_name
    integer                              :: ix

    ix = index(variable%var_name, CFG_category_separator)

    if (ix == 0) then
       category = ""
       var_name = variable%var_name
    else
       category = variable%var_name(1:ix-1)
       var_name = variable%var_name(ix+1:)
    end if

  end subroutine split_category

  !> Resize the storage size of variable, which can be of type integer, logical,
  !> real or character
  subroutine resize_storage(variable)
    type(CFG_var_t), intent(inout) :: variable

    select case (variable%var_type)
    case (CFG_integer_type)
       deallocate( variable%int_data )
       allocate( variable%int_data(variable%var_size) )
    case (CFG_logic_type)
       deallocate( variable%logic_data )
       allocate( variable%logic_data(variable%var_size) )
    case (CFG_real_type)
       deallocate( variable%real_data )
       allocate( variable%real_data(variable%var_size) )
    case (CFG_string_type)
       deallocate( variable%char_data )
       allocate( variable%char_data(variable%var_size) )
    end select
  end subroutine resize_storage

  !> Helper routine to store variables. This is useful because a lot of the same
  !> code is executed for the different types of variables.
  subroutine prepare_store_var(cfg, var_name, var_type, var_size, &
       description, ix, dynamic_size)
    type(CFG_t), intent(inout)    :: cfg
    character(len=*), intent(in)  :: var_name, description
    integer, intent(in)           :: var_type, var_size
    integer, intent(out)          :: ix !< Index of variable
    logical, intent(in), optional :: dynamic_size

    ! Check if variable already exists
    call get_var_index(cfg, var_name, ix)

    if (ix == -1) then ! Create a new variable
       call ensure_free_storage(cfg)
       cfg%sorted               = .false.
       ix                       = cfg%num_vars + 1
       cfg%num_vars             = cfg%num_vars + 1
       cfg%vars(ix)%used        = .false.
       cfg%vars(ix)%stored_data = ""
    else
       ! Only allowed when the variable is not yet created
       if (cfg%vars(ix)%var_type /= CFG_unknown_type) then
          call handle_error("prepare_store_var: variable [" // &
               & trim(var_name) // "] already exists")
       end if
    end if

    cfg%vars(ix)%var_name    = var_name
    cfg%vars(ix)%description = description
    cfg%vars(ix)%var_type    = var_type
    cfg%vars(ix)%var_size    = var_size

    if (present(dynamic_size)) then
       cfg%vars(ix)%dynamic_size = dynamic_size
    else
       cfg%vars(ix)%dynamic_size = .false.
    end if

    select case (var_type)
    case (CFG_integer_type)
       allocate( cfg%vars(ix)%int_data(var_size) )
    case (CFG_real_type)
       allocate( cfg%vars(ix)%real_data(var_size) )
    case (CFG_string_type)
       allocate( cfg%vars(ix)%char_data(var_size) )
    case (CFG_logic_type)
       allocate( cfg%vars(ix)%logic_data(var_size) )
    end select

  end subroutine prepare_store_var

  !> Helper routine to get variables. This is useful because a lot of the same
  !> code is executed for the different types of variables.
  subroutine prepare_get_var(cfg, var_name, var_type, var_size, ix)
    type(CFG_t), intent(inout)    :: cfg
    character(len=*), intent(in)  :: var_name
    integer, intent(in)           :: var_type, var_size
    integer, intent(out)          :: ix
    character(len=CFG_string_len) :: err_string

    call get_var_index(cfg, var_name, ix)

    if (ix == -1) then
       call handle_error("CFG_get: variable ["//var_name//"] not found")
    else if (cfg%vars(ix)%var_type /= var_type) then
       write(err_string, fmt="(A)") "CFG_get: variable [" &
            // var_name // "] has different type (" // &
            trim(CFG_type_names(cfg%vars(ix)%var_type)) // &
            ") than requested (" // trim(CFG_type_names(var_type)) // ")"
       call handle_error(err_string)
    else if (cfg%vars(ix)%var_size /= var_size) then
       write(err_string, fmt="(A,I0,A,I0,A)") "CFG_get: variable [" &
            // var_name // "] has different size (", cfg%vars(ix)%var_size, &
            ") than requested (", var_size, ")"
       call handle_error(err_string)
    else                        ! All good, variable will be used
       cfg%vars(ix)%used = .true.
    end if
  end subroutine prepare_get_var

  !> Add a configuration variable with a real value
  subroutine add_real(cfg, var_name, real_data, comment)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    real(dp), intent(in)         :: real_data
    integer                      :: ix

    call prepare_store_var(cfg, var_name, CFG_real_type, 1, comment, ix)

    if (cfg%vars(ix)%stored_data /= "") then
       call read_variable(cfg%vars(ix))
    else
       cfg%vars(ix)%real_data(1) = real_data
    end if
  end subroutine add_real

  !> Add a configuration variable with an array of type
  !  real
  subroutine add_real_array(cfg, var_name, real_data, comment, dynamic_size)
    type(CFG_t), intent(inout)    :: cfg
    character(len=*), intent(in)  :: var_name, comment
    real(dp), intent(in)          :: real_data(:)
    logical, intent(in), optional :: dynamic_size
    integer                       :: ix

    call prepare_store_var(cfg, var_name, CFG_real_type, &
         size(real_data), comment, ix, dynamic_size)

    if (cfg%vars(ix)%stored_data /= "") then
       call read_variable(cfg%vars(ix))
    else
       cfg%vars(ix)%real_data = real_data
    end if
  end subroutine add_real_array

  !> Add a configuration variable with an integer value
  subroutine add_int(cfg, var_name, int_data, comment)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    integer, intent(in)          :: int_data
    integer                      :: ix

    call prepare_store_var(cfg, var_name, CFG_integer_type, 1, comment, ix)

    if (cfg%vars(ix)%stored_data /= "") then
       call read_variable(cfg%vars(ix))
    else
       cfg%vars(ix)%int_data(1) = int_data
    end if
  end subroutine add_int

  !> Add a configuration variable with an array of type integer
  subroutine add_int_array(cfg, var_name, int_data, comment, dynamic_size)
    type(CFG_t), intent(inout)    :: cfg
    character(len=*), intent(in)  :: var_name, comment
    integer, intent(in)           :: int_data(:)
    logical, intent(in), optional :: dynamic_size
    integer                       :: ix

    call prepare_store_var(cfg, var_name, CFG_integer_type, &
         size(int_data), comment, ix, dynamic_size)

    if (cfg%vars(ix)%stored_data /= "") then
       call read_variable(cfg%vars(ix))
    else
       cfg%vars(ix)%int_data = int_data
    end if
  end subroutine add_int_array

  !> Add a configuration variable with an character value
  subroutine add_string(cfg, var_name, char_data, comment)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment, char_data
    integer                      :: ix

    call prepare_store_var(cfg, var_name, CFG_string_type, 1, comment, ix)
    if (cfg%vars(ix)%stored_data /= "") then
       call read_variable(cfg%vars(ix))
    else
       cfg%vars(ix)%char_data(1) = char_data
    end if
  end subroutine add_string

  !> Add a configuration variable with an array of type character
  subroutine add_string_array(cfg, var_name, char_data, &
       comment, dynamic_size)
    type(CFG_t), intent(inout)    :: cfg
    character(len=*), intent(in)  :: var_name, comment, char_data(:)
    logical, intent(in), optional :: dynamic_size
    integer                       :: ix

    call prepare_store_var(cfg, var_name, CFG_string_type, &
         size(char_data), comment, ix, dynamic_size)

    if (cfg%vars(ix)%stored_data /= "") then
       call read_variable(cfg%vars(ix))
    else
       cfg%vars(ix)%char_data = char_data
    end if
  end subroutine add_string_array

  !> Add a configuration variable with an logical value
  subroutine add_logic(cfg, var_name, logic_data, comment)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    logical, intent(in)          :: logic_data
    integer                      :: ix

    call prepare_store_var(cfg, var_name, CFG_logic_type, 1, comment, ix)

    if (cfg%vars(ix)%stored_data /= "") then
       call read_variable(cfg%vars(ix))
    else
       cfg%vars(ix)%logic_data(1) = logic_data
    end if
  end subroutine add_logic

  !> Add a configuration variable with an array of type logical
  subroutine add_logic_array(cfg, var_name, logic_data, &
       comment, dynamic_size)
    type(CFG_t), intent(inout)    :: cfg
    character(len=*), intent(in)  :: var_name, comment
    logical, intent(in)           :: logic_data(:)
    logical, intent(in), optional :: dynamic_size
    integer                       :: ix

    call prepare_store_var(cfg, var_name, CFG_logic_type, &
         size(logic_data), comment, ix, dynamic_size)

    if (cfg%vars(ix)%stored_data /= "") then
       call read_variable(cfg%vars(ix))
    else
       cfg%vars(ix)%logic_data = logic_data
    end if
  end subroutine add_logic_array

  !> Get a real array of a given name
  subroutine get_real_array(cfg, var_name, real_data)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name
    real(dp), intent(inout)      :: real_data(:)
    integer                      :: ix

    call prepare_get_var(cfg, var_name, CFG_real_type, &
         size(real_data), ix)
    real_data = cfg%vars(ix)%real_data
  end subroutine get_real_array

  !> Get a integer array of a given name
  subroutine get_int_array(cfg, var_name, int_data)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name
    integer, intent(inout)       :: int_data(:)
    integer                      :: ix

    call prepare_get_var(cfg, var_name, CFG_integer_type, &
         size(int_data), ix)
    int_data    = cfg%vars(ix)%int_data
  end subroutine get_int_array

  !> Get a character array of a given name
  subroutine get_string_array(cfg, var_name, char_data)
    type(CFG_t), intent(inout)      :: cfg
    character(len=*), intent(in)    :: var_name
    character(len=*), intent(inout) :: char_data(:)
    integer                         :: ix

    call prepare_get_var(cfg, var_name, CFG_string_type, &
         size(char_data), ix)
    char_data = cfg%vars(ix)%char_data
  end subroutine get_string_array

  !> Get a logical array of a given name
  subroutine get_logic_array(cfg, var_name, logic_data)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name
    logical, intent(inout)       :: logic_data(:)
    integer                      :: ix

    call prepare_get_var(cfg, var_name, CFG_logic_type, &
         size(logic_data), ix)
    logic_data = cfg%vars(ix)%logic_data
  end subroutine get_logic_array

  !> Get a real value of a given name
  subroutine get_real(cfg, var_name, res)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name
    real(dp), intent(out)        :: res
    integer                      :: ix

    call prepare_get_var(cfg, var_name, CFG_real_type, 1, ix)
    res = cfg%vars(ix)%real_data(1)
  end subroutine get_real

  !> Get a integer value of a given name
  subroutine get_int(cfg, var_name, res)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name
    integer, intent(inout)       :: res
    integer                      :: ix

    call prepare_get_var(cfg, var_name, CFG_integer_type, 1, ix)
    res = cfg%vars(ix)%int_data(1)
  end subroutine get_int

  !> Get a logical value of a given name
  subroutine get_logic(cfg, var_name, res)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name
    logical, intent(out)         :: res
    integer                      :: ix

    call prepare_get_var(cfg, var_name, CFG_logic_type, 1, ix)
    res = cfg%vars(ix)%logic_data(1)
  end subroutine get_logic

  !> Get a character value of a given name
  subroutine get_string(cfg, var_name, res)
    type(CFG_t), intent(inout)    :: cfg
    character(len=*), intent(in)  :: var_name
    character(len=*), intent(out) :: res
    integer                       :: ix

    call prepare_get_var(cfg, var_name, CFG_string_type, 1, ix)
    res = cfg%vars(ix)%char_data(1)
  end subroutine get_string

  !> Get or add a real array of a given name
  subroutine add_get_real_array(cfg, var_name, real_data, &
       comment, dynamic_size)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    real(dp), intent(inout)      :: real_data(:)
    logical, intent(in), optional :: dynamic_size

    call add_real_array(cfg, var_name, real_data, comment, dynamic_size)
    call get_real_array(cfg, var_name, real_data)
  end subroutine add_get_real_array

  !> Get or add a integer array of a given name
  subroutine add_get_int_array(cfg, var_name, int_data, &
       comment, dynamic_size)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    integer, intent(inout)       :: int_data(:)
    logical, intent(in), optional :: dynamic_size

    call add_int_array(cfg, var_name, int_data, comment, dynamic_size)
    call get_int_array(cfg, var_name, int_data)
  end subroutine add_get_int_array

  !> Get or add a character array of a given name
  subroutine add_get_string_array(cfg, var_name, char_data, &
       comment, dynamic_size)
    type(CFG_t), intent(inout)      :: cfg
    character(len=*), intent(in)    :: var_name, comment
    character(len=*), intent(inout) :: char_data(:)
    logical, intent(in), optional :: dynamic_size

    call add_string_array(cfg, var_name, char_data, comment, dynamic_size)
    call get_string_array(cfg, var_name, char_data)
  end subroutine add_get_string_array

  !> Get or add a logical array of a given name
  subroutine add_get_logic_array(cfg, var_name, logic_data, &
       comment, dynamic_size)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    logical, intent(inout)       :: logic_data(:)
    logical, intent(in), optional :: dynamic_size

    call add_logic_array(cfg, var_name, logic_data, comment, dynamic_size)
    call get_logic_array(cfg, var_name, logic_data)
  end subroutine add_get_logic_array

  !> Get or add a real value of a given name
  subroutine add_get_real(cfg, var_name, real_data, comment)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    real(dp), intent(inout)      :: real_data

    call add_real(cfg, var_name, real_data, comment)
    call get_real(cfg, var_name, real_data)
  end subroutine add_get_real

  !> Get or add a integer value of a given name
  subroutine add_get_int(cfg, var_name, int_data, comment)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    integer, intent(inout)       :: int_data

    call add_int(cfg, var_name, int_data, comment)
    call get_int(cfg, var_name, int_data)
  end subroutine add_get_int

  !> Get or add a logical value of a given name
  subroutine add_get_logic(cfg, var_name, logical_data, comment)
    type(CFG_t), intent(inout)   :: cfg
    character(len=*), intent(in) :: var_name, comment
    logical, intent(inout)       :: logical_data

    call add_logic(cfg, var_name, logical_data, comment)
    call get_logic(cfg, var_name, logical_data)
  end subroutine add_get_logic

  !> Get a character value of a given name
  subroutine add_get_string(cfg, var_name, string_data, comment)
    type(CFG_t), intent(inout)       :: cfg
    character(len=*), intent(in)     :: var_name, comment
    character(len=*), intent(inout)  :: string_data

    call add_string(cfg, var_name, string_data, comment)
    call get_string(cfg, var_name, string_data)
  end subroutine add_get_string

  !> Get the size of a variable
  subroutine CFG_get_size(cfg, var_name, res)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: var_name
    integer, intent(out)         :: res
    integer                      :: ix

    call get_var_index(cfg, var_name, ix)
    if (ix /= -1) then
       res = cfg%vars(ix)%var_size
    else
       res = -1
       call handle_error("CFG_get_size: variable ["//var_name//"] not found")
    end if
  end subroutine CFG_get_size

  !> Get the type of a given variable of a configuration type
  subroutine CFG_get_type(cfg, var_name, res)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: var_name
    integer, intent(out)         :: res
    integer                      :: ix

    call get_var_index(cfg, var_name, ix)

    if (ix /= -1) then
       res = cfg%vars(ix)%var_type
    else
       res = -1
       call handle_error("CFG_get_type: variable ["//var_name//"] not found")
    end if
  end subroutine CFG_get_type

  !> Routine to ensure that enough storage is allocated for the configuration
  !> type. If not the new size will be twice as much as the current size. If no
  !> storage is allocated yet a minumum amount of starage is allocated.
  subroutine ensure_free_storage(cfg)
    type(CFG_t), intent(inout)   :: cfg
    type(CFG_var_t), allocatable :: cfg_copy(:)
    integer, parameter           :: min_dyn_size = 100
    integer                      :: cur_size, new_size

    if (allocated(cfg%vars)) then
       cur_size = size(cfg%vars)

       if (cur_size < cfg%num_vars + 1) then
          new_size = 2 * cur_size
          allocate(cfg_copy(cur_size))
          cfg_copy = cfg%vars
          deallocate(cfg%vars)
          allocate(cfg%vars(new_size))
          cfg%vars(1:cur_size) = cfg_copy
       end if
    else
       allocate(cfg%vars(min_dyn_size))
    end if

  end subroutine ensure_free_storage

  !> Routine to find the indices of entries in a string
  subroutine get_fields_string(line, delims, n_max, n_found, ixs_start, ixs_end)
    !> The line from which we want to read
    character(len=*), intent(in)  :: line
    !> A string with delimiters. For example delims = " ,'"""//char(9)
    character(len=*), intent(in)  :: delims
    !> Maximum number of entries to read in
    integer, intent(in)           :: n_max
    !> Number of entries found
    integer, intent(inout)        :: n_found
    !> On return, ix_start(i) holds the starting point of entry i
    integer, intent(inout)        :: ixs_start(n_max)
    !> On return, ix_end(i) holds the end point of entry i
    integer, intent(inout)        :: ixs_end(n_max)

    integer                       :: ix, ix_prev

    ix_prev = 0
    n_found = 0

    do while (n_found < n_max)

       ! Find the starting point of the next entry (a non-delimiter value)
       ix = verify(line(ix_prev+1:), delims)
       if (ix == 0) exit

       n_found            = n_found + 1
       ixs_start(n_found) = ix_prev + ix ! This is the absolute position in 'line'

       ! Get the end point of the current entry (next delimiter index minus one)
       ix = scan(line(ixs_start(n_found)+1:), delims) - 1

       if (ix == -1) then              ! If there is no last delimiter,
          ixs_end(n_found) = len(line) ! the end of the line is the endpoint
       else
          ixs_end(n_found) = ixs_start(n_found) + ix
       end if

       ix_prev = ixs_end(n_found) ! We continue to search from here
    end do

  end subroutine get_fields_string

  !> Performa a binary search for the variable 'var_name'
  subroutine binary_search_variable(cfg, var_name, ix)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: var_name
    integer, intent(out)         :: ix
    integer                      :: i_min, i_max, i_mid

    i_min = 1
    i_max = cfg%num_vars
    ix    = - 1

    do while (i_min < i_max)
       i_mid = i_min + (i_max - i_min) / 2
       if ( llt(cfg%vars(i_mid)%var_name, var_name) ) then
          i_min = i_mid + 1
       else
          i_max = i_mid
       end if
    end do

    ! If not found, binary_search_variable is not set here, and stays -1
    if (i_max == i_min .and. cfg%vars(i_min)%var_name == var_name) then
       ix = i_min
    else
       ix = -1
    end if
  end subroutine binary_search_variable

  !> Sort the variables for faster lookup
  subroutine CFG_sort(cfg)
    type(CFG_t), intent(inout) :: cfg

    call qsort_config(cfg%vars(1:cfg%num_vars))
    cfg%sorted = .true.
  end subroutine CFG_sort

  !> Simple implementation of quicksort algorithm to sort the variable list alphabetically.
  recursive subroutine qsort_config(list)
    type(CFG_var_t), intent(inout) :: list(:)
    integer                        :: split_pos

    if (size(list) > 1) then
       call parition_var_list(list, split_pos)
       call qsort_config( list(:split_pos-1) )
       call qsort_config( list(split_pos:) )
    end if
  end subroutine qsort_config

  !> Helper routine for quicksort, to perform partitioning
  subroutine parition_var_list(list, marker)
    type(CFG_var_t), intent(inout) :: list(:)
    integer, intent(out)           :: marker
    integer                        :: left, right, pivot_ix
    type(CFG_var_t)                :: temp
    character(len=CFG_name_len)    :: pivot_value

    left        = 0
    right       = size(list) + 1

    ! Take the middle element as pivot
    pivot_ix    = size(list) / 2
    pivot_value = list(pivot_ix)%var_name

    do while (left < right)

       right = right - 1
       do while (lgt(list(right)%var_name, pivot_value))
          right = right - 1
       end do

       left = left + 1
       do while (lgt(pivot_value, list(left)%var_name))
          left = left + 1
       end do

       if (left < right) then
          temp = list(left)
          list(left) = list(right)
          list(right) = temp
       end if
    end do

    if (left == right) then
       marker = left + 1
    else
       marker = left
    end if
  end subroutine parition_var_list

end module mod_config

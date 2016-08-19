MODULE SF_IOTOOLS
  USE IOFILE
  USE IOPLOT
  USE IOREAD
  private

  public :: splot
  public :: splot3d
  public :: store_data

  public :: sread
  public :: read_data

  public :: set_store_size
  public :: txtfy
  public :: str
  public :: file_size
  public :: file_length
  public :: file_info
  public :: free_unit
  public :: free_units
  public :: data_open
  public :: data_store
  public :: reg_filename,reg,txtfit,txtcut
  public :: create_data_dir,create_dir
  public :: close_file
  public :: get_filename
  public :: get_filepath

END MODULE SF_IOTOOLS

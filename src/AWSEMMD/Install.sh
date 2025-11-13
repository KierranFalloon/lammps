# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

# enforce using portable C locale
LC_ALL=C
export LC_ALL

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# list of files without dependencies
action  atom_vec_awsemmd.cpp
action  atom_vec_awsemmd.h
action  pair_ex_gauss_coul_cut.cpp
action  pair_ex_gauss_coul_cut.h
action  pair_excluded_volume.cpp
action  pair_excluded_volume.h
action  pair_go-contacts.cpp
action  pair_go-contacts.h
action  fix_backbone.cpp
action  fix_backbone.h
action  fix_go-model.cpp
action  fix_go-model.h
action  fix_print_wzero.cpp
action  fix_print_wzero.h
action  fix_qbias.cpp
action  fix_qbias.h
action  fix_spring_rg_papoian.cpp
action  fix_spring_rg_papoian.h
action  compute_contactmap.cpp
action  compute_contactmap.h
action  compute_pairdistmat.cpp
action  compute_pairdistmat.h
action  compute_q_onuchic.cpp
action  compute_q_onuchic.h
action  compute_q_wolynes.cpp
action  compute_q_wolynes.h
action  compute_totalcontacts.cpp
action  compute_totalcontacts.h
action  fragment_memory.cpp 
action  fragment_memory.h
action  smart_matrix_lib.h

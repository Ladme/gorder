mol new ../files/system.gro
mol addfile ../files/traj_whole.xtc waitfor all

proc orderparam-c3 { result { seltext all } } {
  upvar $result arr
  set n [molinfo top get numframes]
  for { set i 2 } { $i <= 20 } { incr i } {
#    puts $i
    set cp [atomselect top "($seltext) and lipid and name C3$i"]
    if { [$cp num] == 0 } {
#      puts "skipping $i"
      continue
    }
    set hx [atomselect top "($seltext) and lipid and name H${i}X"]
    set hy [atomselect top "($seltext) and lipid and name H${i}Y"]
    set hz [atomselect top "($seltext) and lipid and name H${i}Z"]

    set sum 0.0
    set nh 0
    set nres [$cp num]
    for { set frame 500 } { $frame < $n } { incr frame } {
      $cp frame $frame
      $hx frame $frame
      $hy frame $frame
      $hz frame $frame
      set cpx [$cp get x]
      set cpy [$cp get y]
      set cpz [$cp get z]

      set hxx [vecsub $cpx [$hx get x]]
      set hxy [vecsub $cpy [$hx get y]]
      set hxz [vecsub $cpz [$hx get z]]
      foreach dx $hxx dy $hxy dz $hxz {
        set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
        set sum [expr {$sum + $dz*$dz/$norm2}]
      }
      incr nh $nres

      if { [$hy num] != 0 } {
        set hyx [vecsub $cpx [$hy get x]]
        set hyy [vecsub $cpy [$hy get y]]
        set hyz [vecsub $cpz [$hy get z]]
        foreach dx $hyx dy $hyy dz $hyz {
          set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
          set sum [expr {$sum + $dz*$dz/$norm2}]
        }
        incr nh $nres
      }

      if { [$hz num] != 0 } {
        set hzx [vecsub $cpx [$hz get x]]
        set hzy [vecsub $cpy [$hz get y]]
        set hzz [vecsub $cpz [$hz get z]]
        foreach dx $hzx dy $hzy dz $hzz {
          set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
          set sum [expr {$sum + $dz*$dz/$norm2}]
        }
        incr nh $nres
      }

    }
    set arr($i) [expr {-1.5*$sum/$nh + 0.5}]
   puts [expr {-1.5*$sum/$nh + 0.5}]
  }
}

# c2 tail
proc orderparam-c2 { result { seltext all } } {
  upvar $result arr
  set n [molinfo top get numframes]
  for { set i 2 } { $i <= 20 } { incr i } {
#    puts $i
    set cp [atomselect top "($seltext) and lipid and name C2$i"]
    if { [$cp num] == 0 } {
#      puts "skipping $i"
    }
    set hx [atomselect top "($seltext) and lipid and name H${i}R"]
    set hy [atomselect top "($seltext) and lipid and name H${i}S"]
    set hz [atomselect top "($seltext) and lipid and name H${i}T"]
    set h9 [atomselect top "($seltext) and lipid and name H${i}1"]
    set sum 0.0
    set nh 0
    set nres [$cp num]
    for { set frame 500 } { $frame < $n } { incr frame } {
      $cp frame $frame
      $hx frame $frame
      $hy frame $frame
    $hz frame $frame
    $h9 frame $frame
      set cpx [$cp get x]
      set cpy [$cp get y]
      set cpz [$cp get z]

      if { [$hx num] != 0 } {
        set hxx [vecsub $cpx [$hx get x]]
        set hxy [vecsub $cpy [$hx get y]]
        set hxz [vecsub $cpz [$hx get z]]
        foreach dx $hxx dy $hxy dz $hxz {
          set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
          set sum [expr {$sum + $dz*$dz/$norm2}]
        }
        incr nh $nres
      }

      if { [$hy num] != 0 } {
        set hyx [vecsub $cpx [$hy get x]]
        set hyy [vecsub $cpy [$hy get y]]
        set hyz [vecsub $cpz [$hy get z]]
        foreach dx $hyx dy $hyy dz $hyz {
          set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
          set sum [expr {$sum + $dz*$dz/$norm2}]
        }
        incr nh $nres
      }

      if { [$hz num] != 0 } {
        set hzx [vecsub $cpx [$hz get x]]
        set hzy [vecsub $cpy [$hz get y]]
        set hzz [vecsub $cpz [$hz get z]]
        foreach dx $hzx dy $hzy dz $hzz {
          set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
          set sum [expr {$sum + $dz*$dz/$norm2}]
        }
        incr nh $nres
      }


      if { ($i == 9) || ($i == 10) } {	
      if { [$h9 num] != 0 } {
        set h9x [vecsub $cpx [$h9 get x]]
        set h9y [vecsub $cpy [$h9 get y]]
        set h9z [vecsub $cpz [$h9 get z]]
        foreach dx $h9x dy $h9y dz $h9z {
          set norm2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
          set sum [expr {$sum + $dz*$dz/$norm2}]
        }
        incr nh $nres
      }
      }



    }
    set arr($i) [expr {-1.5*$sum/$nh + 0.5}]
    puts [expr {-1.5*$sum/$nh + 0.5}]
  }
}

orderparam-c3 arr3
orderparam-c2 arr2

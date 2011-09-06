for f in height.gfs height1.gfs height2.gfs height3.gfs height4.gfs; do
    if gerris2D $f > ref.gfs; then :
    else
	echo "  FAIL: gerris2D $f"
	exit 1
    fi
    
    if mpirun -np 2 gerris2D $f > run.gfs; then :
    else
	echo "  FAIL: mpirun -np 2 gerris2D $f"
	exit 1
    fi
    
    for v in T_Hbx T_Hby T_Htx T_Hty K; do
	if gfscompare2D -v run.gfs ref.gfs $v 2> log; then :
	else
	    cat log
	    echo "  FAIL: $v $f"
	    exit 1
	fi
	if awk '{ if ($1 == "total" && $8 > 1e-10) exit 1; }' < log; then :
	else
	    cat log
	    echo "  FAIL: $v $f"
	    exit 1
	fi
    done
done

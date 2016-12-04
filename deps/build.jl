try
    run(`gfortran -fPIC -O2 -ggdb -shared -o daskr.so ../src/DASKR/ddaskr.f ../src/DASKR/dlinpk.f ../src/DASKR/daux.f`) 
end

try
    lib = @static is_apple() ? "dylib" : "so"
    run(`gfortran -fPIC -O2 -ggdb -shared -o daskr.$lib ../src/DASKR/ddaskr.f ../src/DASKR/dlinpk.f ../src/DASKR/daux.f`)
catch
end

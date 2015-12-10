import pyterpol
fitter = pyterpol.Fitter()
fitter.choose_fitter('nlopt_nelder_mead', xtol=1e-5, ftol=1e-5)
fitter.save('fitter_save.txt')
print fitter
fitter.clear_all()
print fitter
fitter.load('fitter_save.txt')
print fitter
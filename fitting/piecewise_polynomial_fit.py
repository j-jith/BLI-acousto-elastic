import numpy as np

def curve_fit(func, degree, f_0, f_1, n_steps, pieces):

    w_g = np.linspace(f_0, f_1, n_steps) # radians
    piece_len = (f_1 - f_0)/pieces

    index_pieces = [0]
    f_pieces = [f_0]
    while True:
        if f_pieces[-1] >= f_1:
            break
        else:
            index_pieces.append(np.argmin(np.abs(w_g - (f_pieces[-1]+piece_len))))
            f_pieces.append(w_g[index_pieces[-1]])

    ## Fitting loop
    cf_fit = []
    ww_fit = []
    for n_ii, ii in enumerate(index_pieces[:-1]):

        if n_ii < len(index_pieces)-2:
            ii_1 = index_pieces[n_ii+1] - 1
        else:
            ii_1 = index_pieces[n_ii+1]

        w = w_g[ii:ii_1]
        func_i = func[ii:ii_1]

        ## Curve fit

        #std = np.std(s.imag)
        #mean = np.mean(s.imag)
        #weight = 1/sqrt(2*pi*std**2)*np.exp(-(s.imag - mean)**2/2/std**2)
        #weight = weight/np.max(weight)

        #weight = [1 if x>150*2*pi and x<750*2*pi else 0 for x in s.imag]
        weight = [1 for x in w]

        cf_fit.append(np.polyfit(w, func_i, degree, w=weight))
        ww_fit.append(w)
    
    return ww_fit, cf_fit



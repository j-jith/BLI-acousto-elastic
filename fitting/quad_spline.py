def quad_spline(t, y, n, x=np.array([])):

    step = int(len(t)/n)
    
    if step*n != len(t):
        t = np.append(t[::step], t[-1])
        y = np.append(y[::step], y[-1])
    else:
        t = t[::step]
        y = y[::step]

    # derivatives
    z = np.array([0]) #arbitrary initial value
    for k in range(1, len(t)):
        z = np.append(z, -z[-1] + 2*(y[k] - y[k-1])/(t[k] - t[k-1]))

    alpha = (z[1] - z[0])/2/(t[1] - t[0])
    coeffs = np.array([[alpha, (z[0] - 2*alpha*t[0]), (alpha*t[0]**2 - z[0]*t[0] + y[0])]])
    for k in range(1, len(t)-1):
        alpha = (z[k+1] - z[k])/2/(t[k+1] - t[k])
        coeffs = np.append(coeffs, [[alpha, (z[k] - 2*alpha*t[k]), (alpha*t[k]**2 - z[k]*t[k] + y[k])]], axis=0)


    if len(x) != 0:
        y_new = np.array([])
        for xx in x:
            ii = np.argmin(np.abs(t-xx))
            if t[ii] > xx:
                ii = ii-1
            if ii > len(t)-2:
                ii = ii-1

            y_new = np.append(y_new, coeffs[ii][0]*xx**2 + coeffs[ii][1]*xx + coeffs[ii][2])

    return coeffs, y_new

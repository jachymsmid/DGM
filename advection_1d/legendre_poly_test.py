import matplotlib.pyplot as plt

####################################
#               EXACT              #
####################################

# exact Legendre polynomial evaluation
# up to order 10
def exact(n, x):
    if n > 10:
        print("n must be less or equal to 10")
        return
    match n:
        case 0:
            return 1.0
        case 1:
            return x
        case 2:
            return 1/2*(3*x**2 - 1)
        case 3:
            return 1/2*(5*x**3 - 3*x)
        case 4:
            return 1/8*(35*x**4 - 30*x**2 + 3)
        case 5:
            return 1/8*(63*x**5 - 70*x**3 + 15*x)
        case 6:
            return 1/16*(231*x**6 - 315*x**4 + 105*x**2 - 5)
        case 7:
            return 1/16*(429*x**7 - 693*x**5 + 315*x**3 - 35*x)
        case 8:
            return 1/128*(6435*x**8 - 12012*x**6 + 6930*x**4 - 1260*x**2 + 35)
        case 9:
            return 1/128*(12155*x**9 - 25740*x**7 + 18018*x**5 - 4620*x**3 + 315*x)
        case 10:
            return 1/256*(46189*x**10 - 109395*x**8 + 90090*x**6 - 30030*x**4 + 3465*x**2 - 63)


# exact evaluation of the Legendre polynomial
def exact_deriv(n, x):
    if n > 10:
        print("n must be less or equal to 10")
        return
    match n:
        case 0:
            return 0.0
        case 1:
            return 1.0
        case 2:
            return 1/2*(6*x)
        case 3:
            return 1/2*(15*x**2 - 3)
        case 4:
            return 1/8*(140*x**3 - 60*x)
        case 5:
            return 1/8*(315*x**4 - 210*x**2 + 15)
        case 6:
            return 1/16*(1386*x**5 - 1260*x**3 + 210*x)
        case 7:
            return 1/16*(7*429*x**6 - 5*693*x**4 + 3*315*x**2 - 35)
        case 8:
            return 1/128*(8*6435*x**7 - 6*12012*x**5 + 4*6930*x**3 - 2*1260*x)
        case 9:
            return 1/128*(9*12155*x**8 - 7*25740*x**6 + 5*18018*x**4 - 3*4620*x**2 + 315)
        case 10:
            return 1/256*(10*46189*x**9 - 8*109395*x**7 + 6*90090*x**5 - 4*30030*x**3 + 2*3465*x)


def normalized_exact(n, x):
    return sqrt((2*n+1)/2)*exact(n,x)


def normalized_exact_deriv(n, x):
    return sqrt((2*n+1)/2)*exact_deriv(n,x)


####################################
#            RECCURENT             #
####################################

# bonnet's formula, from wikipedia
def bonnet(n,x):
    if n == 0:
        return 1.0
    if n == 1:
        return x

    p0 = 1.0
    p1 = x

    for i in range(1,n):
        p2 = ((2*i+1)*x*p1-i*p0)/(i+1)
        p0 = p1
        p1 = p2
    return p2


# bonnet's formula from the .jpynb
def bonnet2(n,x):
    if n == 0:
        return 1.0
    if n == 1:
        return x

    p0 = 1
    p1 = x

    for i in range(2,n+1):
        p2 = ((2*i-1)*x*p1-(i-1)*p0)/i
        p0 = p1
        p1 = p2
    return p2


# hesthavens formula
def hesthaven(n, x):
    pass

def spectral(n, x):
    if n == 0:
        return 1.0
    if n == 1:
        return x
    p0 = 1.0
    p1 = x
    for i in range(1, n):
        a_n = 2*(i+1)**2 * 2*i
        c_n = 2*i*(2*i+1)*(2*i+2)
        d_n = 2*i**2 * (2*i + 2)
        p2 = (c_n*x*p1-d_n*p0)/a_n

        p0 = p1
        p1 = p2
    return p2


def main():
    x = 1.0
    eps = 1e-15
    max_err = 0.0

    print(f"\nTesting Legendre polynomial evaluation at x = {x}")


    print("- Bonnet v1")
    for i in range(11):
        max_err = max(exact(i,x) - bonnet(i,x), max_err)
    if max_err < eps:
        print("\tOK")
    else:
        print(f"\tERROR: maximal err = {max_err}")

    max_err = 0.0
    print("- Bonnet v2")
    for i in range(11):
        max_err = max(exact(i,x) - bonnet(i,x), max_err)
    if max_err < eps:
        print("\tOK")
    else:
        print(f"\tERROR: maximal err = {max_err}")

    max_err = 0.0
    print("- Spectral")
    for i in range(11):
        max_err = max(exact(i,x) - spectral(i,x), max_err)
    if max_err < eps:
        print("\tOK")
    else:
        print(f"\tERROR: maximal err = {max_err}")

    # TODO:
    print("- Hesthaven")
    for i in range(11):
        max_err = max(exact(i,x) - spectral(i,x), max_err)
    if max_err < eps:
        print("\tOK")
    else:
        print(f"\tERROR: maximal err = {max_err}")

    print(f"\nTesting the evaluation of the derivative of the Legendre polynomial at x = {x}")

    max_err = 0.0
    print("- Wiki")
    for i in range(11):
        max_err = max(exact(i,x) - bonnet(i,x), max_err)
    if max_err < eps:
        print("\tOK")
    else:
        print(f"\tERROR: maximal err = {max_err}")

    max_err = 0.0
    print("- jpynb")
    for i in range(11):
        max_err = max(exact(i,x) - bonnet(i,x), max_err)
    if max_err < eps:
        print("\tOK")
    else:
        print(f"\tERROR: maximal err = {max_err}")

if __name__ == "__main__":
    main()

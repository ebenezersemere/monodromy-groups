def Gamma1(n):
    """
    Gamma1 returns the lattes map for a given n (int) using the elliptic curve y^2 = x^3 + 1.
    """
    x, y, z = var('x, y, z') # initialize variables x, y, z
    E = EllipticCurve(y^2 - x^3 - 1).base_extend(QuadraticField(-3)) # initialize elliptical curve y^2 = x^3 + 1.
    a = E.division_polynomial(2*n, two_torsion_multiplicity=1) # numerator
    b = E.division_polynomial(n, two_torsion_multiplicity=1) # denominator
    c = ((1/2)*((1)-((a)/(2*(b**4))))) # equation
    d = c.substitute(x=(4*z*(z-1))**(1/3), y=(1-2*z)).factor() # complete
    return d

def NGamma1(n):
    """
    NGamma1 calls Gamma1 to print all lattes map with n up to the given n (int).
    """
    for val in range(2, n + 1): # loops through values
        print(f" Lattes map: {Gamma1(val)} \n\t") # prints lattes map for specific n

def DegSeq1(n):

    # First define the Division Polynomials for the elliptic curve y^2 = x^3 + 1.
    E = EllipticCurve([0,1])

    z = var('z'); x=(4*z*(z-1))**(1/3); y=(1-2*z)
    a = E.division_polynomial(2*n, (x,y), two_torsion_multiplicity=1) # numerator
    b = E.division_polynomial(n, (x,y), two_torsion_multiplicity=1) # denominator
    c = ((1/2)*((1)-((a)/(2*(b**4))))).factor()

    # Compute the Belyi maps by turning d into a string and converting back
    R.<z> = FunctionField(QQbar)
    gamma = sage_eval( str(c), locals={'z':z})

    # Compute those divisors which give information about the Black Vertices, White Vertices, and Faces
    DivList = [ gamma.divisor_of_zeros(), (1-gamma).divisor_of_zeros(), (1/gamma).divisor_of_zeros() ]

    # Pull out information about the degrees of these vertices
    DegreeSequence = [sorted([ D.multiplicity(p) for p in D.support() ]) for D in DivList ]

    return DegreeSequence

def NDegSeq1(n):
    """
    NDegSeq1 calls DegSeq1 to print all degree sequences whose n is contained in first, last (ints).
    """
    for val in range(2, n + 1): # loops through values
        print(f" Degree sequence for n={val}: {DegSeq1(val)} \n\t") # prints degree sequence for specific n

def Combined1(n):
    for val in range(2, n + 1):
        print(f" N = {val} \n Lattes map: {Gamma1(val)} \n Degree sequence: {DegSeq1(val)} \n\t")
        # prints lattes map & degree sequence for specific n

###########
###########
###########

def Gamma2(n):
    """
    Gamma2 returns the lattes map for a given n (int) using the elliptic curve y^2 = x^3 - x.
    """
    x, y, z = var('x, y, z') # initialize variables x, y, z
    E = EllipticCurve(y^2 - x^3 + x).base_extend(QuadraticField(-3)) # initialize elliptical curve y^2 = x^3 - x.
    a = E.division_polynomial(n+1, two_torsion_multiplicity=1) # numerator
    b = E.division_polynomial(n-1, two_torsion_multiplicity=1) # numerator
    c = E.division_polynomial(n, two_torsion_multiplicity=1) # denominator
    d = (x - ((a*b)/(c**2))) ** 2 # equation
    e = d.substitute(x=(z**(1/2)), y=(z**(1/4))*((z-1)**(1/2))).factor() # complete
    return e

def NGamma2(n):
    """
    NGamma2 calls Gamma2 to print all lattes map with n up to the given n (int).
    """
    for val in range(2, n + 1): # loops through values
        print(f" Lattes map for n={val}: {Gamma2(val)} \n\t") # prints lattes map for specific n

def DegSeq2(n):

    # First define the Division Polynomials for the elliptic curve y^2 = x^3 - x.
    E = EllipticCurve([-1,0])

    z = var('z'); x = (z^(1/2)); y = (z^(1/4))*((z-1)^(1/2))
    a = E.division_polynomial(n+1, (x,y), two_torsion_multiplicity=1)
    b = E.division_polynomial(n-1, (x,y), two_torsion_multiplicity=1)
    c = E.division_polynomial(n, (x,y), two_torsion_multiplicity=1)
    d = ((x - (a*b/c^2))^2).factor()

    # Compute the Belyi maps by turning d into a string and converting back
    R.<z> = FunctionField(QQbar)
    gamma = sage_eval( str(d), locals={'z':z})

    # Compute those divisors which give information about the Black Vertices, White Vertices, and Faces
    DivList = [ gamma.divisor_of_zeros(), (1-gamma).divisor_of_zeros(), (1/gamma).divisor_of_zeros() ]

    # Pull out information about the degrees of these vertices
    DegreeSequence = [sorted([ D.multiplicity(p) for p in D.support() ]) for D in DivList ]

    return DegreeSequence

def NDegSeq2(n):
    """
    NDegSeq2 calls DegSeq2 to print all degree sequences whose n is contained in first, last (ints).
    """
    for val in range(2, n + 1): # loops through values
        print(f" Degree sequence for n={val}: {DegSeq3(val)} \n\t") # prints degree sequence for specific n

def Combined2(n):
    for val in range(2, n + 1):
        print(f" N = {val} \n Lattes map for n={val}: {Gamma2(val)} \n Degree sequence for n={val}: {DegSeq2(val)} \n\t")
        # prints lattes map & degree sequence for specific n

###########
###########
###########

def Gamma3(n):
    """
    Gamma3 returns the lattes map for a given n (int) using the elliptic curve y^2 = x^3 + 1.
    """
    x, y, z = var('x, y, z') # initialize variables x, y, z
    E = EllipticCurve(y^2 - x^3 - 1).base_extend(QuadraticField(-3)) # initialize elliptic curve y^2 = x^3 + 1.
    a = E.division_polynomial(n-1, two_torsion_multiplicity=1) # numerator
    b = E.division_polynomial(n+1, two_torsion_multiplicity=1) # numerator
    c = E.division_polynomial(n, two_torsion_multiplicity=1) # denominator
    d = - (( x - ((a * b) / (c ** 2))) ** 3) # equation
    e = d.substitute(x=(-z**(1/3)), y=((1-z)**(1/2))).factor() # complete
    return e

def NGamma3(n):
    """
    NGamma3 calls Gamma3 to print all lattes map with n up to the given n (int).
    """
    for val in range(2, n + 1): # loops through values
        print(f" Lattes map for n={val}: {Gamma3(val)} \n\t") # prints lattes map for specific n

def DegSeq3(n):

    # First define the Division Polynomials for the elliptic curve y^2 = x^3 + 1.
    E = EllipticCurve([0,1])
    z = var('z'); x =(-z**(1/3)); y=((1-z)**(1/2))
    a = E.division_polynomial(n-1, (x,y), two_torsion_multiplicity=1) # numerator
    b = E.division_polynomial(n+1, (x,y), two_torsion_multiplicity=1) # numerator
    c = E.division_polynomial(n, (x,y), two_torsion_multiplicity=1) # denominator
    d = (- (( x - ((a * b) / c ** 2)) ** 3)).factor() # equation.factor()

    # Compute the Belyi maps by turning d into a string and converting back
    R.<z> = FunctionField(QQbar)
    gamma = sage_eval( str(d), locals={'z':z})

    # Compute those divisors which give information about the Black Vertices, White Vertices, and Faces
    DivList = [ gamma.divisor_of_zeros(), (1-gamma).divisor_of_zeros(), (1/gamma).divisor_of_zeros() ]

    # Pull out information about the degrees of these vertices
    DegreeSequence = [sorted([ D.multiplicity(p) for p in D.support() ]) for D in DivList ]

    return DegreeSequence

def NDegSeq3(n):
    """
    NDegSeq3 calls DegSeq3 to print all degree sequences whose n is contained in first, last (ints).
    """
    for val in range(2, n + 1): # loops through values
        print(f" Degree sequence for n={val}: {DegSeq3(val)} \n\t") # prints degree sequence for specific n

def Combined3(n):
    for val in range(2, n + 1):
        print(f" N = {val} \n Lattes map for n={val}: {Gamma3(val)} \n Degree sequence for n={val}: {DegSeq3(val)} \n\t")
        # prints lattes map & degree sequence for specific n

###########
###########
###########

def Gamma_ID(a, b, c):
    """
    Group_ID takes three generators, sigma 0, 1, and infinity (tuple of lists), and returns the group ID.
    """
    G = PermutationGroup([a, b, c])
    return G.group_id()


def NGamma_ID(group, n):
    """
    NGroup_ID calls Group_ID to return a specified group ID given a group and n (ints).
    """
    if int(group) == 1:
        if int(n) == 1:
            return Gamma_ID([1], [1], [1])
        elif int(n) == 2:
            return Gamma_ID([(1, 3, 2)], [(2, 3, 4)], [(1, 2, 4)]) # M = 3, N = 2
        elif int(n) == 3:
            return Gamma_ID([(1, 6, 5), (2, 4, 3), (7, 9, 8)], [(1, 2, 3), (4, 7, 8), (5, 6, 9)], [(1, 4, 9), (2, 5, 7)]) # Slide 19, N = 3
        elif int(n) == 4:
            return Gamma_ID([(1, 8, 7), (2, 9, 5), (3, 6, 10), (11, 13, 12), (14, 16, 15)], [(1, 2, 3), (4, 5, 6), (7, 8, 16), (9, 14, 11), (10, 12, 15)], [(1, 10, 16), (2, 7, 14), (3, 5, 4), (6, 9, 12), (11, 15, 13)]) # Slide 19, N = 4

    elif int(group) == 2:
        if int(n) == 1:
            return Gamma_ID([1], [1], [1])
        elif int(n) == 2:
            return Gamma_ID([(1, 2, 4, 3)], [(1, 2), (3, 4)], [(2, 3)]) # M = 4, N = 2
        elif int(n) == 3:
            return Gamma_ID([(1, 2, 9, 8), (4, 5, 7, 6)], [(1, 5), (2, 6), (3, 4), (8, 9)], [(1, 4, 3, 6), (2, 7, 5, 8)]) # Slide 25, N = 3
        elif int(n) == 4:
            return Gamma_ID([(1, 2, 16, 15), (3, 5, 6, 4), (7, 9, 13, 11), (8, 12, 14, 10)], [(1, 9), (2, 10), (3, 7), (4, 8), (5, 6), (11, 12), (13, 14), (15, 16)], [(12, 13), (1, 7, 4, 10), (2, 14, 9, 15), (3, 11, 8, 6)]) # Slide 25, N = 4

    elif int(group) == 3:
        if int(n) == 1:
            return Gamma_ID([1], [1], [1])
        elif int(n) == 2:
            return Gamma_ID([(1, 4, 3)], [(1, 2), (3, 4)], [(1, 2, 3)]) # M = 6, N = 2
        elif int(n) == 3:
            return Gamma_ID([(1, 9, 8), (2, 3, 4), (5, 7, 6)], [(1, 2), (3, 5), (4, 6), (8, 9)], [(3, 6), (1, 4, 7, 5, 2, 8)]) # Slide 31, N = 3
        elif int(n) == 4:
            return Gamma_ID([(1, 16, 15), (2, 3, 4), (6, 9, 10), (7, 13, 11), (8, 12, 14)], [(1, 2), (3, 7), (4, 8), (5, 6), (9, 11), (10, 12), (13, 14), (15, 16)], [(9, 13, 12), (1, 4, 14, 7, 2, 15), (3, 11, 6, 5, 10, 8)]) # Slide 31, N = 4

def Beta_ID(group, n):
    """
    Match_ID takes a group and n (ints) and returns the ID for the match.
    """
    F.<a,b,c> = FreeGroup()
    if int(group) == 1:
        G1 = F / [a^n, b^n, c^3, a*b*a^(n-1)*b^(n-1), (c*a*c^2)*b^(n-1), (c*b*c^2)*(a*b)]
        H1 = G1.as_permutation_group()
        return H1.group_id()
    elif int(group) == 2:
        G2 = F / [a^n, b^n, c^4, a*b*a^(n-1)*b^(n-1), (c*a*c^3)*b^(n-1), (c*b*c^3)*a]
        H2 = G2.as_permutation_group()
        return H2.group_id()
    elif int(group) == 3:
        G3 = F / [a^n, b^n, c^6, a*b*a^(n-1)*b^(n-1), c*a*c^(5)*b^(n-1), c*b*c^(5)*a*b^(n-1)]
        H3 = G3.as_permutation_group()
        return H3.group_id()

def Label_Check():
    """
    NLabel_Check calls Label_Check to display various pairs of corresponding group IDs.
    """
    for group in range(1, 4):
        for n in range(1, 5):
            print(f"ID Pair for Groups {group} & N = {n}: Gamma {NGamma_ID(group, n)} & Beta {Beta_ID(group, n)} \n")

###########
###########
###########


def main():

if __name__ == '__main__':
    main()

for n in range(2, 100):
    print("""  
    auto y%d = update<mpq_class, %d, m>(my_y0 * (w * y%d + %s));
    mpq_class e%d = integrate(my_y0 * (w * y%d + %s)) / integrate(my_y0 * my_y0);
    std::cout << e%d << '\\n';
    """ % (n, n, n-1,
    '+'.join('mpq_class(-e%d) * y%d' % (i, n - i) for i in range(1, n)) + '+ mpq_class(-e%d) * my_y0' % n, 
    n+1, n, 
    '+'.join('mpq_class(-e%d) * y%d' % (i, n + 1 - i) for i in range(1, n+1)), 
    n+1,
    ))

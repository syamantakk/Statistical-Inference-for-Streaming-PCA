function eta = get_learning_rate(n, alpha, eigengap)
    eta = alpha*log(n)/(n * eigengap);
end
#outputFinal[ipga,] = c(p1, p2, p3, p12, p13, p23, p123)
# yy = c(p1, p2, p3, p12, p13, p23, p123)
p1 = 0.653381
p2 = 0.5
p3 = 0.414935
p12 = 0.797074
p13 = 0.743130
p23 = 0.620224
p123 = 0.825666

yy = c(p1, p2, p3, p12, p13, p23, p123)

##solve ccf element
A <- c(1,0,0,1,1,0,1,
       0,1,0,1,0,1,1,
       0,0,1,0,1,1,1,
       1,1,0,1,1,1,1,
       1,0,1,1,1,1,1,
       0,1,1,1,1,1,1,
       1,1,1,1,1,1,1)
A <- matrix(A, nrow = 7, byrow = TRUE)

B <- log(1-yy)

sol <- solve(A, B)

X <- 1-exp(sol)

Q1 <- X[1]
Q2 <- X[2]
Q3 <- X[3]
Q12 <- X[4]
Q13 <- X[5]
Q23 <- X[6]
Q123 <- X[7]

parallel_S = Q1*Q2*Q3 + Q123 + Q1*Q23 + Q2*Q13 + Q3*Q12 + Q12*Q13 + Q13*Q23 + Q12*Q23 
serial_S = 1-(1-Q1)*(1-Q2)*(1-Q3)*(1-Q12)*(1-Q13)*(1-Q23)*(1-Q123)

print(parallel_S)
print(serial_S)
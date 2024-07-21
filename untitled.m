

Path1 = readmatrix('OscProbPath.csv')
Path2 = readmatrix('ModelPath.csv')
Path2bkp =  readmatrix('ModelPath.csv')

test(1:15,1) = Path2(1:15,2)
test(16,1) = Path2(16,2) + Path2(17,2)
test(17:85,1) = Path2(18:86,2)
